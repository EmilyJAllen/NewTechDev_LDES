%% ldes_main.m
% Annual charge/discharge model of a thermal (carbon block + TPV) LDES serving a data center.
%
% Pipeline
%   1. ldes_calibrate   runs the 3D block model (thermal_block_run) and builds a fast lumped
%                       surrogate of one block: TPV power, coolant heat and insulation loss vs.
%                       stored energy. Done once, cached in ldes_surrogate.mat.
%   2. ldes_cycle       the charge/discharge efficiency model: kWh in -> electricity, coolant
%                       heat, insulation loss, energy left, for a chosen discharge window.
%   3. ldes_dispatch    a year of 15-minute steps: price-threshold charge/discharge plus
%                       backup during grid outages, repeated over random outage years.
%
% Outputs: annual cost of charging electricity, unserved load, CO2 from charging electricity,
%          plus capex, round trip, losses and avoided grid purchases.
%
% SITE ASSUMPTION: the data center is in western Pennsylvania (PJM market, eGRID subregion
% RFCW). Price, emissions and outage inputs below are sourced for that location.
%
% SOURCES (inputs below cite these by tag)
%   [LOAD]  Load_Profile_1000_MWitcapacity_1.2_PUE.csv: 1 week (Jan 1-7 2026) of 15-min data,
%           1000 MW IT capacity, PUE 1.2; total_MW = IT (incl. distribution/UPS) + non-IT.
%   [SOM]   Monitoring Analytics (PJM IMM), 2025 State of the Market Report for PJM, Vol. 1:
%           2025 real-time load-weighted average LMP $50.73/MWh, median $38.11, std $58.74.
%   [eGRID] US EPA eGRID2023 Summary Tables (rev2, 2025): RFCW total output CO2 rate 911.4 lb/MWh.
%   [EIA]   EIA Electric Power Annual, Table 11.3 (2025 data), Pennsylvania, all events incl.
%           major event days: SAIFI 1.408 interruptions/yr, SAIDI 458.1 min/yr.
%   [ANT]   Antora Energy (antora.com/system-launch): energy stored above 1800 C at full scale,
%           design above 1500 C; Contrary Research profile: >40 % heat-to-electric TPV,
%           <1 %/day heat loss (company claims).
%   [TEA]   Kelsall, Buznitsky & Henry, "Technoeconomic Analysis of Thermal Energy Grid Storage
%           Using Graphite and Tin", arXiv:2106.07624v2 (2021): system cost ~$50/kWh-e at
%           >1 GWh scale (energy-only CPE <$20/kWh-e), lifetime >30 yr.
%   [OMB]   OMB M-25-23 (2025) reinstated Circular A-94 (1992): 7 % real discount rate.

clear; clc; close all;

%% ============================== INPUTS ==============================
% --- Load profile [LOAD] ----------------------------------------------------
in.load.file       = 'Load_Profile_1000_MWitcapacity_1.2_PUE.csv';
in.load.column     = 'total_MW';  % whole-facility load (IT x PUE): what the grid / battery must serve
in.load.units      = 'MW';
in.load.scaleTo_MW = [];          % use the CSV as is (mean ~907 MW)
in.load.scaleBasis = 'mean';
in.load.startDow   = 3;           % Jan 1 2026 is a Thursday (Mon = 0); the week is tiled to 2026
in.dt_h            = 0.25;        % 15-minute steps, as in the CSV

% --- Outages [EIA] -----------------------------------------------------------
in.outage.freq_per_yr = 1.408;           % PA SAIFI 2025, all events
in.outage.dur_h       = 458.1/1.408/60;  % PA SAIDI/SAIFI 2025 = 5.42 h average outage (CAIDI)
in.outage.mode        = 'random';        % Poisson count, uniform start times (see ldes_outages)
in.outage.startHours  = [];
in.outage.nMC         = 100;             % random outage years averaged
in.outage.seed        = 1;

% --- Electricity price & dispatch rule [SOM] ------------------------------------
in.price.mode       = 'pjm2025';    % calibrated to PJM 2025 statistics (see ldes_price.m);
                                    % 'csv' + a real PJM LMP export replaces it
in.price.mean_MWh   = 50.73;        % load-weighted average real-time LMP, 2025
in.price.median_MWh = 38.11;
in.price.std_MWh    = 58.74;
in.price.seed       = 2025;
in.price.file       = '';  in.price.column = [];  in.price.units = '$/MWh';   % for 'csv'
in.Pmax             = in.price.median_MWh/1000;   % $/kWh: price > Pmax discharge, < Pmax charge.
                                    % The 2025 median splits the year evenly between the two.
in.Pcharge_below    = [];           % [] = single threshold, as specified
in.reserveFrac      = 0.0;          % 0 = the stated rule (no energy held back for outages)
in.allowExport      = false;        % price-driven discharge only offsets the data-center load

% --- Emissions [eGRID] ---------------------------------------------------------
in.EF_t_per_kWh = 911.4 * 0.45359237e-3 / 1000;  % 911.4 lb/MWh -> 4.134e-4 t CO2/kWh (RFCW)

% --- Battery --------------------------------------------------------------------
in.batt.capacity_kWh = [];      % [] = one block exactly as modelled in thermal_block_3D_full.m
                                %      (2.44 x 6.1 x 2.6 m graphite, 1000 -> 2000 K = 35.3 MWh_th)
in.batt.Pcharge_kW   = [];      % [] = the 3D model's heater rating, 1e5 W/m3 x 38.7 m3 = 3.87 MW
in.batt.etaHeater    = 1.0;     % Joule heating deposits all electrical input in the block;
                                % transformer/rectifier losses are outside the model boundary
in.batt.E0frac       = 1.0;
in.batt.warmStart    = true;    % start the year from the end state of an outage-free year
in.batt.cost_per_kWh = 50;      % $/kWh-e, total system [TEA] (GWh-scale; a single block costs more)
in.batt.costBasis    = 'electric';  % kWh-e = electricity from one full discharge (Tmax -> empty)
in.batt.life_yr      = 30;      % [TEA] (>30 yr); also NREL ATB's default 30-yr cost recovery
in.batt.discount     = 0.07;    % real [OMB]

% --- Surrogate / calibration ------------------------------------------------------
in.sur.file        = 'ldes_surrogate.mat';
in.sur.recalibrate = false;     % true after you change thermal_block_params.m
in.sur.opts        = struct('Tmax', 2000, 'Tmin', 1000, 'tDis_h', 60);
                                % Tmax 2000 K = 1727 C: inside Antora's >1500 C design range and
                                % below the >1800 C demonstrated at full scale [ANT].
                                % Tmin 1000 K: where the 3D model's TPV efficiency curve reaches 0.

% --- Efficiency-model demo (ldes_cycle) ---------------------------------------------
in.cycle.Echarge_kWh  = [];     % [] = 80% of capacity
in.cycle.tDischarge_h = [4 8 12 24];
in.cycle.tStandby_h   = 0;
in.cycle.validate3D   = false;  % true: also run the full 3D model for the first window (slow)

in.makePlots = true;
%% ====================================================================

%% 1) Surrogate of one block
if in.sur.recalibrate || exist(in.sur.file, 'file') ~= 2
    o = in.sur.opts;  o.file = in.sur.file;
    sur = ldes_calibrate(thermal_block_params(), o);
else
    S = load(in.sur.file);  sur = S.sur;
    fprintf('Loaded surrogate %s (mesh %dx%dx%d, %s)\n', in.sur.file, sur.mesh, sur.created);
end
if isempty(in.batt.capacity_kWh), in.batt.capacity_kWh = sur.Ecap_kWh; end
scale = in.batt.capacity_kWh / sur.Ecap_kWh;
if isempty(in.batt.Pcharge_kW), in.batt.Pcharge_kW = sur.Pcharge_kW*scale; end
if abs(scale - 1) > 1e-6
    fprintf(['Battery = %.2f block-equivalents: power, losses and capacity scale linearly.\n' ...
             '  (Conservative for losses: one bigger block has less surface per kWh.)\n'], scale);
end

%% 2) Charge/discharge efficiency model
Ech = in.cycle.Echarge_kWh;
if isempty(Ech), Ech = 0.8*in.batt.capacity_kWh; end
co = struct('sur', sur, 'scale', scale, 'Pcharge_kW', in.batt.Pcharge_kW, ...
            'etaHeater', in.batt.etaHeater, 'tStandby_h', in.cycle.tStandby_h);
fprintf('\nEfficiency model: charge %.0f kWh, standby %g h\n', Ech, in.cycle.tStandby_h);
fprintf('%10s %12s %12s %12s %12s %12s %8s\n', 'window(h)', 'elec(kWh)', 'coolant', 'insul.', 'heater', 'left', 'RT');
cyc = cell(numel(in.cycle.tDischarge_h), 1);
for j = 1:numel(in.cycle.tDischarge_h)
    c = ldes_cycle(Ech, in.cycle.tDischarge_h(j), co);  cyc{j} = c;
    fprintf('%10g %12.0f %12.0f %12.0f %12.0f %12.0f %8.3f\n', c.tDischarge_h, c.Eel_kWh, ...
            c.Ecool_kWh, c.Eins_kWh, c.Eheater_kWh, c.Eleft_kWh, c.roundTrip);
end
if in.cycle.validate3D
    c3 = ldes_cycle(Ech, in.cycle.tDischarge_h(1), setfield(co, 'mode', '3d')); %#ok<SFLD>
    fprintf('%10s %12.0f %12.0f %12.0f %12.0f %12.0f %8.3f   <- full 3D check\n', '3D', c3.Eel_kWh, ...
            c3.Ecool_kWh, c3.Eins_kWh, c3.Eheater_kWh, c3.Eleft_kWh, c3.roundTrip);
end
% Electric capacity (kWh-e): electricity from one full block discharge (Tmax uniform -> cold),
% integrated from the 3D calibration run; used for the $/kWh-e capex basis.
Ee_full_kWh = scale*trapz(sur.raw.A.t_h, sur.raw.A.P_W(:,2))/1e3;

%% 3) Annual inputs
[load_kW, linfo] = ldes_load_profile(struct('file', in.load.file, 'column', in.load.column, ...
    'units', in.load.units, 'scaleTo_MW', in.load.scaleTo_MW, 'scaleBasis', in.load.scaleBasis, 'dt_h', in.dt_h));
n = numel(load_kW);
in.price.startDow = in.load.startDow;
[price, pinfo] = ldes_price(in.price, n, in.dt_h, load_kW);
fprintf('\nLoad: %s [%s] | mean %.1f MW, peak %.1f MW, %.0f MWh/yr\n', linfo.source, linfo.column, linfo.mean_MW, linfo.peak_MW, linfo.annual_MWh);
if isfield(pinfo, 'mean_MWh')
    fprintf('Price (%s): load-wtd mean $%.2f, median $%.2f, std $%.2f, max $%.0f /MWh; Pmax $%.2f/MWh\n', ...
            pinfo.mode, pinfo.mean_MWh, pinfo.median_MWh, pinfo.std_MWh, pinfo.max_MWh, 1000*in.Pmax);
end

cfg.dt_h = in.dt_h;  cfg.scale = scale;  cfg.Pcharge_kW = in.batt.Pcharge_kW;
cfg.etaHeater = in.batt.etaHeater;  cfg.priceDischarge = in.Pmax;
cfg.priceCharge = in.Pmax;  if ~isempty(in.Pcharge_below), cfg.priceCharge = in.Pcharge_below; end
cfg.reserveFrac = in.reserveFrac;  cfg.EF_t_per_kWh = in.EF_t_per_kWh;  cfg.allowExport = in.allowExport;
cfg.E0_kWh = in.batt.E0frac*in.batt.capacity_kWh;

if in.batt.warmStart                                  % periodic start: end-of-year state, no outages
    w = ldes_dispatch(load_kW, price, false(n,1), sur, cfg);
    cfg.E0_kWh = w.E_kWh(end);
end

%% 4) Outage years (Monte Carlo)
if exist('rng', 'file') || exist('rng', 'builtin'), rng(in.outage.seed); else, rand('seed', in.outage.seed); end %#ok<RAND>
nMC = in.outage.nMC;  if ~strcmpi(in.outage.mode, 'random'), nMC = 1; end
f = {'chargeCost_USD','Echarge_kWh','CO2_t','unserved_kWh','outageLoad_kWh','Eel_kWh', ...
     'Eel_backup_kWh','Eel_arbitrage_kWh','avoidedPurchase_USD','Ecool_kWh','Eins_kWh', ...
     'Eheater_kWh','roundTrip','CO2_net_t','balance_kWh','unservedHours'};
R = zeros(nMC, numel(f));  nOut = zeros(nMC,1);  runs = cell(nMC,1);
for m = 1:nMC
    [mask, ev] = ldes_outages(n, in.dt_h, in.outage.freq_per_yr, in.outage.dur_h, in.outage.mode, in.outage.startHours);
    r = ldes_dispatch(load_kW, price, mask, sur, cfg);
    r.mask = mask;  r.events = ev;  runs{m} = r;  nOut(m) = size(ev,1);
    for q = 1:numel(f), R(m,q) = r.(f{q}); end
end
mu = mean(R, 1);  col = @(name) R(:, strcmp(f, name));
pct    = @(v, q) subsref(sort(v(:)), struct('type','()','subs',{{1 + round(q/100*(numel(v)-1))}}));
fmt    = @(x) regexprep(sprintf('%.0f', x), '\d(?=(\d{3})+$)', '$0,');
plural = @(k) repmat('s', 1, k ~= 1);

%% 5) Report
Cap = in.batt.capacity_kWh;
if strcmpi(in.batt.costBasis, 'electric'), capex = in.batt.cost_per_kWh*Ee_full_kWh;
else,                                      capex = in.batt.cost_per_kWh*Cap; end
crf = in.batt.discount*(1+in.batt.discount)^in.batt.life_yr / ((1+in.batt.discount)^in.batt.life_yr - 1);

fprintf('\n=================== ANNUAL RESULTS (%d outage year%s) ===================\n', nMC, plural(nMC));
fprintf('Battery: %.0f kWh_th (%.2f blocks), charger %.0f kW, peak TPV %.0f kW, self-discharge %.2f %%/day at full\n', ...
        Cap, scale, in.batt.Pcharge_kW, max(sur.Pel_kW)*scale, sur.selfDischarge_pct_per_day_full);
fprintf('Outages: %.2f/yr x %.1f h (mean %.2f simulated/yr)\n', in.outage.freq_per_yr, in.outage.dur_h, mean(nOut));
fprintf('\n--- Requested outputs ---\n');
fprintf('Annual cost of charging electricity:   $%12s   (P5-P95 $%s - $%s)\n', fmt(mu(strcmp(f,'chargeCost_USD'))), ...
        fmt(pct(col('chargeCost_USD'),5)), fmt(pct(col('chargeCost_USD'),95)));
fprintf('Unserved load:                          %12s kWh/yr (P95 %s kWh, worst %s kWh)\n', fmt(mu(strcmp(f,'unserved_kWh'))), ...
        fmt(pct(col('unserved_kWh'),95)), fmt(max(col('unserved_kWh'))));
fprintf('   = %.1f %% of load during outages (%s kWh/yr demanded during outages)\n', ...
        100*mu(strcmp(f,'unserved_kWh'))/max(mu(strcmp(f,'outageLoad_kWh')),eps), fmt(mu(strcmp(f,'outageLoad_kWh'))));
fprintf('Annual CO2 from charging electricity:   %12.1f t\n', mu(strcmp(f,'CO2_t')));
fprintf('\n--- Supporting numbers ---\n');
fprintf('Charging electricity drawn:            %12s kWh\n', fmt(mu(strcmp(f,'Echarge_kWh'))));
fprintf('TPV electricity delivered:             %12s kWh  (backup %s, price-driven %s)\n', fmt(mu(strcmp(f,'Eel_kWh'))), ...
        fmt(mu(strcmp(f,'Eel_backup_kWh'))), fmt(mu(strcmp(f,'Eel_arbitrage_kWh'))));
fprintf('Heat to TPV coolant:                   %12s kWh\n', fmt(mu(strcmp(f,'Ecool_kWh'))));
fprintf('Heat lost via insulation:              %12s kWh\n', fmt(mu(strcmp(f,'Eins_kWh'))));
fprintf('Heater loss:                           %12s kWh\n', fmt(mu(strcmp(f,'Eheater_kWh'))));
fprintf('Annual round trip (elec out / in):     %12.3f\n', mu(strcmp(f,'roundTrip')));
fprintf('Grid purchases avoided by TPV output:  $%12s\n', fmt(mu(strcmp(f,'avoidedPurchase_USD'))));
fprintf('Net energy-arbitrage cost:             $%12s   (charge cost - avoided purchases)\n', ...
        fmt(mu(strcmp(f,'chargeCost_USD')) - mu(strcmp(f,'avoidedPurchase_USD'))));
fprintf('Net CO2 (charging - displaced grid):   %12.1f t\n', mu(strcmp(f,'CO2_net_t')));
fprintf('Capex: $%s (%s basis: %s kWh-e per full discharge), annualized $%s/yr (%g yr, %.0f%% real)\n', fmt(capex), in.batt.costBasis, ...
        fmt(Ee_full_kWh), fmt(capex*crf), in.batt.life_yr, 100*in.batt.discount);
fprintf('Energy-balance residual (should be ~0): %.3g kWh\n\n', max(abs(col('balance_kWh'))));

%% 6) Plots (the realization closest to the mean unserved load)
if in.makePlots
    [~, mi] = min(abs(col('unserved_kWh') - mu(strcmp(f,'unserved_kWh'))));
    r = runs{mi};  t_d = (0:n-1)'*in.dt_h/24;
    figure('Name','Annual state of charge');
    subplot(2,1,1);
    plot(t_d, 100*r.E_kWh/Cap, 'LineWidth',1); hold on;
    yl = [min(0, min(100*r.E_kWh/Cap)) 105];
    for j = 1:size(r.events,1)
        patch(t_d(r.events(j,[1 2 2 1])), yl([1 1 2 2]), [1 0.6 0.6], 'EdgeColor','none', 'FaceAlpha',0.4);
    end
    hold off; ylim(yl); grid on; ylabel('State of charge (%)'); title('Battery state of charge (red = outage)');
    subplot(2,1,2);
    plot(t_d, price, 'LineWidth',1); hold on; plot(t_d([1 end]), [1 1]*in.Pmax, 'k--'); hold off;
    grid on; xlabel('Day of year'); ylabel('Price ($/kWh)'); legend('Price','P_{max}','Location','best');

    if ~isempty(r.events)
        k0 = max(r.events(1,1) - round(48/in.dt_h), 1);  k1 = min(r.events(1,2) + round(48/in.dt_h), n);
        kk = k0:k1;  th = (kk - r.events(1,1))*in.dt_h;
        figure('Name','First outage');
        subplot(2,1,1);
        plot(th, load_kW(kk)/1e3, th, r.Pel_kW(kk)/1e3, th, r.Punserved_kW(kk)/1e3, th, r.Pgrid_kW(kk)/1e3, 'LineWidth',1.3);
        grid on; ylabel('MW'); legend('Load','TPV electric','Unserved','Grid -> charger','Location','best');
        title('Around the first outage (t = 0 at outage start)');
        subplot(2,1,2);
        plot(th, 100*r.E_kWh(kk)/Cap, 'LineWidth',1.3); grid on;
        xlabel('Hours from outage start'); ylabel('State of charge (%)');
    end

    if nMC > 1
        figure('Name','Unserved load across outage years');
        hist(col('unserved_kWh')/1e3, max(10, round(nMC/5)));
        xlabel('Unserved load (MWh/yr)'); ylabel('Count'); grid on;
        title(sprintf('%d random outage years', nMC));
    end
end
