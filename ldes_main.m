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
%                       backup during grid outages, repeated over random outage years for
%                       each outage-duration scenario (median / bad / catastrophic), for two
%                       dispatch cases: Arbitrage (price rule) and Backup-only (held full).
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
%           major event days: SAIFI 1.408 interruptions/yr (used for outage frequency).
%   [ANT]   Antora Energy (antora.com/system-launch): energy stored above 1800 C at full scale,
%           design above 1500 C; Contrary Research profile: >40 % heat-to-electric TPV,
%           <1 %/day heat loss (company claims).
%   [TEA]   Kelsall, Buznitsky & Henry, "Technoeconomic Analysis of Thermal Energy Grid Storage
%           Using Graphite and Tin", arXiv:2106.07624v2 (2021): system cost ~$50/kWh-e at
%           >1 GWh scale (energy-only CPE <$20/kWh-e), lifetime >30 yr.
%   [OMB]   OMB M-25-23 (2025) reinstated Circular A-94 (1992): 7 % real discount rate.
%   [TEAM]  Outage duration scenarios (median / bad / catastrophic), provided by the project team.

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

% --- Outages: three duration scenarios ------------------------------------------
% Each scenario is run as its own set of nMC random outage years. The outage COUNT and
% START TIMES are drawn identically across scenarios (same seed), so the scenarios differ
% only in how long each outage lasts.
in.outage.freq_per_yr = 1.408;           % [EIA] PA SAIFI 2025, all events (outages / yr)
in.outage.scenarios   = struct( ...
    'name',        {'Median', 'Bad',  'Catastrophic'}, ...
    'dur_h',       { 0.473,    0.605,  70.3         }, ...   % [TEAM] hours per outage
    'freq_per_yr', { [],       [],     []           });      % [] = in.outage.freq_per_yr
in.outage.mode        = 'random';        % Poisson count, uniform start times (see ldes_outages)
in.outage.startHours  = [];              % for 'given' mode
in.outage.nMC         = 100;             % random outage years per scenario
in.outage.seed        = 1;
in.outage.resultsCSV  = 'ldes_results_by_case.csv';   % '' = don't write

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
% Dispatch cases, each run against every outage scenario. reserveFrac = share of capacity
% that price-driven discharge may not touch (kept for outages).
in.cases = struct( ...
    'name',        {'Arbitrage', 'Backup-only'}, ...
    'reserveFrac', { 0,           1           });   % 0: the stated price rule as is
                                                   % 1: never discharge for price; blocks are
                                                   %    topped up when price < Pmax and held full
in.allowExport      = false;        % price-driven discharge only offsets the data-center load

% --- Emissions [eGRID] ---------------------------------------------------------
in.EF_t_per_kWh = 911.4 * 0.45359237e-3 / 1000;  % 911.4 lb/MWh -> 4.134e-4 t CO2/kWh (RFCW)

% --- Battery: sized to the data center ---------------------------------------------
% The battery is N identical blocks, each exactly as modelled in thermal_block_3D_full.m
% (2.44 x 6.1 x 2.6 m graphite, 1000 -> 2000 K = 35.3 MWh_th, 3.87 MW heater).
% N = the fewest blocks that carry the data center's PEAK load through one full outage of
% the sizing scenario's duration with zero unserved load, starting fully charged
% (ldes_size_blocks). TPV power drops steeply as blocks drain, so this is power-limited.
in.batt.sizing.mode     = 'outage';  % 'outage' (size as above) | 'fixed' (use nBlocks)
in.batt.sizing.scenario = 'Median';  % which in.outage.scenarios entry to size for
in.batt.sizing.loadBasis = 'peak';   % 'peak' or 'mean' facility load
in.batt.sizing.startSOC = 1.0;       % SOC assumed at outage start when sizing
in.batt.nBlocks         = [];        % used when mode = 'fixed'
in.batt.Pcharge_cap_kW  = 'peak';    % total charging power limit: 'peak' = facility peak load
                                     % (grid interconnection sized for the site), a number in kW,
                                     % or [] for no cap (each block at its own 3.87 MW heater)
in.batt.etaHeater    = 1.0;     % Joule heating deposits all electrical input in the block;
                                % transformer/rectifier losses are outside the model boundary
in.batt.E0frac       = 1.0;
in.batt.warmStart    = true;    % start the year from the end state of an outage-free year
in.batt.cost_per_kWh = 50;      % $/kWh-e, total system [TEA] (GWh scale, which this battery is)
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
in.plotDir   = 'plots';        % PNGs saved here ('' = don't save)
%% ====================================================================

%% 1) Surrogate of one block
if in.sur.recalibrate || exist(in.sur.file, 'file') ~= 2
    o = in.sur.opts;  o.file = in.sur.file;
    sur = ldes_calibrate(thermal_block_params(), o);
else
    S = load(in.sur.file);  sur = S.sur;
    fprintf('Loaded surrogate %s (mesh %dx%dx%d, %s)\n', in.sur.file, sur.mesh, sur.created);
end

%% 1b) Load profile and battery sizing
[load_kW, linfo] = ldes_load_profile(struct('file', in.load.file, 'column', in.load.column, ...
    'units', in.load.units, 'scaleTo_MW', in.load.scaleTo_MW, 'scaleBasis', in.load.scaleBasis, 'dt_h', in.dt_h));
n = numel(load_kW);
fprintf('Load: %s [%s] | mean %.1f MW, peak %.1f MW, %.0f MWh/yr\n', linfo.source, linfo.column, ...
        linfo.mean_MW, linfo.peak_MW, linfo.annual_MWh);

switch lower(in.batt.sizing.mode)
    case 'outage'
        k = find(strcmpi({in.outage.scenarios.name}, in.batt.sizing.scenario), 1);
        if isempty(k), error('Sizing scenario "%s" is not in in.outage.scenarios', in.batt.sizing.scenario); end
        if strcmpi(in.batt.sizing.loadBasis, 'mean'), Lsz = mean(load_kW); else, Lsz = max(load_kW); end
        [nBlk, szi] = ldes_size_blocks(sur, Lsz, in.outage.scenarios(k).dur_h, in.batt.sizing.startSOC);
        fprintf(['Sizing: %d blocks carry %.1f MW (%s load) for %.3f h (%s outage) from %.0f%% SOC\n' ...
                 '        (power-only bound %d blocks, energy-only bound %d blocks)\n'], nBlk, Lsz/1e3, ...
                in.batt.sizing.loadBasis, szi.D_h, in.outage.scenarios(k).name, 100*szi.startSOC, szi.N_power, szi.N_energy);
    case 'fixed'
        nBlk = in.batt.nBlocks;
    otherwise
        error('Unknown sizing mode "%s"', in.batt.sizing.mode);
end
scale = nBlk;
if strcmpi(in.batt.sizing.mode,'outage'), sizeNote = sprintf('%s outage (%.3f h) at %s load from %.0f%% SOC', in.outage.scenarios(k).name, szi.D_h, in.batt.sizing.loadBasis, 100*szi.startSOC);
else, sizeNote = 'fixed block count'; end
in.batt.capacity_kWh = nBlk*sur.Ecap_kWh;
Pheaters = nBlk*sur.Pcharge_kW;
if ischar(in.batt.Pcharge_cap_kW) && strcmpi(in.batt.Pcharge_cap_kW, 'peak'), cap = max(load_kW);
elseif isempty(in.batt.Pcharge_cap_kW), cap = Inf;
else, cap = in.batt.Pcharge_cap_kW; end
in.batt.Pcharge_kW = min(Pheaters, cap);
fprintf('Battery: %d blocks = %.2f GWh_th; heaters %.0f MW, charging limited to %.0f MW\n\n', ...
        nBlk, in.batt.capacity_kWh/1e6, Pheaters/1e3, in.batt.Pcharge_kW/1e3);

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

%% 3) Prices
in.price.startDow = in.load.startDow;
[price, pinfo] = ldes_price(in.price, n, in.dt_h, load_kW);
if isfield(pinfo, 'mean_MWh')
    fprintf('Price (%s): load-wtd mean $%.2f, median $%.2f, std $%.2f, max $%.0f /MWh; Pmax $%.2f/MWh\n', ...
            pinfo.mode, pinfo.mean_MWh, pinfo.median_MWh, pinfo.std_MWh, pinfo.max_MWh, 1000*in.Pmax);
end

cfg.dt_h = in.dt_h;  cfg.scale = scale;  cfg.Pcharge_kW = in.batt.Pcharge_kW;
cfg.etaHeater = in.batt.etaHeater;  cfg.priceDischarge = in.Pmax;
cfg.priceCharge = in.Pmax;  if ~isempty(in.Pcharge_below), cfg.priceCharge = in.Pcharge_below; end
cfg.EF_t_per_kWh = in.EF_t_per_kWh;  cfg.allowExport = in.allowExport;

%% 4) Dispatch cases x outage scenarios (Monte Carlo per scenario)
f = {'chargeCost_USD','Echarge_kWh','CO2_t','unserved_kWh','outageLoad_kWh','outageHours', ...
     'Eel_kWh','Eel_backup_kWh','Eel_arbitrage_kWh','avoidedPurchase_USD','Ecool_kWh','Eins_kWh', ...
     'Eheater_kWh','roundTrip','CO2_net_t','balance_kWh','unservedHours'};
sc_ = in.outage.scenarios;  nS = numel(sc_);  nC = numel(in.cases);
nMC = in.outage.nMC;  if ~strcmpi(in.outage.mode, 'random'), nMC = 1; end
res = struct('case', {}, 'reserveFrac', {}, 'name', {}, 'dur_h', {}, 'freq', {}, 'R', {}, 'nOut', {}, 'E0', {}, 'events', {});
for c = 1:nC
    cfg.reserveFrac = in.cases(c).reserveFrac;
    cfg.E0_kWh = in.batt.E0frac*in.batt.capacity_kWh;
    if in.batt.warmStart                              % periodic start: end-of-year state, no outages
        w = ldes_dispatch(load_kW, price, false(n,1), sur, cfg);
        cfg.E0_kWh = w.E_kWh(end);
    end
    for s = 1:nS
        fr = sc_(s).freq_per_yr;  if isempty(fr), fr = in.outage.freq_per_yr; end
        rng(in.outage.seed);                          % common random numbers across scenarios and cases
        R = zeros(nMC, numel(f));  nOut = zeros(nMC,1);  evs = cell(nMC,1);
        for m = 1:nMC
            [ofrac, ev] = ldes_outages(n, in.dt_h, fr, sc_(s).dur_h, in.outage.mode, in.outage.startHours);
            r = ldes_dispatch(load_kW, price, ofrac, sur, cfg);
            nOut(m) = size(ev,1);  evs{m} = ev;     % outage times, so any year can be replotted
            for q = 1:numel(f), R(m,q) = r.(f{q}); end
        end
        res(end+1) = struct('case', in.cases(c).name, 'reserveFrac', cfg.reserveFrac, 'name', sc_(s).name, ...
                            'dur_h', sc_(s).dur_h, 'freq', fr, 'R', R, 'nOut', nOut, 'E0', cfg.E0_kWh, 'events', {evs}); %#ok<SAGROW>
        fprintf('%-12s | %-13s %6.3f h x %.3f/yr: done (%d outage years)\n', in.cases(c).name, ...
                sc_(s).name, sc_(s).dur_h, fr, nMC);
    end
end
col    = @(R, name) R(:, strcmp(f, name));
avg    = @(R, name) mean(col(R, name));
pct    = @(v, q) subsref(sort(v(:)), struct('type','()','subs',{{1 + round(q/100*(numel(v)-1))}}));
fmt    = @(x) regexprep(sprintf('%.0f', round(x) + 0), '\d(?=(\d{3})+$)', '$0,');

%% 5) Report
Cap = in.batt.capacity_kWh;
if strcmpi(in.batt.costBasis, 'electric'), capex = in.batt.cost_per_kWh*Ee_full_kWh;
else,                                      capex = in.batt.cost_per_kWh*Cap; end
crf = in.batt.discount*(1+in.batt.discount)^in.batt.life_yr / ((1+in.batt.discount)^in.batt.life_yr - 1);

fprintf('\n================== ANNUAL RESULTS BY DISPATCH CASE AND OUTAGE SCENARIO ==================\n');
fprintf('Battery: %d blocks, %.2f GWh_th, charging %.0f MW, peak TPV %.0f MW, self-discharge %.2f %%/day at full\n', ...
        nBlk, Cap/1e6, in.batt.Pcharge_kW/1e3, max(sur.Pel_kW)*scale/1e3, sur.selfDischarge_pct_per_day_full);
fprintf('Sized for: %s\n', sizeNote);
fprintf('Each column: mean of %d random outage years (P95 / worst where shown)\n\n', nMC);

rows = { ...  % label, value function, format ('fmt' = integer with thousands separators)
  'Reserve for outages (% of capacity)',    @(x) 100*x.reserveFrac,                '%.0f'
  'Outage duration (h/outage)',             @(x) x.dur_h,                          '%.3f'
  'Outage frequency (outages/yr)',          @(x) x.freq,                           '%.3f'
  'Outages simulated (mean/yr)',            @(x) mean(x.nOut),                     '%.2f'
  'Grid-out hours (h/yr)',                  @(x) avg(x.R,'outageHours'),           '%.2f'
  '--- Requested outputs ---',              [],                                    ''
  'Annual cost of charging elec. ($)',      @(x) avg(x.R,'chargeCost_USD'),        'fmt'
  'Unserved load (kWh/yr)',                 @(x) avg(x.R,'unserved_kWh'),          'fmt'
  '   P95 (kWh/yr)',                        @(x) pct(col(x.R,'unserved_kWh'),95),  'fmt'
  '   worst year (kWh/yr)',                 @(x) max(col(x.R,'unserved_kWh')),     'fmt'
  '   % of load during outages',            @(x) 100*avg(x.R,'unserved_kWh')/max(avg(x.R,'outageLoad_kWh'),eps), '%.2f'
  'Annual CO2 from charging elec. (t)',     @(x) avg(x.R,'CO2_t'),                 '%.1f'
  '--- Supporting ---',                     [],                                    ''
  'Load during outages (kWh/yr)',           @(x) avg(x.R,'outageLoad_kWh'),        'fmt'
  'Charging electricity (kWh/yr)',          @(x) avg(x.R,'Echarge_kWh'),           'fmt'
  'TPV electricity, backup (kWh/yr)',       @(x) avg(x.R,'Eel_backup_kWh'),        'fmt'
  'TPV electricity, price-driven (kWh/yr)', @(x) avg(x.R,'Eel_arbitrage_kWh'),     'fmt'
  'Heat to TPV coolant (kWh/yr)',           @(x) avg(x.R,'Ecool_kWh'),             'fmt'
  'Heat lost via insulation (kWh/yr)',      @(x) avg(x.R,'Eins_kWh'),              'fmt'
  'Annual round trip (elec out / in)',      @(x) avg(x.R,'roundTrip'),             '%.3f'
  'Grid purchases avoided ($/yr)',          @(x) avg(x.R,'avoidedPurchase_USD'),   'fmt'
  'Net energy-arbitrage cost ($/yr)',       @(x) avg(x.R,'chargeCost_USD') - avg(x.R,'avoidedPurchase_USD'), 'fmt'
  'Net CO2, charging - displaced (t)',      @(x) avg(x.R,'CO2_net_t'),             '%.1f'
  'Energy-balance residual, max (kWh)',     @(x) max(abs(col(x.R,'balance_kWh'))), '%.1e' };

cw = 14;  nR = numel(res);
fprintf('%-40s', '');  for j = 1:nR, fprintf('%*s', cw, res(j).case); end;  fprintf('\n');
fprintf('%-40s', '');  for j = 1:nR, fprintf('%*s', cw, res(j).name); end;  fprintf('\n');
tbl = nan(size(rows,1), nR);
for i = 1:size(rows,1)
    fprintf('%-40s', rows{i,1});
    if isempty(rows{i,2}), fprintf('\n'); continue; end
    for s = 1:nR
        v = rows{i,2}(res(s));  tbl(i,s) = v;
        if strcmp(rows{i,3}, 'fmt'), str = fmt(v); else, str = sprintf(rows{i,3}, v); end
        fprintf('%*s', cw, str);
    end
    fprintf('\n');
end
fprintf('\nCapex: $%s (%s basis: %s kWh-e per full discharge), annualized $%s/yr (%g yr, %.0f%% real)\n\n', ...
        fmt(capex), in.batt.costBasis, fmt(Ee_full_kWh), fmt(capex*crf), in.batt.life_yr, 100*in.batt.discount);

if ~isempty(in.outage.resultsCSV)                     % same table, machine-readable
    fid = fopen(in.outage.resultsCSV, 'w');
    fprintf(fid, 'metric');  for s = 1:nR, fprintf(fid, ',%s: %s', res(s).case, res(s).name); end;  fprintf(fid, '\n');
    for i = 1:size(rows,1)
        if isempty(rows{i,2}), continue; end
        fprintf(fid, '"%s"', strtrim(rows{i,1}));  fprintf(fid, ',%.10g', tbl(i,:));  fprintf(fid, '\n');
    end
    fprintf(fid, '"Capex ($)"');                fprintf(fid, ',%.10g', capex*ones(1,nR));      fprintf(fid, '\n');
    fprintf(fid, '"Annualized capex ($/yr)"');  fprintf(fid, ',%.10g', capex*crf*ones(1,nR));  fprintf(fid, '\n');
    fclose(fid);
    fprintf('Results table written to %s\n\n', in.outage.resultsCSV);
end

%% 6) Plots: one figure per dispatch case (a row per outage scenario) + a summary figure
% For each case x scenario the plotted year is the outage year whose unserved load is
% closest to that scenario's mean (among years with at least one outage). It is re-run
% from its stored outage start times, so it matches the Monte Carlo exactly.
if in.makePlots
    if ~isempty(in.plotDir) && ~exist(in.plotDir, 'dir'), mkdir(in.plotDir); end
    t_d   = (0:n-1)'*in.dt_h/24;
    clrs  = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
    for c = 1:nC
        cfg.reserveFrac = in.cases(c).reserveFrac;
        fig = figure('Name', ['Dispatch case: ' in.cases(c).name], 'Position', [50 50 1500 300*nS]);
        for s = 1:nS
            x  = res((c-1)*nS + s);
            uv = col(x.R, 'unserved_kWh');  cand = find(x.nOut > 0);
            if isempty(cand), continue; end
            [~, j] = min(abs(uv(cand) - mean(uv)));  m = cand(j);
            ev = x.events{m};
            ofrac = ldes_outages(n, in.dt_h, 0, x.dur_h, 'given', ev(:,3));
            cfg.E0_kWh = x.E0;
            r = ldes_dispatch(load_kW, price, ofrac, sur, cfg);
            soc = 100*r.E_kWh/Cap;

            % (a) the year: state of charge, outages shaded
            subplot(nS, 3, 3*(s-1) + 1);
            yl = [min(0, min(soc)) 135];                       % strip above 100% holds the legend
            for e = 1:size(ev,1)
                hp = patch(ev(e,[3 4 4 3])/24 + [-0.5 0.5 0.5 -0.5].*(x.dur_h < 24), [yl(1) yl(1) 100 100], ...
                           [1 0.75 0.75], 'EdgeColor', [0.9 0.3 0.3]);  hold on;
                if e == 1, hpo = hp; end
            end
            hl = plot(t_d, soc, 'Color', clrs(1,:), 'LineWidth', 0.8);  hold off;
            xlim([0 t_d(end)]);  ylim(yl);  set(gca, 'YTick', 0:20:100);  grid on;  ylabel('State of charge (%)');
            legend([hpo hl], {'Grid outage', 'State of charge'}, 'Location', 'north', ...
                   'Orientation', 'horizontal', 'FontSize', 7);
            if s == nS, xlabel('Day of year'); end
            title(sprintf('%s outages (%.3g h): year shown has %d outage%s, %s MWh unserved', ...
                  x.name, x.dur_h, size(ev,1), repmat('s', 1, size(ev,1) ~= 1), fmt(uv(m)/1e3)));

            % (b) zoom on the outage with the most unserved load: power flows
            [~, e1] = max(arrayfun(@(q) sum(r.Punserved_kW(ev(q,1):ev(q,2))), 1:size(ev,1)));
            pre = 6;  post = max(12, 0.3*x.dur_h);
            k0 = max(floor((ev(e1,3) - pre)/in.dt_h) + 1, 1);
            k1 = min(ceil((ev(e1,4) + post)/in.dt_h), n);
            kk = k0:k1;  th = (kk' - 1)*in.dt_h - ev(e1,3);
            subplot(nS, 3, 3*(s-1) + 2);
            ytop = 1.8*max(load_kW)/1e3;                        % headroom for the legend
            patch([0 x.dur_h x.dur_h 0], [0 0 1 1]*1.15*max(load_kW)/1e3, [1 0.75 0.75], 'EdgeColor', 'none');  hold on;
            stairs(th, load_kW(kk)/1e3,        'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
            stairs(th, r.Pel_kW(kk)/1e3,       'Color', clrs(3,:), 'LineWidth', 1.4);
            stairs(th, r.Punserved_kW(kk)/1e3, 'Color', clrs(2,:), 'LineWidth', 1.4);
            stairs(th, r.Pgrid_kW(kk)/1e3,     'Color', clrs(4,:), 'LineWidth', 1.0);  hold off;
            xlim([th(1) th(end)]);  ylim([0 ytop]);  grid on;  ylabel('MW (15-min average)');
            legend({'Outage','Load','TPV output','Unserved','Grid to chargers'}, 'Location', 'northeast', 'FontSize', 7);
            if s == nS, xlabel('Hours from outage start'); end
            title(sprintf('Outage with most unserved load: %s MWh unserved', ...
                  fmt(sum(r.Punserved_kW(ev(e1,1):ev(e1,2)))*in.dt_h/1e3)));

            % (c) same window: state of charge
            subplot(nS, 3, 3*(s-1) + 3);
            hp = patch([0 x.dur_h x.dur_h 0], [0 0 100 100], [1 0.75 0.75], 'EdgeColor', 'none');  hold on;
            hl = plot(th, soc(kk), 'Color', clrs(1,:), 'LineWidth', 1.4);  hold off;
            xlim([th(1) th(end)]);  ylim([0 135]);  set(gca, 'YTick', 0:20:100);  grid on;  ylabel('State of charge (%)');
            legend([hp hl], {'Grid outage', 'State of charge'}, 'Location', 'north', ...
                   'Orientation', 'horizontal', 'FontSize', 7);
            if s == nS, xlabel('Hours from outage start'); end
            title(sprintf('SOC %.0f%% at outage start, %.0f%% at end', ...
                  soc(max(ev(e1,1)-1,1)), soc(ev(e1,2))));
        end
        if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
            sgtitle(sprintf('%s dispatch (reserve %.0f%%): %d blocks, %.1f GWh_{th}', ...
                    in.cases(c).name, 100*in.cases(c).reserveFrac, nBlk, Cap/1e6));
        end
        if ~isempty(in.plotDir)
            print(fig, fullfile(in.plotDir, sprintf('case_%s.png', regexprep(lower(in.cases(c).name), '[^a-z0-9]+', '_'))), '-dpng', '-r130');
        end
    end

    % Summary: the three requested outputs, case vs case, per scenario
    fig = figure('Name', 'Requested outputs by case', 'Position', [50 50 1500 420]);
    M = @(name, k) reshape(arrayfun(@(x) avg(x.R, name), res), nS, nC)/k;
    pan = {'chargeCost_USD', 1e6, 'Annual charging cost ($M)',         'linear', '$%.1fM'
           'unserved_kWh',   1e3, 'Unserved load (MWh/yr, log scale)',  'log',    'num'
           'CO2_t',          1e3, 'Charging CO_2 (kt/yr)',             'linear', '%.0f kt'};
    for p = 1:3
        subplot(1, 3, p);
        v = M(pan{p,1}, pan{p,2});
        islog = strcmp(pan{p,4}, 'log');
        if islog                                       % bars rise from a floor below the smallest value
            base = 10^floor(log10(max(min(v(v > 0)), 1e-3)) - 1);
            vb = max(v, base);
            hb = bar(vb, 'BaseValue', base);
        else
            hb = bar(v);
        end
        for c = 1:nC, set(hb(c), 'FaceColor', clrs(c,:)); end
        set(gca, 'XTick', 1:nS, 'XTickLabel', {sc_.name}, 'YScale', pan{p,4});  grid on;  ylabel(pan{p,3});
        if islog, ylim([base, 10^ceil(log10(max(vb(:))) + 0.5)]); else, ylim([0 1.4*max(v(:))]); end
        for c = 1:nC                                    % value labels on top of each bar
            for k = 1:nS
                xb = k + (c - (nC + 1)/2)*0.8/nC;
                yv = v(k,c);  if islog, yv = max(yv, base); end
                if islog, yt = yv*1.6; else, yt = yv + 0.02*max(v(:)); end
                if strcmp(pan{p,5}, 'num')
                    if v(k,c) >= 100, lbl = fmt(v(k,c)); else, lbl = sprintf('%.2g', v(k,c)); end
                else
                    lbl = sprintf(pan{p,5}, v(k,c));
                end
                text(xb, yt, lbl, 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'bottom', 'FontSize', 7);
            end
        end
        if p == 1, legend({in.cases.name}, 'Location', 'northeast'); end
        title(pan{p,3});
    end
    if ~isempty(in.plotDir), print(fig, fullfile(in.plotDir, 'summary_requested_outputs.png'), '-dpng', '-r130'); end

    % Price series vs Pmax (shared by every case)
    fig = figure('Name', 'Electricity price', 'Position', [50 50 1200 320]);
    plot(t_d, 1000*price, 'Color', clrs(1,:), 'LineWidth', 0.5);  hold on;
    plot(t_d([1 end]), [1 1]*1000*in.Pmax, 'k--', 'LineWidth', 1.2);  hold off;
    xlim([0 t_d(end)]);  ylim([0 pct(1000*price, 99.9)]);  grid on;
    xlabel('Day of year');  ylabel('$/MWh');  legend({'Price', 'P_{max}'}, 'Location', 'northeast');
    title('Electricity price (PJM-2025-calibrated, y-axis cut at the 99.9th percentile)');
    if ~isempty(in.plotDir), print(fig, fullfile(in.plotDir, 'price.png'), '-dpng', '-r130'); end
    if ~isempty(in.plotDir), fprintf('Plots saved to %s%s\n', in.plotDir, filesep); end
end
