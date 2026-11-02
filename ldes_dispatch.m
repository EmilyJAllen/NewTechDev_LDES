function out = ldes_dispatch(load_kW, price, outage, sur, cfg)
% LDES_DISPATCH  Annual 15-minute charge/discharge simulation of the thermal LDES.
%
% Each step, in priority order:
%   1) Grid outage      -> battery discharges to serve the load (TPV power permitting).
%                          Whatever it cannot cover is UNSERVED LOAD. No charging.
%   2) price > priceDischarge and energy above the reserve
%                       -> discharge to offset grid purchases (up to the load; no export
%                          unless cfg.allowExport), never dipping below the reserve.
%   3) price < priceCharge and not full
%                       -> charge from the grid at up to the charger rating.
%   4) otherwise        -> idle, shutters closed (standby self-discharge only).
% Discharge power is modulated with the shutter opening x in [0,1]; the maximum
% available TPV power depends on the block's state of charge (from the surrogate).
%
% cfg fields (all required; ldes_main sets them):
%   dt_h, scale, Pcharge_kW, etaHeater, priceDischarge, priceCharge,
%   reserveFrac, E0_kWh, EF_t_per_kWh, allowExport

n   = numel(load_kW);  dt = cfg.dt_h;  sc = cfg.scale;
Ecap = sur.Ecap_kWh*sc;   Efl = sur.Efloor_kWh*sc;   Eres = cfg.reserveFrac*Ecap;
eta  = cfg.etaHeater;     Pchmax = cfg.Pcharge_kW;

% local copies of the surrogate grid (per block) for a fast inline lookup
gE0 = sur.E_kWh(1);  gdE = sur.dE;  gn = numel(sur.E_kWh);
tEl = sur.Pel_kW;  tCo = sur.Pcool_kW;  tLo = sur.PlossOpen_kW;  tLc = sur.PlossClosed_kW;

E = cfg.E0_kWh;
Es = zeros(n,1);  Pgrid = zeros(n,1);  Pel = zeros(n,1);  Pco = zeros(n,1);
Plo = zeros(n,1); Pun = zeros(n,1);    mode = zeros(n,1,'int8');

for k = 1:n
    % --- surrogate lookup at the start-of-step energy
    s = (E/sc - gE0)/gdE;
    i = floor(s);  if i < 0, i = 0; elseif i > gn-2, i = gn-2; end
    w = s - i;     if w < 0, w = 0; elseif w > 1, w = 1; end
    i = i + 1;
    pel = sc*(tEl(i) + w*(tEl(i+1)-tEl(i)));
    pco = sc*(tCo(i) + w*(tCo(i+1)-tCo(i)));
    plo = sc*(tLo(i) + w*(tLo(i+1)-tLo(i)));
    plc = sc*(tLc(i) + w*(tLc(i+1)-tLc(i)));

    x = 0;  pin = 0;
    if outage(k)                                           % 1) backup
        if pel > 0, x = min(1, load_kW(k)/pel); end
        Pun(k) = load_kW(k) - x*pel;
        mode(k) = 2;
    elseif price(k) > cfg.priceDischarge && E > Eres && pel > 0   % 2) price-driven discharge
        if cfg.allowExport, x = 1; else, x = min(1, load_kW(k)/pel); end
        dOut = pel + pco + plo - plc;                      % extra drain per unit x
        if dOut > 0
            xmax = ((E - Eres)/dt - plc)/dOut;             % stay above the reserve
            x = max(0, min(x, xmax));
        end
        mode(k) = -1;
    elseif price(k) < cfg.priceCharge && E < Ecap          % 3) charge
        pin = min(Pchmax, max(0, ((Ecap - E)/dt + plc)/eta));
        mode(k) = 1;
    end

    ploss = x*plo + (1-x)*plc;
    Pgrid(k) = pin;  Pel(k) = x*pel;  Pco(k) = x*pco;  Plo(k) = ploss;
    E = E + (eta*pin - x*(pel + pco) - ploss)*dt;
    if E < Efl, E = Efl; end
    Es(k) = E;
end

kWh = @(v) sum(v)*dt;
out.E_kWh   = Es;  out.Pgrid_kW = Pgrid;  out.Pel_kW = Pel;  out.Pcool_kW = Pco;
out.Ploss_kW = Plo; out.Punserved_kW = Pun; out.mode = mode;  out.Ecap_kWh = Ecap;

out.Echarge_kWh   = kWh(Pgrid);                            % grid electricity into the charger
out.chargeCost_USD = sum(Pgrid.*price)*dt;
out.CO2_t         = cfg.EF_t_per_kWh * out.Echarge_kWh;
out.unserved_kWh  = kWh(Pun);
out.outageLoad_kWh = kWh(load_kW(outage));
out.unservedHours = nnz(Pun > 1e-6)*dt;
out.Eel_kWh       = kWh(Pel);
out.Eel_backup_kWh = kWh(Pel(outage));
out.Eel_arbitrage_kWh = out.Eel_kWh - out.Eel_backup_kWh;
out.avoidedPurchase_USD = sum(Pel(~outage).*price(~outage))*dt;   % grid energy the TPV output displaced
out.Ecool_kWh     = kWh(Pco);
out.Eins_kWh      = kWh(Plo);
out.Eheater_kWh   = (1-eta)*out.Echarge_kWh;
out.dE_kWh        = Es(end) - cfg.E0_kWh;
out.balance_kWh   = out.Echarge_kWh - out.Eheater_kWh - out.Eel_kWh - out.Ecool_kWh - out.Eins_kWh - out.dE_kWh;
out.roundTrip     = out.Eel_kWh / max(out.Echarge_kWh, eps);
out.nCycles_equiv = out.Eel_kWh / max(Ecap, eps);
out.CO2_net_t     = cfg.EF_t_per_kWh * (out.Echarge_kWh - out.Eel_arbitrage_kWh);
end
