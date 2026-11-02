function out = ldes_cycle(Echarge_kWh, tDischarge_h, opts)
% LDES_CYCLE  Charge/discharge efficiency model for one charge -> (standby) -> discharge cycle.
%
%   out = ldes_cycle(Echarge_kWh, tDischarge_h)
%   out = ldes_cycle(Echarge_kWh, tDischarge_h, opts)
%
% Takes in : Echarge_kWh   electricity drawn from the grid to charge (kWh)
%            tDischarge_h  length of the discharge window, TPV faces fully open (h)
% Outputs  : out.Eel_kWh    electricity delivered by the TPV cells (kWh)
%            out.Ecool_kWh  heat rejected to the TPV coolant (kWh)
%            out.Eins_kWh   heat lost into the insulation over the whole cycle (kWh)
%            out.Eheater_kWh  heater/power-electronics loss, (1-etaHeater)*charge (kWh)
%            out.Eleft_kWh  energy left in the block above T_min at the end (kWh)
%            out.roundTrip  Eel / Echarge
%            out.balance_kWh  Echarge - (Eel+Ecool+Eins+Eheater+dE_block); ~0
%
% The block starts uniform at T_min (empty), charges at the heater rating until
% Echarge has been drawn (or the block is full), sits shuttered for tStandby_h,
% then discharges with all TPV faces open for tDischarge_h.
%
% opts.mode       'surrogate' (default, fast lumped model from ldes_calibrate)
%                 '3d'        (runs the full 3D model; use it to validate the surrogate)
% opts.sur        surrogate struct (default: loaded from opts.surFile)
% opts.surFile    'ldes_surrogate.mat'
% opts.scale      battery size in block-equivalents (default 1)
% opts.Pcharge_kW charger rating (default: the 3D model heater rating x scale)
% opts.etaHeater  grid -> block heat efficiency (default 1, i.e. block-only)
% opts.tStandby_h shuttered hold between charge and discharge (default 0)
% opts.dt_h       surrogate time step (default 1/60 h)
% opts.p, opts.mesh, opts.dtOut_h   3D-mode parameters, mesh override, output spacing

if nargin < 3, opts = struct(); end
opts = setdef(opts, 'mode', 'surrogate');
opts = setdef(opts, 'surFile', 'ldes_surrogate.mat');
if ~isfield(opts, 'sur') || isempty(opts.sur)
    S = load(opts.surFile);  opts.sur = S.sur;
end
sur  = opts.sur;
opts = setdef(opts, 'scale', 1);
opts = setdef(opts, 'Pcharge_kW', sur.Pcharge_kW*opts.scale);
opts = setdef(opts, 'etaHeater', 1);
opts = setdef(opts, 'tStandby_h', 0);
opts = setdef(opts, 'dt_h', 1/60);
opts = setdef(opts, 'dtOut_h', 5/60);
opts = setdef(opts, 'mesh', []);

tCh  = Echarge_kWh / opts.Pcharge_kW;          % h
tOp  = tCh + opts.tStandby_h;
tEnd = tOp + tDischarge_h;

switch lower(opts.mode)
    %% ------------------------------------------------------------------ surrogate
    case 'surrogate'
        sc = opts.scale;  Ecap = sur.Ecap_kWh*sc;  Efl = sur.Efloor_kWh*sc;
        dt = opts.dt_h;   n = max(ceil(tEnd/dt - 1e-9), 1);
        E  = 0;  Ein = 0;  Eel = 0;  Ecool = 0;  Eins = 0;  Eh = 0;
        ts = zeros(n+1,1);  Es = zeros(n+1,1);  Pels = zeros(n+1,1);
        for k = 1:n
            t0 = (k-1)*dt;  h = min(dt, tEnd - t0);  tm = t0 + h/2;
            [pel, pco, plo, plc] = ldes_lookup(sur, E, sc);
            Pin = 0;  x = 0;
            if tm < tCh                                  % charging
                Pin = min(opts.Pcharge_kW, (Echarge_kWh - Ein)/h);
                Pin = min(Pin, max(0, ((Ecap - E)/h + plc)/opts.etaHeater));   % block full
            elseif tm >= tOp                             % discharging, faces open
                x = 1;
            end
            Ploss = x*plo + (1-x)*plc;
            Ein   = Ein   + Pin*h;
            Eh    = Eh    + (1-opts.etaHeater)*Pin*h;
            Eel   = Eel   + x*pel*h;
            Ecool = Ecool + x*pco*h;
            Eins  = Eins  + Ploss*h;
            E = E + (opts.etaHeater*Pin - x*(pel + pco) - Ploss)*h;
            E = max(E, Efl);
            ts(k+1) = t0 + h;  Es(k+1) = E;  Pels(k+1) = x*pel;
        end
        out.Echarge_kWh = Ein;
        out.Eel_kWh  = Eel;   out.Ecool_kWh = Ecool;  out.Eins_kWh = Eins;
        out.Eheater_kWh = Eh; out.Eleft_kWh = E;
        out.t_h = ts;  out.E_kWh = Es;  out.Pel_kW = Pels;

    %% ------------------------------------------------------------------ full 3D
    case '3d'
        if isfield(opts, 'p') && ~isempty(opts.p), p = opts.p; else, p = thermal_block_params(); end
        if ~isempty(opts.mesh)
            p.Nx = opts.mesh(1); p.Ny = opts.mesh(2); p.Nz = opts.mesh(3);
            p.dx = p.Lx/p.Nx; p.dy = p.Ly/p.Ny; p.dz = p.Lz/p.Nz;
        end
        sc = opts.scale;                                 % 3D model is one block: scale results
        V  = p.Lx*p.Ly*p.Lz;
        Pblk_W = opts.Pcharge_kW*1e3/sc;                 % heater power per block (grid side)
        p.T0   = sur.Tmin;
        tChs = tCh*3600;  tOps = tOp*3600;  tEnds = tEnd*3600;
        p.S    = @(t) opts.etaHeater*Pblk_W/V * (t < tChs);
        p.open = @(t) double(t >= tOps);
        p.tEnd = tEnds;
        p.nOut = max(ceil(tEnd/opts.dtOut_h), 2) + 1;
        p.breaks = [tChs tOps];
        r = thermal_block_run(p);
        Eabs0 = sur.M*p.h(sur.Tmin);                     % uniform block at T_min
        out.Echarge_kWh = Echarge_kWh;                   % drawn from grid (no Tmax cap in 3D mode)
        out.Eel_kWh   = sc*r.Eel;
        out.Ecool_kWh = sc*(r.Egap - r.Eel);             % block-side heat into the TPV stacks
        out.Ecool_coolantSide_kWh = sc*r.Ecool;          % what actually reached the coolant
        out.Eins_kWh  = sc*r.EinsBlock;                  % block-side insulation loss
        out.Eins_ambientSide_kWh = sc*r.Eins;
        out.Eheater_kWh = (1-opts.etaHeater)*Echarge_kWh;
        out.Eleft_kWh = sc*(r.Emean_J(end) - Eabs0)/3.6e6;
        out.t_h = r.t/3600;  out.E_kWh = sc*(r.Emean_J - Eabs0)/3.6e6;  out.Pel_kW = sc*r.P(:,2)/1e3;
        out.res3d = r;
    otherwise
        error('ldes_cycle: unknown mode "%s"', opts.mode);
end

out.roundTrip   = out.Eel_kWh / max(out.Echarge_kWh, eps);
out.balance_kWh = out.Echarge_kWh - (out.Eel_kWh + out.Ecool_kWh + out.Eins_kWh + out.Eheater_kWh + out.Eleft_kWh);
out.mode = opts.mode;  out.tDischarge_h = tDischarge_h;  out.tStandby_h = opts.tStandby_h;
end

function s = setdef(s, f, v)
    if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
