function sur = ldes_calibrate(p, opts)
% LDES_CALIBRATE  Build a fast lumped surrogate of ONE carbon block from the 3D model.
%
% The annual dispatch has 35,040 fifteen-minute steps, so it cannot call the 3D ODE
% model every step. Instead the block is reduced to one state, E = block energy in kWh
% above a uniform block at Tmin, with these tabulated flows:
%
%   Pel(E), Pcool(E)   from a 3D DISCHARGE run (block uniform at Tmax, TPV faces open):
%                      electricity out, and heat into the TPV stacks (gap flow - electric)
%                      that ends up in the coolant. Block-side flows, so the lumped
%                      balance tracks the 3D block enthalpy.
%   PlossClosed(E)     standby (shutters closed) loss through ALL six insulation stacks,
%   PlossOpen(E)       loss through the insulation-only faces while the TPV faces are open,
%                      both from the STEADY 1D solution through the same insulation stack
%                      (same k(T), thicknesses, faceScale, Rc, hOut as the 3D model).
%
% Why steady insulation? The 3D model starts the insulation on a linear profile, which is
% far from steady at 2000 K: the block then pours heat into the insulation for days
% (soak transient) and a short 3D standby run overstates self-discharge ~2-4x. A battery
% that cycles daily sits close to the steady insulation state, because the insulation
% time constant (~L^2/alpha, days) is much longer than a cycle. Set opts.standby3D = true
% to also run the 3D standby case and print it next to the steady value as a check.
%
% Lumped balance used by ldes_cycle and ldes_dispatch (x = shutter opening, 0..1):
%   dE/dt = eta_heater*P_charge - x*(Pel + Pcool) - x*PlossOpen - (1-x)*PlossClosed
%
% Usage:  sur = ldes_calibrate();                              % defaults, full 3D mesh
%         sur = ldes_calibrate([], struct('mesh',[8 8 12]))    % coarser, faster
%
% opts (defaults): Tmax 2000 K, Tmin 1000 K (the 3D model's TPV eta is zero below 1000 K),
%   tDis_h 60, dtOut_h 0.25, S_charge 1e5 W/m3 (heater rating in thermal_block_3D_full.m),
%   mesh [], insNodes [], standby3D false, tStandby_h 168, nGrid 400,
%   file 'ldes_surrogate.mat', plot true.

if nargin < 1 || isempty(p),    p = thermal_block_params(); end
if nargin < 2,                   opts = struct(); end
opts = setdef(opts, 'Tmax', 2000);          % K   max block temperature (1727 C)
opts = setdef(opts, 'Tmin', 1000);          % K   "empty": TPV output reaches zero
opts = setdef(opts, 'tDis_h', 60);
opts = setdef(opts, 'dtOut_h', 0.25);
opts = setdef(opts, 'S_charge', 1e5);       % W/m3 heater rating
opts = setdef(opts, 'mesh', []);            % [Nx Ny Nz] override (even numbers)
opts = setdef(opts, 'insNodes', []);        % insulation nodes override, e.g. [3 2]
opts = setdef(opts, 'standby3D', false);
opts = setdef(opts, 'tStandby_h', 168);
opts = setdef(opts, 'nGrid', 400);
opts = setdef(opts, 'file', 'ldes_surrogate.mat');
opts = setdef(opts, 'plot', true);

if ~isempty(opts.mesh)
    p.Nx = opts.mesh(1); p.Ny = opts.mesh(2); p.Nz = opts.mesh(3);
    p.dx = p.Lx/p.Nx; p.dy = p.Ly/p.Ny; p.dz = p.Lz/p.Nz;
end
if ~isempty(opts.insNodes), p.ins.nodes = opts.insNodes; end

V = p.Lx*p.Ly*p.Lz;   M = p.rho*V;
Eu     = @(T) M*(p.h(T) - p.h(opts.Tmin))/3.6e6;     % uniform-block energy above Tmin (kWh)
Ecap   = Eu(opts.Tmax);
Efloor = Eu(p.Tamb);                                   % negative: block cooled to ambient
Eabs0  = M*p.h(opts.Tmin);
Tg  = linspace(p.Tamb, opts.Tmax, 2000)';  EgU = Eu(Tg);
T_of_E = @(E) interp1(EgU, Tg, min(max(E, EgU(1)), EgU(end)));

fprintf('Calibrating surrogate: block %.1f m3, %.1f t, capacity %.0f kWh_th (%.0f-%.0f K)\n', ...
        V, M/1e3, Ecap, opts.Tmin, opts.Tmax);

%% 3D discharge run: faces open from a uniform block at Tmax
pA = p;  pA.T0 = opts.Tmax;  pA.S = @(t) 0;  pA.open = @(t) 1;
pA.tEnd = opts.tDis_h*3600;  pA.nOut = round(opts.tDis_h/opts.dtOut_h) + 1;
pA.breaks = [1 2 5 10 15 20 30 45 60 90]*60;          % extra saved points while the surface drops fast
fprintf('3D discharge run (%g h)...\n', opts.tDis_h);
rA = thermal_block_run(pA);
kA = rA.t > 30;                                        % t = 0 has cold TPV cells; drop it
EA    = (rA.Emean_J(kA) - Eabs0)/3.6e6;
PelA  = rA.P(kA,2)/1e3;   PgapA = rA.P(kA,5)/1e3;
[EAs, iA] = uniqueSorted(EA);   PelA = PelA(iA);   PgapA = PgapA(iA);

%% Tables on a uniform E grid
Eg = linspace(Efloor, Ecap, opts.nGrid)';
Tm = T_of_E(Eg);
% Inside the run: interpolate. Above it (first ~30 s): hold. Below it (block colder than
% the run reached): Pel ramps to 0 at E = 0, the gap heat ramps to 0 at Efloor.
Pel   = tab(Eg, EAs, PelA,  0,      0);
Pgap  = tab(Eg, EAs, PgapA, Efloor, 0);
Pcool = max(Pgap - Pel, 0);

% Steady insulation losses, block surface taken at the mean block temperature
fl  = [p.radXm p.radXp p.radYm p.radYp p.radZm p.radZp];
PlC = steadyInsLoss(p, Tm, true(1,6))/1e3;
PlO = steadyInsLoss(p, Tm, ~fl)/1e3;

%% Optional 3D standby check
chk = [];
if opts.standby3D
    pB = p;  pB.T0 = opts.Tmax;  pB.S = @(t) 0;  pB.open = @(t) 0;
    pB.tEnd = opts.tStandby_h*3600;  pB.nOut = round(opts.tStandby_h/opts.dtOut_h) + 1;
    fprintf('3D standby check (%g h)...\n', opts.tStandby_h);
    rB  = thermal_block_run(pB);
    EBe = (rB.Emean_J(end) - Eabs0)/3.6e6;
    chk = [rB.P(end,6) rB.P(end,4) steadyInsLoss(p, T_of_E(EBe), true(1,6))]/1e3;
    fprintf('  after %g h: 3D block-side %.1f kW, 3D ambient-side %.1f kW, steady model %.1f kW\n', ...
            opts.tStandby_h, chk);
    fprintf('  (3D block- and ambient-side converge to the steady value as the insulation soaks)\n');
end

%% Pack
sur.E_kWh          = Eg;
sur.Pel_kW         = Pel;
sur.Pcool_kW       = Pcool;
sur.PlossOpen_kW   = PlO;
sur.PlossClosed_kW = PlC;
sur.Tmean_K        = Tm;
sur.Ecap_kWh       = Ecap;
sur.Efloor_kWh     = Efloor;
sur.dE             = Eg(2) - Eg(1);
sur.Pcharge_kW     = opts.S_charge*V/1e3;
sur.Tmin = opts.Tmin;  sur.Tmax = opts.Tmax;  sur.V = V;  sur.M = M;
sur.selfDischarge_pct_per_day_full = 100*24*PlC(end)/Ecap;
sur.selfDischarge_pct_per_day_half = 100*24*interp1(Eg, PlC, Ecap/2)/Ecap;
sur.raw.A = struct('t_h', rA.t/3600, 'E_kWh', (rA.Emean_J - Eabs0)/3.6e6, 'P_W', rA.P);
sur.standbyCheck_kW = chk;
sur.opts = opts;  sur.mesh = [p.Nx p.Ny p.Nz];  sur.created = datestr(now);

fprintf('  Peak TPV electric power (one block): %.0f kW; at 50%% energy: %.0f kW\n', ...
        max(Pel), interp1(Eg, Pel, Ecap/2));
fprintf('  Standby self-discharge: %.2f %%/day at full (%.1f kW), %.2f %%/day at half\n', ...
        sur.selfDischarge_pct_per_day_full, PlC(end), sur.selfDischarge_pct_per_day_half);
if ~isempty(opts.file)
    save(opts.file, 'sur');
    fprintf('  Saved surrogate to %s\n', opts.file);
end

if opts.plot
    figure('Name','Surrogate calibration');
    subplot(1,2,1);
    plot(Eg, Pel, Eg, Pcool, Eg, PlO, Eg, PlC, 'LineWidth',1.5); grid on;
    legend('TPV electric','Heat to coolant','Insulation loss (open)','Insulation loss (closed)','Location','northwest');
    xlabel('Block energy above T_{min} (kWh)'); ylabel('Power (kW, one block)');
    title('Surrogate tables (one block)');
    subplot(1,2,2);
    plot(rA.t/3600, (rA.Emean_J-Eabs0)/3.6e6, 'LineWidth',1.5); grid on;
    xlabel('Time (h)'); ylabel('Block energy above T_{min} (kWh)'); title('3D discharge run');
end
end

%% ---------------- helpers ----------------
function P = steadyInsLoss(p, Tsurf, faces)
% Steady 1D heat loss (W) through the insulation stacks of the selected faces, block
% surface at Tsurf (vector). Per face: Rc, then each layer with k(T) (Kirchhoff transform),
% then the outer film hOut to Tamb. Solved for the flux q by bisection.
    area = [p.Ly*p.Lz p.Ly*p.Lz p.Lx*p.Lz p.Lx*p.Lz p.Lx*p.Ly p.Lx*p.Ly];
    Tk = linspace(p.Tamb - 50, max(Tsurf(:)) + 50, 4000)';
    nL = numel(p.ins.thick);
    K  = zeros(numel(Tk), nL);                          % K_L(T) = int k_L dT
    for L = 1:nL
        kk = p.ins.k{L}(Tk);
        K(:,L) = [0; cumsum(0.5*(kk(1:end-1) + kk(2:end)).*diff(Tk))];
    end
    P = zeros(size(Tsurf));
    for f = find(faces)
        th = p.ins.thick * p.ins.faceScale(f);
        qlo = zeros(size(Tsurf));
        qhi = max(Tsurf - p.Tamb, 0) ./ (p.ins.Rc + 1/p.hOut);    % upper bound: no insulation
        for it = 1:60
            q = 0.5*(qlo + qhi);
            T = Tsurf - q*p.ins.Rc;
            for L = 1:nL
                KL = interp1(Tk, K(:,L), T) - q*th(L);
                T  = interp1(K(:,L), Tk, max(KL, K(1,L)));
            end
            r = T - p.Tamb - q/p.hOut;              % > 0: q too small
            qlo(r > 0) = q(r > 0);   qhi(r <= 0) = q(r <= 0);
        end
        P = P + area(f)*0.5*(qlo + qhi);
    end
end

function s = setdef(s, f, v)
    if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function [x, idx] = uniqueSorted(x)
    [x, idx] = sort(x(:));
    keep = [true; diff(x) > 0];
    x = x(keep);  idx = idx(keep);
end

function y = tab(Eg, Ed, Pd, Ezero, Pzero)
    % Interpolate inside [Ed(1), Ed(end)], hold the end value above, ramp linearly
    % to Pzero at Ezero below (zero below Ezero).
    y = zeros(size(Eg));
    in = Eg >= Ed(1) & Eg <= Ed(end);
    y(in) = interp1(Ed, Pd, Eg(in));
    y(Eg > Ed(end)) = Pd(end);
    lo = Eg < Ed(1);
    if any(lo)
        w = (Eg(lo) - Ezero) / max(Ed(1) - Ezero, eps);
        y(lo) = Pzero + (Pd(1) - Pzero) * max(w, 0);
    end
end
