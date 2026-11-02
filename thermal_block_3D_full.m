%% thermal_block_3D_full_tpv.m
% 3D transient model of a FULL carbon (graphite) thermal storage block with
%   (1) explicit multilayer insulation on all six faces,
%   (2) an explicit gap between the block and the TPV cells on "open" faces, and
%   (3) a layered TPV cell stack (cell + substrate + cold plate) cooled by a coolant.
%
% Heat path on an open face during the exposure window (p.open):
%   block surface cell --[radiation across gap + gas conduction]--> TPV cell front node
%   --[1D conduction through the TPV layers]--> coolant (film coefficient p.tpv.hCool)
%   A fraction eta(T_emitter, T_cell) of the absorbed radiation leaves as ELECTRICITY;
%   the rest, plus gas conduction, heats the cell stack.
% When the window is closed (p.open = 0) the face is insulated (shutter closed): the
%   insulation stack is connected and the TPV stack just relaxes toward the coolant.
% Set p.open = @(t) 1 to model permanent exposure (block loses heat to the cells always).
%
% Radiation uses the surface-resistance network for two facing gray surfaces:
%   q = sigma*(Tb^4 - Tc^4) / ((1-eb)/eb + 1/F + (1-ec)/ec)
% with F from gap width and face size (set p.gap.F to a number to override).
%
% Face order everywhere: [x-  x+  y-  y+  z-  z+]
% All properties and dimensions are PLACEHOLDERS: replace with your materials/design.

clear; clc;

%% Block geometry (FULL dimensions) and mesh. Even N so the center is a cell boundary.
p.Lx = 2.44;  p.Ly = 6.1;  p.Lz = 2.6;     % m
p.Nx = 16;    p.Ny = 16;    p.Nz = 24;
p.dx = p.Lx/p.Nx;  p.dy = p.Ly/p.Ny;  p.dz = p.Lz/p.Nz;
xc = ((1:p.Nx)-0.5)*p.dx - p.Lx/2;
yc = ((1:p.Ny)-0.5)*p.dy - p.Ly/2;

%% Which faces carry a gap + TPV stack (1) or are insulated only (0)
p.radXm = 1;  p.radXp = 1;
p.radYm = 1;  p.radYp = 1;
p.radZm = 0;  p.radZp = 0;

%% Block constants
p.rho   = 1800;                     % kg/m3
p.sigma = 5.670e-8;                 % W/m2-K4
p.epsB  = 0.85;                     % block surface emissivity (placeholder)
p.Tamb  = 300;                      % K
p.T0    = 1000;                     % initial block temperature (K)

%% Block properties
p.cp = @(T) 2000 - 1900*exp(-T/600);            % J/kg-K
p.h  = @(T) 2000*T + 1.14e6*exp(-T/600);        % enthalpy (integral of cp), J/kg
p.k  = @(T) 120*(300./T).^0.8;                  % W/m-K

%% Insulation stack (layer 1 touches the block; last layer faces ambient)
p.ins.thick = [0.3 0.3];                      % layer thicknesses (m)
p.ins.nodes = [6 4];                            % finite-volume nodes per layer
p.ins.rho   = [150 250];                        % kg/m3
p.ins.cp    = {p.cp, @(T) 1000*ones(size(T))};  % J/kg-K per layer
p.ins.k     = {@(T) 0.08 + 4e-11*T.^3, ...      % hot-side layer (carbon-fiber-like), W/m-K
               @(T) 0.03 + 1.5e-11*T.^3};       % backup layer (microporous-like), W/m-K
p.ins.faceScale = [1 1 1 1 1 1];                % per-face multiplier on every layer's thickness
p.ins.Rc    = 0;                                % block/insulation contact resistance (m2-K/W)
p.hOut      = 10;                               % outer surface film coefficient (W/m2-K)

%% Gap between block face and TPV cell front surface
p.gap.width = 0.005;                            % m
p.gap.k     = @(T) zeros(size(T));         % W/m-K, air (rough); use zeros(size(T)) for vacuum @(T) 0.0263*(T/300).^0.8
p.gap.F     = NaN;                              % view factor: NaN = compute from gap and face size

%% TPV cell stack (layer 1 faces the gap; last layer touches the coolant)
p.tpv.thick = [0.0005 0.003];                   % m: [cell + substrate, copper cold plate]
p.tpv.nodes = [2 3];
p.tpv.rho   = [5300 8900];                      % kg/m3
p.tpv.cp    = {@(T) 330*ones(size(T)), @(T) 385*ones(size(T))};   % J/kg-K
p.tpv.k     = {@(T) 40*ones(size(T)),  @(T) 390*ones(size(T))};   % W/m-K
p.tpv.epsC  = 0.3;                              % effective absorptivity of the cell front surface
p.tpv.hCool = 2000;                             % coolant film coefficient (W/m2-K)
p.tpv.Tcool = 300;                              % coolant temperature (K)
% TPV electrical efficiency of absorbed radiation (PLACEHOLDER curve): rises with emitter
% temperature, derates with cell temperature.
p.tpv.eta   = @(Tem,Tc) 0.40*(1 - exp(-max(Tem-1000,0)/500)) .* max(1 - 0.002*(Tc-300), 0);

%% Schedule
p.S    = @(t) 1.0e5 * (t < 6*3600);                    % Joule heating, W/m3 (charge 0-6 h)
p.open = @(t) double(t >= 8*3600 && t < 30*3600);      % TPV faces exposed 8-30 h

%% Run
res = runModel(p);
t = res.t;  Tfield = res.Tfield;

%% Block temperature histories
ic = p.Nx/2 + [0 1];  jc = p.Ny/2 + [0 1];  kc = p.Nz/2 + [0 1];   % the 8 cells around the center
Tcenter = squeeze(mean(mean(mean(Tfield(ic,jc,kc,:),1),2),3));
TxpFace = squeeze(mean(mean(Tfield(end,jc,kc,:),2),3));
TxmFace = squeeze(mean(mean(Tfield(1,jc,kc,:),2),3));
Tcorner = squeeze(Tfield(end,end,end,:));
figure('Name','Block temperatures');
plot(t/3600, Tcenter, t/3600, TxpFace, t/3600, TxmFace, t/3600, Tcorner, 'LineWidth',1.5);
legend('Center','x+ face center','x- face center','x+/y+/z+ corner','Location','best');
xlabel('Time (h)'); ylabel('T (K)'); grid on; title('Block temperatures');

%% TPV stack temperatures at the center of the first open face
if ~isempty(res.tpv)
    figure('Name','TPV stack');
    plot(t/3600, res.tpv(:,[1 end]), 'LineWidth',1.5);
    legend('Cell front node','Cold-plate back node','Location','best');
    xlabel('Time (h)'); ylabel('T (K)'); grid on;
    title(sprintf('TPV stack, face %d center', res.tpvFace));
end

%% Power flows and energy summary
figure('Name','Power flows');
plot(t/3600, res.P/1e3, 'LineWidth',1.5);
legend('Joule heating in','Electric out (TPV)','Heat to TPV coolant','Heat lost through insulation','Location','best');
xlabel('Time (h)'); ylabel('Power (kW, full block)'); grid on;
fprintf('\nJoule energy in:        %.2f kWh\n', res.Ejoule);
fprintf('TPV electric out:       %.2f kWh\n', res.Eel);
fprintf('Heat to TPV coolant:    %.2f kWh\n', res.Ecool);
fprintf('Heat lost via insul.:   %.2f kWh\n', res.Eins);
fprintf('Round trip (electric out / Joule in, within the simulated window): %.3f\n', res.rt);
fprintf('Block energy left above initial at end: %.2f kWh\n\n', res.E_kWh(end));

%% Block stored energy
figure('Name','Block stored energy');
plot(t/3600, res.E_kWh, 'LineWidth',1.5); grid on;
xlabel('Time (h)'); ylabel('Block energy above initial (kWh)');

%% Temperatures inside the y+ face insulation (column near the face center)
Ni = res.p.ins.N;  nodeIdx = [1, res.p.ins.nodes(1), Ni];
figure('Name','Insulation temperatures');
plot(t/3600, res.insYp(:,nodeIdx), 'LineWidth',1.5);
legend('Insulation hot face','Layer 1/2 interface','Outer surface','Location','best');
xlabel('Time (h)'); ylabel('T (K)'); grid on; title('y+ face insulation');

%% Mid-plane (z = 0) snapshots
snapTimes = [6 18 30]*3600;
figure('Name','Mid-plane snapshots');
clim = [min(Tfield(:)) max(Tfield(:))];
for s = 1:3
    idx  = find(t >= snapTimes(s), 1);
    Fmid = 0.5*(Tfield(:,:,kc(1),idx) + Tfield(:,:,kc(2),idx));
    subplot(1,3,s);
    imagesc(xc, yc, Fmid', clim); axis image; set(gca,'YDir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('z = 0 plane, t = %.0f h', t(idx)/3600));
end

%% Optional sweeps
doSweepGap = false;                              % gap width
if doSweepGap
    gaps = [0.002 0.005 0.010 0.020 0.050];
    rt = zeros(size(gaps));  Eel = rt;
    for s = 1:numel(gaps)
        ps = p;  ps.gap.width = gaps(s);
        r  = runModel(ps);  rt(s) = r.rt;  Eel(s) = r.Eel;
        fprintf('gap %5.1f mm: electric out %.2f kWh, round trip %.3f\n', gaps(s)*1e3, Eel(s), rt(s));
    end
end
doSweepTPV = false;                              % TPV cell + substrate thickness (layer 1)
if doSweepTPV
    thk = [0.0002 0.0005 0.001 0.002];
    for s = 1:numel(thk)
        ps = p;  ps.tpv.thick(1) = thk(s);
        r  = runModel(ps);
        fprintf('TPV layer 1 %.2f mm: electric out %.2f kWh, round trip %.3f\n', thk(s)*1e3, r.Eel, r.rt);
    end
end

%% ---------------- local functions ----------------
function res = runModel(p)
    p.ins = buildStack(p.ins);   p.ins.dxnF = p.ins.dxn * p.ins.faceScale;
    p.tpv = buildStack(p.tpv);
    Nx = p.Nx; Ny = p.Ny; Nz = p.Nz; Ni = p.ins.N; Nt = p.tpv.N;

    fl = [p.radXm p.radXp p.radYm p.radYp p.radZm p.radZp];
    p.tpv.faces = find(fl);                              % faces that carry a TPV stack

    % View factor from gap width and face size (equal, directly opposed rectangles)
    dimsFace = [p.Ly p.Lz; p.Ly p.Lz; p.Lx p.Lz; p.Lx p.Lz; p.Lx p.Ly; p.Lx p.Ly];
    p.gap.Ff = p.gap.F*ones(1,6);
    if isnan(p.gap.F)
        for f = 1:6
            p.gap.Ff(f) = viewFactorRect(dimsFace(f,1), dimsFace(f,2), p.gap.width);
        end
    end

    % State layout: [block; 6 insulation stacks (x- x+ y- y+ z- z+); TPV stacks on open faces]
    ncol   = [Ny*Nz Ny*Nz Nx*Nz Nx*Nz Nx*Ny Nx*Ny];      % stack columns per face
    d.ncol = ncol;
    d.nb   = Nx*Ny*Nz;
    d.nI   = Ni*ncol;
    d.off  = d.nb + [0 cumsum(d.nI(1:5))];
    endI   = d.nb + sum(d.nI);
    d.nT   = Nt*ncol(p.tpv.faces);
    d.offT = endI + [0 cumsum(d.nT(1:end-1))];
    p.dim  = d;

    % Initial condition
    Y0 = zeros(endI + sum(d.nT), 1);
    Y0(1:d.nb) = p.T0;
    for f = 1:6
        dxf   = p.ins.dxnF(:,f);
        xn    = cumsum(dxf) - dxf/2;
        Iprof = p.T0 + (p.Tamb - p.T0)*(xn/sum(dxf));
        Y0(d.off(f) + (1:d.nI(f))) = repmat(Iprof, ncol(f), 1);
    end
    for j = 1:numel(p.tpv.faces)
        Y0(d.offT(j) + (1:d.nT(j))) = p.tpv.Tcool;
    end

    opts  = odeset('JPattern', sparsityAll(p), 'RelTol',1e-4, 'AbsTol',1e-3);
    tspan = linspace(0, 40*3600, 201);
    tic; [t,Y] = ode15s(@(t,Y) rhsAll(t,Y,p), tspan, Y0, opts); toc

    nt = numel(t);
    res.t      = t;
    res.Tfield = reshape(Y(:,1:d.nb)', Nx, Ny, Nz, nt);
    c4 = centerCol(4, Nx, Ny, Nz);
    res.insYp  = Y(:, d.off(4) + (c4-1)*Ni + (1:Ni));    % y+ insulation near the face center
    res.tpv = [];  res.tpvFace = 0;
    if ~isempty(p.tpv.faces)
        f1 = p.tpv.faces(1);  c1 = centerCol(f1, Nx, Ny, Nz);
        res.tpv     = Y(:, d.offT(1) + (c1-1)*Nt + (1:Nt));
        res.tpvFace = f1;
    end

    % Block stored energy and power flows (re-evaluate fluxes at the saved times)
    E = zeros(nt,1);  Pm = zeros(nt,4);
    for m = 1:nt
        Tm = res.Tfield(:,:,:,m);
        E(m) = sum(p.rho*p.h(Tm(:)))*p.dx*p.dy*p.dz;
        [~,P] = rhsAll(t(m), Y(m,:)', p);
        Pm(m,:) = [P.joule P.el P.cool P.ins];           % W, full block
    end
    res.E_kWh = (E - E(1))/3.6e6;
    res.P     = Pm;
    kWh = @(x) trapz(t, x)/3.6e6;
    res.Ejoule = kWh(Pm(:,1));  res.Eel = kWh(Pm(:,2));
    res.Ecool  = kWh(Pm(:,3));  res.Eins = kWh(Pm(:,4));
    res.rt     = res.Eel / res.Ejoule;
    res.p = p;
end

function S = buildStack(S)
    % Per-node thickness and layer membership from a layer table
    nL = numel(S.thick);
    dxn = [];  rows = cell(nL,1);  idx = 0;
    for L = 1:nL
        nn = S.nodes(L);
        dxn = [dxn; (S.thick(L)/nn)*ones(nn,1)]; %#ok<AGROW>
        rows{L} = idx + (1:nn);
        idx = idx + nn;
    end
    S.dxn = dxn;  S.rows = rows;  S.N = idx;
end

function [dY, P] = rhsAll(t,Y,p)
    Nx = p.Nx; Ny = p.Ny; Nz = p.Nz; Ni = p.ins.N; Nt = p.tpv.N; d = p.dim;
    ncol = d.ncol;
    Tb = reshape(Y(1:d.nb), Nx, Ny, Nz);
    I  = cell(1,6);
    for f = 1:6
        I{f} = reshape(Y(d.off(f) + (1:d.nI(f))), Ni, ncol(f));
    end
    TP = cell(1,6);
    for j = 1:numel(p.tpv.faces)
        f = p.tpv.faces(j);
        TP{f} = reshape(Y(d.offT(j) + (1:d.nT(j))), Nt, ncol(f));
    end

    open = p.open(t);
    kk   = p.k(Tb);
    Ts   = {reshape(Tb(1,:,:),1,[]),   reshape(Tb(end,:,:),1,[]), ...
            reshape(Tb(:,1,:),1,[]),   reshape(Tb(:,end,:),1,[]), ...
            reshape(Tb(:,:,1),1,[]),   reshape(Tb(:,:,end),1,[])};
    dxb  = [p.dx p.dx p.dy p.dy p.dz p.dz];             % block cell size normal to each face
    area = [p.dy*p.dz p.dy*p.dz p.dx*p.dz p.dx*p.dz p.dx*p.dy p.dx*p.dy];   % surface-cell face area
    fl   = [p.radXm p.radXp p.radYm p.radYp p.radZm p.radZp];

    qIn = cell(1,6); dI = cell(1,6); qOutIns = cell(1,6);
    qGap = cell(1,6); qEl = cell(1,6); dTP = cell(1,6); qBack = cell(1,6);
    for f = 1:6
        gate = open*fl(f);                               % 1 while this face is exposed to its TPV cells
        [qIn{f}, dI{f}, qOutIns{f}] = insColumn(Ts{f}, p.k(Ts{f}), dxb(f), I{f}, p.ins.dxnF(:,f), p, 1 - gate);
        if fl(f)
            [qGap{f}, qEl{f}, dTP{f}, qBack{f}] = tpvColumn(Ts{f}, TP{f}, gate, f, p);
        else
            qGap{f} = zeros(1, ncol(f));
        end
    end
    Qt = cell(1,6);                                      % total heat leaving the block per face (W/m2)
    for f = 1:6
        Qt{f} = qIn{f} + qGap{f};
    end

    % Block conduction (flux positive in +direction; at low faces the outflow is minus the flux)
    qx = zeros(Nx+1, Ny, Nz);
    qx(2:Nx,:,:) = -0.5*(kk(1:end-1,:,:) + kk(2:end,:,:)) .* (Tb(2:end,:,:) - Tb(1:end-1,:,:)) / p.dx;
    qx(1,:,:)    = -reshape(Qt{1}, 1, Ny, Nz);
    qx(end,:,:)  =  reshape(Qt{2}, 1, Ny, Nz);

    qy = zeros(Nx, Ny+1, Nz);
    qy(:,2:Ny,:) = -0.5*(kk(:,1:end-1,:) + kk(:,2:end,:)) .* (Tb(:,2:end,:) - Tb(:,1:end-1,:)) / p.dy;
    qy(:,1,:)    = -reshape(Qt{3}, Nx, 1, Nz);
    qy(:,end,:)  =  reshape(Qt{4}, Nx, 1, Nz);

    qz = zeros(Nx, Ny, Nz+1);
    qz(:,:,2:Nz) = -0.5*(kk(:,:,1:end-1) + kk(:,:,2:end)) .* (Tb(:,:,2:end) - Tb(:,:,1:end-1)) / p.dz;
    qz(:,:,1)    = -reshape(Qt{5}, Nx, Ny);
    qz(:,:,end)  =  reshape(Qt{6}, Nx, Ny);

    divq = (qx(2:end,:,:) - qx(1:end-1,:,:))/p.dx ...
         + (qy(:,2:end,:) - qy(:,1:end-1,:))/p.dy ...
         + (qz(:,:,2:end) - qz(:,:,1:end-1))/p.dz;
    dTb = (-divq + p.S(t)) ./ (p.rho*p.cp(Tb));

    dY = zeros(numel(Y), 1);
    dY(1:d.nb) = dTb(:);
    for f = 1:6
        dY(d.off(f) + (1:d.nI(f))) = dI{f}(:);
    end
    for j = 1:numel(p.tpv.faces)
        f = p.tpv.faces(j);
        dY(d.offT(j) + (1:d.nT(j))) = dTP{f}(:);
    end

    if nargout > 1                                       % diagnostics: total powers (W, full block)
        P.joule = p.S(t)*p.Lx*p.Ly*p.Lz;
        P.ins = 0;  P.el = 0;  P.cool = 0;
        for f = 1:6
            P.ins = P.ins + sum(qOutIns{f})*area(f);
            if fl(f)
                P.el   = P.el   + sum(qEl{f})*area(f);
                P.cool = P.cool + sum(qBack{f})*area(f);
            end
        end
    end
end

function [qGap, qEl, dT, qBack] = tpvColumn(Ts, T, gate, f, p)
    % T is Nt x ncols: node 1 faces the gap, node Nt touches the coolant.
    % Ts: block surface-cell temperature (1 x ncols). gate: 1 = exposed, 0 = shuttered.
    Nt = p.tpv.N;  dxn = p.tpv.dxn;
    Tc = T(1,:);
    Rrad = (1 - p.epsB)/p.epsB + 1/p.gap.Ff(f) + (1 - p.tpv.epsC)/p.tpv.epsC;
    qRad = p.sigma*(Ts.^4 - Tc.^4) / Rrad;                       % radiation across the gap
    qCnd = p.gap.k(0.5*(Ts + Tc)) .* (Ts - Tc) / p.gap.width;    % gas conduction across the gap
    qEl  = gate * p.tpv.eta(Ts, Tc) .* max(qRad, 0);             % electricity leaving the cell
    qGap = gate * (qRad + qCnd);                                 % block -> cell (0 when shuttered)

    [K, C] = stackProps(T, p.tpv);
    Rint  = dxn(1:end-1)./(2*K(1:end-1,:)) + dxn(2:end)./(2*K(2:end,:));
    Rout  = dxn(end)./(2*K(end,:)) + 1/p.tpv.hCool;
    qBack = (T(end,:) - p.tpv.Tcool) ./ Rout;                    % heat to coolant

    Q = zeros(Nt+1, size(T,2));
    Q(1,:)    = qGap - qEl;                                      % heat deposited in the cell
    Q(2:Nt,:) = (T(1:end-1,:) - T(2:end,:)) ./ Rint;
    Q(end,:)  = qBack;
    dT = (Q(1:end-1,:) - Q(2:end,:)) ./ (dxn .* C);
end

function [qIn, dI, qOut] = insColumn(Ts, kb, dxb, I, dxn, p, conn)
    % I is Ni x ncols: node 1 touches the block, node Ni faces ambient.
    % conn: 1 = stack connected to the block, 0 = adiabatic inner side.
    Ni = numel(dxn);
    [K, C] = stackProps(I, p.ins);

    Rin  = dxb./(2*kb) + dxn(1)./(2*K(1,:)) + p.ins.Rc;          % block cell center -> node 1
    qIn  = conn * (Ts - I(1,:)) ./ Rin;                          % W/m2 into the stack
    Rint = dxn(1:end-1)./(2*K(1:end-1,:)) + dxn(2:end)./(2*K(2:end,:));
    Rout = dxn(end)./(2*K(end,:)) + 1/p.hOut;                    % last node -> ambient
    qOut = (I(end,:) - p.Tamb) ./ Rout;

    Q = zeros(Ni+1, size(I,2));
    Q(1,:)    = qIn;
    Q(2:Ni,:) = (I(1:end-1,:) - I(2:end,:)) ./ Rint;
    Q(end,:)  = qOut;
    dI = (Q(1:end-1,:) - Q(2:end,:)) ./ (dxn .* C);              % C = rho*cp (J/m3-K)
end

function [K, C] = stackProps(T, S)
    % Layer-wise conductivity and volumetric heat capacity at each node
    K = zeros(size(T));  C = zeros(size(T));
    for L = 1:numel(S.thick)
        r = S.rows{L};
        K(r,:) = S.k{L}(T(r,:));
        C(r,:) = S.rho(L) * S.cp{L}(T(r,:));
    end
end

function F = viewFactorRect(X, Y, D)
    % View factor between two identical, directly opposed X-by-Y rectangles separated by D
    x = X/D;  y = Y/D;  a = sqrt(1 + x^2);  b = sqrt(1 + y^2);
    F = 2/(pi*x*y) * ( log(sqrt((1 + x^2)*(1 + y^2)/(1 + x^2 + y^2))) ...
        + x*b*atan(x/b) + y*a*atan(y/a) - x*atan(x) - y*atan(y) );
end

function c = centerCol(f, Nx, Ny, Nz)
    % Stack column index of the surface cell nearest the center of face f
    switch f
        case {1,2},  c = Ny/2 + (Nz/2 - 1)*Ny;
        case {3,4},  c = Nx/2 + (Nz/2 - 1)*Nx;
        otherwise,   c = Nx/2 + (Ny/2 - 1)*Nx;
    end
end

function b = surfIdx(f, Nx, Ny, Nz)
    % Linear indices of the block surface cells on face f, in stack-column order
    switch f
        case {1,2}
            [J,K] = ndgrid(1:Ny, 1:Nz);  ix = 1 + (f-1)*(Nx-1);
            b = ix + (J(:)-1)*Nx + (K(:)-1)*Nx*Ny;
        case {3,4}
            [I,K] = ndgrid(1:Nx, 1:Nz);  jy = 1 + (f-3)*(Ny-1);
            b = I(:) + (jy-1)*Nx + (K(:)-1)*Nx*Ny;
        otherwise
            [I,J] = ndgrid(1:Nx, 1:Ny);  kz = 1 + (f-5)*(Nz-1);
            b = I(:) + (J(:)-1)*Nx + (kz-1)*Nx*Ny;
    end
end

function S = sparsityAll(p)
    % Jacobian sparsity: block neighbors, all stack chains, and block <-> node 1 couplings
    Nx = p.Nx; Ny = p.Ny; Nz = p.Nz; Ni = p.ins.N; Nt = p.tpv.N; d = p.dim;
    Dn = spdiags(ones(Ni,3), -1:1, Ni, Ni);
    Dt = spdiags(ones(Nt,3), -1:1, Nt, Nt);
    blocks = {sparsity3D(Nx,Ny,Nz)};
    for f = 1:6
        blocks{end+1} = kron(speye(d.ncol(f)), Dn); %#ok<AGROW>
    end
    for j = 1:numel(p.tpv.faces)
        blocks{end+1} = kron(speye(d.ncol(p.tpv.faces(j))), Dt); %#ok<AGROW>
    end
    S = blkdiag(blocks{:});
    N = size(S,1);

    b = [];  n = [];
    for f = 1:6                                           % block <-> insulation node 1
        c = (1:d.ncol(f))';
        b = [b; surfIdx(f,Nx,Ny,Nz)];                     %#ok<AGROW>
        n = [n; d.off(f) + 1 + (c-1)*Ni];                 %#ok<AGROW>
    end
    for j = 1:numel(p.tpv.faces)                          % block <-> TPV node 1
        f = p.tpv.faces(j);  c = (1:d.ncol(f))';
        b = [b; surfIdx(f,Nx,Ny,Nz)];                     %#ok<AGROW>
        n = [n; d.offT(j) + 1 + (c-1)*Nt];                %#ok<AGROW>
    end
    S = S + sparse([b; n], [n; b], 1, N, N);
end

function S = sparsity3D(Nx,Ny,Nz)
    % Jacobian sparsity for the block: each cell couples to its 6 neighbors
    Dx = spdiags(ones(Nx,3), -1:1, Nx, Nx);
    Dy = spdiags(ones(Ny,3), -1:1, Ny, Ny);
    Dz = spdiags(ones(Nz,3), -1:1, Nz, Nz);
    S  = kron(speye(Nz), kron(speye(Ny), Dx)) ...
       + kron(speye(Nz), kron(Dy, speye(Nx))) ...
       + kron(Dz, speye(Nx*Ny));
end