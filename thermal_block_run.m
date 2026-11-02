function res = thermal_block_run(p)
% THERMAL_BLOCK_RUN  Run the 3D carbon-block + insulation + TPV model.
% This is runModel() from thermal_block_3D_full.m with three changes only:
%   1) time span comes from p.tEnd / p.nOut instead of a hard-coded 0-40 h,
%   2) the power diagnostics also return the BLOCK-SIDE flows
%        P.gap      heat leaving the block across the TPV gaps (W)
%        P.insBlock heat leaving the block into the insulation (W)
%      so the block energy balance closes exactly:
%        Joule = dE/dt + gap + insBlock,  gap = electric + heat into the TPV stacks
%   3) res.Emean_J / res.Tmean are returned for the surrogate calibration.
% The physics (rhsAll, tpvColumn, insColumn, ...) is unchanged.
%
% res.P columns: [joule  el  cool(coolant side)  ins(ambient side)  gap  insBlock]  (W)

    p.ins = buildStack(p.ins);   p.ins.dxnF = p.ins.dxn * p.ins.faceScale;
    p.tpv = buildStack(p.tpv);
    Nx = p.Nx; Ny = p.Ny; Nz = p.Nz; Ni = p.ins.N; Nt = p.tpv.N;

    fl = [p.radXm p.radXp p.radYm p.radYp p.radZm p.radZp];
    p.tpv.faces = find(fl);

    dimsFace = [p.Ly p.Lz; p.Ly p.Lz; p.Lx p.Lz; p.Lx p.Lz; p.Lx p.Ly; p.Lx p.Ly];
    p.gap.Ff = p.gap.F*ones(1,6);
    if isnan(p.gap.F)
        for f = 1:6
            p.gap.Ff(f) = viewFactorRect(dimsFace(f,1), dimsFace(f,2), p.gap.width);
        end
    end

    ncol   = [Ny*Nz Ny*Nz Nx*Nz Nx*Nz Nx*Ny Nx*Ny];
    d.ncol = ncol;
    d.nb   = Nx*Ny*Nz;
    d.nI   = Ni*ncol;
    d.off  = d.nb + [0 cumsum(d.nI(1:5))];
    endI   = d.nb + sum(d.nI);
    d.nT   = Nt*ncol(p.tpv.faces);
    d.offT = endI + [0 cumsum(d.nT(1:end-1))];
    p.dim  = d;

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
    tspan = linspace(0, p.tEnd, p.nOut);
    if ~isfield(p, 'breaks'), p.breaks = []; end
    brk   = p.breaks(p.breaks > 0 & p.breaks < p.tEnd);
    tspan = unique([tspan, brk(:)']);                    % saved points land on every switch time
    tic; [t,Y] = ode15s(@(t,Y) rhsAll(t,Y,p), tspan, Y0, opts);
    if ~p.quiet, fprintf('  3D model: %.0f h simulated in %.1f s\n', p.tEnd/3600, toc); end

    nt = numel(t);
    res.t      = t;
    res.Tfield = reshape(Y(:,1:d.nb)', Nx, Ny, Nz, nt);

    E = zeros(nt,1);  Pm = zeros(nt,6);
    for m = 1:nt
        Tm = res.Tfield(:,:,:,m);
        E(m) = sum(p.rho*p.h(Tm(:)))*p.dx*p.dy*p.dz;
        [~,P] = rhsAll(t(m), Y(m,:)', p);
        Pm(m,:) = [P.joule P.el P.cool P.ins P.gap P.insBlock];
    end
    % One-sided values at switch times (heater off, shutters open/close) so the
    % energy integrals do not smear a step change across a whole output interval.
    Pminus = Pm;  Pplus = Pm;  dlt = 1e-6;
    for b = brk(:)'
        m = find(abs(t - b) < 1e-6, 1);
        if isempty(m), continue; end
        [~,P] = rhsAll(b - dlt, Y(m,:)', p);  Pminus(m,:) = [P.joule P.el P.cool P.ins P.gap P.insBlock];
        [~,P] = rhsAll(b + dlt, Y(m,:)', p);  Pplus(m,:)  = [P.joule P.el P.cool P.ins P.gap P.insBlock];
    end
    Eint = sum(diff(t) .* (Pplus(1:end-1,:) + Pminus(2:end,:))/2, 1) / 3.6e6;   % kWh

    res.Emean_J = E;                                         % absolute block enthalpy (J)
    res.Tmean   = squeeze(mean(mean(mean(res.Tfield,1),2),3));
    res.E_kWh = (E - E(1))/3.6e6;
    res.P     = Pm;
    res.Ejoule = Eint(1);  res.Eel = Eint(2);
    res.Ecool  = Eint(3);  res.Eins = Eint(4);
    res.Egap   = Eint(5);  res.EinsBlock = Eint(6);
    res.balance_kWh = res.Ejoule - res.Egap - res.EinsBlock - res.E_kWh(end);   % ~0 if converged
    res.rt     = res.Eel / max(res.Ejoule, eps);
    res.p = p;
end

function S = buildStack(S)
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
    dxb  = [p.dx p.dx p.dy p.dy p.dz p.dz];
    area = [p.dy*p.dz p.dy*p.dz p.dx*p.dz p.dx*p.dz p.dx*p.dy p.dx*p.dy];
    fl   = [p.radXm p.radXp p.radYm p.radYp p.radZm p.radZp];

    qIn = cell(1,6); dI = cell(1,6); qOutIns = cell(1,6);
    qGap = cell(1,6); qEl = cell(1,6); dTP = cell(1,6); qBack = cell(1,6);
    for f = 1:6
        gate = open*fl(f);
        [qIn{f}, dI{f}, qOutIns{f}] = insColumn(Ts{f}, p.k(Ts{f}), dxb(f), I{f}, p.ins.dxnF(:,f), p, 1 - gate);
        if fl(f)
            [qGap{f}, qEl{f}, dTP{f}, qBack{f}] = tpvColumn(Ts{f}, TP{f}, gate, f, p);
        else
            qGap{f} = zeros(1, ncol(f));
        end
    end
    Qt = cell(1,6);
    for f = 1:6
        Qt{f} = qIn{f} + qGap{f};
    end

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

    if nargout > 1
        P.joule = p.S(t)*p.Lx*p.Ly*p.Lz;
        P.ins = 0;  P.el = 0;  P.cool = 0;  P.gap = 0;  P.insBlock = 0;
        for f = 1:6
            P.ins      = P.ins      + sum(qOutIns{f})*area(f);
            P.insBlock = P.insBlock + sum(qIn{f})*area(f);
            if fl(f)
                P.el   = P.el   + sum(qEl{f})*area(f);
                P.cool = P.cool + sum(qBack{f})*area(f);
                P.gap  = P.gap  + sum(qGap{f})*area(f);
            end
        end
    end
end

function [qGap, qEl, dT, qBack] = tpvColumn(Ts, T, gate, f, p)
    Nt = p.tpv.N;  dxn = p.tpv.dxn;
    Tc = T(1,:);
    Rrad = (1 - p.epsB)/p.epsB + 1/p.gap.Ff(f) + (1 - p.tpv.epsC)/p.tpv.epsC;
    qRad = p.sigma*(Ts.^4 - Tc.^4) / Rrad;
    qCnd = p.gap.k(0.5*(Ts + Tc)) .* (Ts - Tc) / p.gap.width;
    qEl  = gate * p.tpv.eta(Ts, Tc) .* max(qRad, 0);
    qGap = gate * (qRad + qCnd);

    [K, C] = stackProps(T, p.tpv);
    Rint  = dxn(1:end-1)./(2*K(1:end-1,:)) + dxn(2:end)./(2*K(2:end,:));
    Rout  = dxn(end)./(2*K(end,:)) + 1/p.tpv.hCool;
    qBack = (T(end,:) - p.tpv.Tcool) ./ Rout;

    Q = zeros(Nt+1, size(T,2));
    Q(1,:)    = qGap - qEl;
    Q(2:Nt,:) = (T(1:end-1,:) - T(2:end,:)) ./ Rint;
    Q(end,:)  = qBack;
    dT = (Q(1:end-1,:) - Q(2:end,:)) ./ (dxn .* C);
end

function [qIn, dI, qOut] = insColumn(Ts, kb, dxb, I, dxn, p, conn)
    Ni = numel(dxn);
    [K, C] = stackProps(I, p.ins);

    Rin  = dxb./(2*kb) + dxn(1)./(2*K(1,:)) + p.ins.Rc;
    qIn  = conn * (Ts - I(1,:)) ./ Rin;
    Rint = dxn(1:end-1)./(2*K(1:end-1,:)) + dxn(2:end)./(2*K(2:end,:));
    Rout = dxn(end)./(2*K(end,:)) + 1/p.hOut;
    qOut = (I(end,:) - p.Tamb) ./ Rout;

    Q = zeros(Ni+1, size(I,2));
    Q(1,:)    = qIn;
    Q(2:Ni,:) = (I(1:end-1,:) - I(2:end,:)) ./ Rint;
    Q(end,:)  = qOut;
    dI = (Q(1:end-1,:) - Q(2:end,:)) ./ (dxn .* C);
end

function [K, C] = stackProps(T, S)
    K = zeros(size(T));  C = zeros(size(T));
    for L = 1:numel(S.thick)
        r = S.rows{L};
        K(r,:) = S.k{L}(T(r,:));
        C(r,:) = S.rho(L) * S.cp{L}(T(r,:));
    end
end

function F = viewFactorRect(X, Y, D)
    x = X/D;  y = Y/D;  a = sqrt(1 + x^2);  b = sqrt(1 + y^2);
    F = 2/(pi*x*y) * ( log(sqrt((1 + x^2)*(1 + y^2)/(1 + x^2 + y^2))) ...
        + x*b*atan(x/b) + y*a*atan(y/a) - x*atan(x) - y*atan(y) );
end

function b = surfIdx(f, Nx, Ny, Nz)
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
    for f = 1:6
        c = (1:d.ncol(f))';
        b = [b; surfIdx(f,Nx,Ny,Nz)];                     %#ok<AGROW>
        n = [n; d.off(f) + 1 + (c-1)*Ni];                 %#ok<AGROW>
    end
    for j = 1:numel(p.tpv.faces)
        f = p.tpv.faces(j);  c = (1:d.ncol(f))';
        b = [b; surfIdx(f,Nx,Ny,Nz)];                     %#ok<AGROW>
        n = [n; d.offT(j) + 1 + (c-1)*Nt];                %#ok<AGROW>
    end
    S = S + sparse([b; n], [n; b], 1, N, N);
end

function S = sparsity3D(Nx,Ny,Nz)
    Dx = spdiags(ones(Nx,3), -1:1, Nx, Nx);
    Dy = spdiags(ones(Ny,3), -1:1, Ny, Ny);
    Dz = spdiags(ones(Nz,3), -1:1, Nz, Nz);
    S  = kron(speye(Nz), kron(speye(Ny), Dx)) ...
       + kron(speye(Nz), kron(Dy, speye(Nx))) ...
       + kron(Dz, speye(Nx*Ny));
end
