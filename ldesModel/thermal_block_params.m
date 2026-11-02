function p = thermal_block_params()
% THERMAL_BLOCK_PARAMS  Default parameters of the 3D carbon-block + TPV model.
% Same values as thermal_block_3D_full.m, packaged as a function so the
% calibration and cycle tools can call the 3D model repeatedly.
% Values are those in thermal_block_3D_full.m (that file labels its material properties as
% design placeholders; update both together).
% Face order everywhere: [x-  x+  y-  y+  z-  z+]

%% Block geometry (FULL dimensions) and mesh. Even N so the center is a cell boundary.
p.Lx = 2.44;  p.Ly = 6.1;  p.Lz = 2.6;     % m
p.Nx = 16;    p.Ny = 16;    p.Nz = 24;
p.dx = p.Lx/p.Nx;  p.dy = p.Ly/p.Ny;  p.dz = p.Lz/p.Nz;

%% Which faces carry a gap + TPV stack (1) or are insulated only (0)
p.radXm = 1;  p.radXp = 1;
p.radYm = 1;  p.radYp = 1;
p.radZm = 0;  p.radZp = 0;

%% Block constants
p.rho   = 1800;                     % kg/m3
p.sigma = 5.670e-8;                 % W/m2-K4
p.epsB  = 0.85;                     % block surface emissivity
p.Tamb  = 300;                      % K
p.T0    = 1000;                     % initial block temperature (K)

%% Block properties
p.cp = @(T) 2000 - 1900*exp(-T/600);            % J/kg-K
p.h  = @(T) 2000*T + 1.14e6*exp(-T/600);        % enthalpy (integral of cp), J/kg
p.k  = @(T) 120*(300./T).^0.8;                  % W/m-K

%% Insulation stack (layer 1 touches the block; last layer faces ambient)
p.ins.thick = [0.3 0.3];
p.ins.nodes = [6 4];
p.ins.rho   = [150 250];
p.ins.cp    = {p.cp, @(T) 1000*ones(size(T))};
p.ins.k     = {@(T) 0.08 + 4e-11*T.^3, ...
               @(T) 0.03 + 1.5e-11*T.^3};
p.ins.faceScale = [1 1 1 1 1 1];
p.ins.Rc    = 0;
p.hOut      = 10;

%% Gap between block face and TPV cell front surface
p.gap.width = 0.005;
p.gap.k     = @(T) zeros(size(T));               % vacuum; air: @(T) 0.0263*(T/300).^0.8
p.gap.F     = NaN;

%% TPV cell stack (layer 1 faces the gap; last layer touches the coolant)
p.tpv.thick = [0.0005 0.003];
p.tpv.nodes = [2 3];
p.tpv.rho   = [5300 8900];
p.tpv.cp    = {@(T) 330*ones(size(T)), @(T) 385*ones(size(T))};
p.tpv.k     = {@(T) 40*ones(size(T)),  @(T) 390*ones(size(T))};
p.tpv.epsC  = 0.3;
p.tpv.hCool = 2000;
p.tpv.Tcool = 300;
p.tpv.eta   = @(Tem,Tc) 0.40*(1 - exp(-max(Tem-1000,0)/500)) .* max(1 - 0.002*(Tc-300), 0);

%% Schedule (overwritten by the callers)
p.S    = @(t) 1.0e5 * (t < 6*3600);                    % Joule heating, W/m3
p.open = @(t) double(t >= 8*3600 && t < 30*3600);      % TPV faces exposed

%% Solver output (new: runModel used to hard-code 0-40 h, 201 points)
p.tEnd  = 40*3600;      % s
p.nOut  = 201;          % saved time points
p.quiet = false;        % true = suppress the tic/toc print
p.breaks = [];          % times (s) where p.S or p.open switch; makes the energy integrals exact
                        % e.g. for the schedule above: p.breaks = [6 8 30]*3600
end
