function cfg = teaching_config(mode)
%TEACHING_CONFIG One place for all parameters used in the classroom examples.
%
% Start here before reading the solver.  Every quantity is in SI units.
% Students can change the values marked "student experiment" without changing
% the equations or numerical implementation.  The default model is a 2-D,
% laminar, low-Reynolds-number demonstration; it is not a stall prediction.
if nargin == 0, mode = 'demo'; end

%% 1. Prescribed flow and airfoil (student experiment)
cfg.mode = mode;
cfg.alphaDeg = [-4 0 4 8];
cfg.rho = 1.225;             % Fluid density [kg/m^3]
cfg.U = 1;                   % Inlet speed, U_inf [m/s]
cfg.c = 1;                   % Airfoil chord [m]
cfg.Re = 100;                % Reynolds number, U_inf*c/nu [-]
cfg.nu = cfg.U*cfg.c/cfg.Re;

%% 2. Cartesian CFD domain and MAC grid (student experiment)
cfg.Lx = 6*cfg.c;            % Domain length [m]
cfg.Ly = 4*cfg.c;            % Domain height [m]
cfg.nx = 240;                % Number of pressure cells in x
cfg.ny = 160;                % Number of pressure cells in y
cfg.center = [2*cfg.c 2*cfg.c]; % Airfoil-section centre [m m]
cfg.thickness = 0.12;        % NACA 0012 thickness-to-chord ratio [-]
cfg.subcells = 4;            % Mask sampling: subcells-by-subcells points per face

%% 3. Time marching and Brinkman immersed boundary (numerical controls)
cfg.dtScale = 1;             % Multiplies the stable time-step estimate
cfg.etaRatio = 0.01;         % eta/dt: smaller values impose no-slip more strongly
cfg.endTime = 12*cfg.c/cfg.U;   % End of transient calculation [s]
cfg.avgStart = 8*cfg.c/cfg.U;   % Start of force averaging window [s]

%% 4. One-way structural model (student experiment)
cfg.span = 0.5;              % Cantilever span; the 2-D load is repeated along it [m]
cfg.EI = 50;                 % Effective bending stiffness [N m^2]
cfg.GJ = 10;                 % Effective torsional stiffness [N m^2]
cfg.elasticAxis = 0.25;      % Elastic-axis location from leading edge [x/c]
cfg.massPerLength = 10;      % Effective spanwise mass [kg/m]
cfg.polarInertiaPerLength = 1; % Polar inertia per span [kg m]
cfg.dampingRatio = [0.03 0.04]; % Bending and torsion damping ratios [-]

%% 5. Files and classroom visualisation
cfg.makePlots = true;
cfg.makeGridPreview = true; % Save a mesh/mask preview before each time march
cfg.gridPreviewResolution = 180;
cfg.saveSnapshots = false; % Video output requires transient snapshots
cfg.snapshotDt = 0.1*cfg.c/cfg.U;
cfg.motionDisplayScale = [100 20]; % Animation-only scale factors: displacement and twist

%% 6. Predefined runtime levels
switch lower(mode)
    case 'smoke'
        cfg.alphaDeg = 4; cfg.nx = 120; cfg.ny = 80;
        cfg.endTime = 0.1; cfg.avgStart = 0.05;
    case 'demo'
    case 'study'
        cfg.nx = 360; cfg.ny = 240;
        cfg.endTime = 24; cfg.avgStart = 16;
    otherwise
        error('Mode must be ''smoke'', ''demo'', or ''study''.');
end
end
