function cfg = teaching_config(mode)
% Teaching parameters in SI units. This is a low-Reynolds-number laminar case.
if nargin == 0, mode = 'demo'; end
cfg.mode = mode;
cfg.alphaDeg = [-4 0 4 8];
cfg.rho = 1.225; cfg.U = 1; cfg.c = 1; cfg.Re = 100;
cfg.nu = cfg.U*cfg.c/cfg.Re;
cfg.Lx = 6*cfg.c; cfg.Ly = 4*cfg.c;
cfg.nx = 240; cfg.ny = 160;
cfg.center = [2*cfg.c 2*cfg.c];
cfg.thickness = 0.12; % NACA 0012
cfg.subcells = 4; % 4-by-4 sub-sampling for each velocity-control-volume fraction
cfg.dtScale = 1; cfg.etaRatio = 0.01;
cfg.endTime = 12*cfg.c/cfg.U;
cfg.avgStart = 8*cfg.c/cfg.U;
cfg.span = 0.5; % Cantilever length; the same 2-D load is repeated along span
cfg.EI = 50; cfg.GJ = 10; % N m^2; effective teaching stiffnesses
cfg.elasticAxis = 0.25; % x/c measured from the leading edge
cfg.makePlots = true;
cfg.makeGridPreview = true; % Save a mesh/mask preview before each time march
cfg.gridPreviewResolution = 180;
cfg.saveSnapshots = false; % Video output requires transient snapshots
cfg.snapshotDt = 0.1*cfg.c/cfg.U;
cfg.massPerLength = 10; % kg/m; effective spanwise mass density
cfg.polarInertiaPerLength = 1; % kg m; polar inertia per unit span about elastic axis
cfg.dampingRatio = [0.03 0.04]; % bending and torsion modal damping ratios
cfg.motionDisplayScale = [100 20]; % Animation-only scale factors: displacement and twist
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
