function cfg = student_case_config(mode)
%STUDENT_CASE_CONFIG Parameters students may change for a new physical case.
%
% This is the only file beginners need to edit. Leave the solver files
% unchanged and run RUN_STUDENT_CASE afterwards.
%
% Example:
%   cfg = student_case_config('demo');
%   cfg.alphaDeg = 6;
%   cfg.U = 2;
%   cfg.nu = cfg.U*cfg.c/cfg.Re;  % Keep Re fixed after changing U.

if nargin == 0, mode = 'demo'; end
cfg = teaching_config(mode);

%% A. FLOW AND BOUNDARY CONDITIONS
cfg.U = 1;                    % Uniform horizontal inlet velocity U_inf [m/s]
cfg.Re = 100;                 % Reynolds number U_inf*c/nu [-]
cfg.nu = cfg.U*cfg.c/cfg.Re;  % Kinematic viscosity [m^2/s], derived from U, c, Re

%% B. AIRFOIL GEOMETRY AND PRESCRIBED INCIDENCE
cfg.alphaDeg = 8;             % One prescribed geometric angle of attack [deg]
cfg.thickness = 0.12;         % NACA 00xx thickness ratio: 0.12 gives NACA 0012
cfg.center = [2*cfg.c 2*cfg.c]; % Fixed-section centre [m m]

%% C. DOMAIN AND RESOLUTION
% Keep dx=dy. Therefore scale nx and ny in the same 3:2 ratio.
cfg.Lx = 6*cfg.c; cfg.Ly = 4*cfg.c;
cfg.nx = 240; cfg.ny = 160;

%% D. ONE-WAY STRUCTURAL PARAMETERS
cfg.EI = 50;                  % Bending stiffness [N m^2]
cfg.GJ = 10;                  % Torsional stiffness [N m^2]
cfg.dampingRatio = [0.03 0.04];

%% E. OPTIONAL OUTPUTS
cfg.makeGridPreview = true;   % See grid + immersed-boundary mask before marching.
cfg.makePlots = true;
cfg.saveSnapshots = false;    % Set true only when exporting a transient video.
end
