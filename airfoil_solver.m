function r = airfoil_solver(cfg)
%AIRFOIL_SOLVER Public interface to the fixed-geometry one-way-FSI solver.
%
% This function plays the role of a C/C++ header plus public implementation:
% it defines the stable input/output contract used by student scripts.
%
% Input
%   cfg.U, cfg.Re, cfg.nu       inlet flow and viscosity
%   cfg.alphaDeg                one prescribed airfoil incidence [deg]
%   cfg.Lx, cfg.Ly, cfg.nx, cfg.ny
%                               Cartesian domain and resolution
%   cfg.EI, cfg.GJ, cfg.dampingRatio
%                               one-way structural properties
%
% Output r
%   r.history                   [t, CL, CD, CM, diagnostics ...]
%   r.CL, r.CD, r.CM            late-time mean aerodynamic coefficients
%   r.U, r.V, r.p               late-time mean CFD fields
%   r.structure                 static cantilever response to mean CFD load
%   r.dynamics                  transient structural response to CFD load history
%
% Physical contract
%   The airfoil is fixed in SOLVE_AIRFOIL_NS. Its CFD load history drives
%   WING_DYNAMICS after the CFD run. The returned deformation never feeds
%   back into the CFD calculation.

required = {'U','Re','nu','alphaDeg','Lx','Ly','nx','ny','EI','GJ'};
for k = 1:numel(required)
    assert(isfield(cfg,required{k}), 'Missing cfg.%s.', required{k});
end
assert(isscalar(cfg.alphaDeg), ...
    'AIRFOIL_SOLVER accepts one incidence. Use RUN_AIRFOIL_LAB for a sweep.');
assert(cfg.U > 0 && cfg.nu > 0 && cfg.EI > 0 && cfg.GJ > 0, ...
    'U, nu, EI, and GJ must be positive.');

r = solve_airfoil_ns(cfg,cfg.alphaDeg);
r.dynamics = wing_dynamics(cfg,r);
end
