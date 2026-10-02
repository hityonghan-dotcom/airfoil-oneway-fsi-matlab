function [polar,cases] = lesson_02_angle_sweep(mode)
%LESSON_02_ANGLE_SWEEP Sweep prescribed incidence after the single-case lesson.
%
% Edit cfg.alphaDeg below. The solver runs each fixed geometry separately,
% averages the loads, then sends each load history to the one-way structure.

if nargin == 0, mode = 'demo'; end
cfg = teaching_config(mode);
cfg.alphaDeg = [-4 0 4 8];   % Student experiment: choose the sweep here.
[polar,cases] = run_airfoil_lab(cfg);
end
