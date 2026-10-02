function r = run_student_case(mode)
%RUN_STUDENT_CASE Beginner entry point: edit STUDENT_CASE_CONFIG, then run this.
%
% The workflow is intentionally short:
%   1. Read and edit STUDENT_CASE_CONFIG.
%   2. This function calls the public solver interface.
%   3. Inspect r and the automatically written grid preview.
%
% For an angle-of-attack sweep use RUN_AIRFOIL_LAB('demo') afterwards.

if nargin == 0, mode = 'demo'; end
cfg = student_case_config(mode);
root = fileparts(mfilename('fullpath'));
folder = fullfile(root,'results','student_case');
if ~exist(folder,'dir'), mkdir(folder); end
cfg.gridPreviewFolder = folder;

fprintf('\nStudent case: U = %.3g m/s, Re = %.3g, alpha = %.3g deg\n', ...
    cfg.U,cfg.Re,cfg.alphaDeg);
r = airfoil_solver(cfg);
save(fullfile(folder,sprintf('student_case_alpha_%g.mat',cfg.alphaDeg)),'r','cfg');

fprintf('Mean coefficients: CL = %.5f, CD = %.5f, CM = %.5f\n',r.CL,r.CD,r.CM);
fprintf('One-way tip response: bending = %.3f mm, twist = %.4f deg\n', ...
    1e3*r.structure.tipBending,r.structure.tipTwistDeg);
fprintf('Saved result and grid preview in %s\n',folder);
end
