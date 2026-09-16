function [filename,r] = run_video_demo(alphaDeg)
% Recompute transient flow and export it as MP4; do not relabel mean flow as transient.
% Example: [filename,r] = run_video_demo(8);
if nargin==0, alphaDeg=8; end
cfg=teaching_config('demo'); cfg.makePlots=false; cfg.saveSnapshots=true;
folder=fullfile(fileparts(mfilename('fullpath')),'results','video');
if ~exist(folder,'dir'), mkdir(folder); end
cfg.gridPreviewFolder=folder;
r=solve_airfoil_ns(cfg,alphaDeg);
save(fullfile(folder,sprintf('transient_alpha_%g.mat',alphaDeg)),'r','cfg');
filename=fullfile(folder,sprintf('flow_alpha_%g.mp4',alphaDeg));
info=export_flow_video(r,cfg,filename,15);
disp(info);
end
