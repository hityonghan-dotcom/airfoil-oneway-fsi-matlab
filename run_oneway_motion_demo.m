function [filename,r] = run_oneway_motion_demo(alphaDeg)
% Reuse an existing fixed-geometry transient field where possible.
% Only the structural response is recomputed; structural parameters do not change CFD loads.
if nargin==0, alphaDeg=8; end
cfg=teaching_config('demo'); cfg.makePlots=false; cfg.saveSnapshots=true;
root=fileparts(mfilename('fullpath')); folder=fullfile(root,'results','oneway_motion');
if ~exist(folder,'dir'), mkdir(folder); end
cfg.gridPreviewFolder=folder;
old=fullfile(root,'results','video',sprintf('transient_alpha_%g.mat',alphaDeg));
if exist(old,'file')
    saved=load(old,'r','cfg');
    physical={'rho','U','c','Re','nu','Lx','Ly','nx','ny','center', ...
        'thickness','subcells','dtScale','etaRatio','endTime','avgStart','elasticAxis'};
    for f=physical
        assert(isequal(saved.cfg.(f{1}),cfg.(f{1})), ...
            'CFD config changed; rerun run_video_demo before using cached flow.');
    end
    r=saved.r;
    if cfg.makeGridPreview
        plot_grid_preview(cfg,alphaDeg,r.xb,r.yb,r.chiU,r.chiV, ...
            fullfile(folder,sprintf('grid_preview_alpha_%g.png',alphaDeg)));
    end
else
    r=solve_airfoil_ns(cfg,alphaDeg);
end
r.dynamics=wing_dynamics(cfg,r);
d=r.dynamics;
motion=table(d.time,d.h,d.theta*180/pi,d.Lprime,d.Mprime, ...
    'VariableNames',{'time_s','tip_bending_m','tip_twist_deg','lift_N_per_m','moment_N'});
writetable(motion,fullfile(folder,sprintf('motion_alpha_%g.csv',alphaDeg)));
% Save structure and configuration; large flow snapshots remain in the transient MAT file.
save(fullfile(folder,sprintf('structure_alpha_%g.mat',alphaDeg)),'d','cfg');
filename=fullfile(folder,sprintf('oneway_motion_alpha_%g.mp4',alphaDeg));
info=export_oneway_motion_video(r,cfg,filename);
disp(info);
fprintf('Natural frequencies: %.3f Hz (bending), %.3f Hz (torsion)\n',d.naturalHz);
fprintf('Peak displacement %.3f um; peak twist %.6f deg\n', ...
    max(abs(d.h))*1e6,max(abs(d.theta))*180/pi);
end
