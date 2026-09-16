function [polar,cases] = run_airfoil_lab(mode)
% Run from this folder: run_airfoil_lab('smoke') or run_airfoil_lab('demo').
if nargin==0, mode='demo'; end
cfg=teaching_config(mode);
folder=fullfile(fileparts(mfilename('fullpath')),'results',mode);
if ~exist(folder,'dir'), mkdir(folder); end
cfg.gridPreviewFolder=folder;
cases=cell(numel(cfg.alphaDeg),1); data=zeros(numel(cases),13);
for k=1:numel(cases)
    r=solve_airfoil_ns(cfg,cfg.alphaDeg(k));
    r.dynamics=wing_dynamics(cfg,r); cases{k}=r;
    data(k,:)=[r.alphaDeg,r.CL,r.CD,r.CM,r.CLstd,r.CDstd, ...
        r.structure.tipBending,r.structure.tipTwistDeg,r.maxDiv,r.maxSlip, ...
        r.pressureResidual,r.windowChange,r.maxCFL];
    h=array2table(r.history,'VariableNames',{'time_s','CL','CD','CM', ...
        'div_nondim','pressure_residual','solid_speed_over_U','CFL','max_speed','lift_N_per_m'});
    writetable(h,fullfile(folder,sprintf('history_alpha_%g.csv',r.alphaDeg)));
    d=r.dynamics;
    motion=table(d.time,d.h,d.theta*180/pi,'VariableNames', ...
        {'time_s','tip_bending_m','tip_twist_deg'});
    writetable(motion,fullfile(folder,sprintf('motion_alpha_%g.csv',r.alphaDeg)));
    save(fullfile(folder,sprintf('case_alpha_%g.mat',r.alphaDeg)),'r','cfg');
end
polar=array2table(data,'VariableNames',{'alpha_deg','CL','CD','CM','CL_std','CD_std', ...
    'tip_bending_m','tip_twist_deg','max_div_nondim','solid_speed_over_U', ...
    'pressure_residual','window_change','max_CFL'});
writetable(polar,fullfile(folder,'polar.csv')); disp(polar);
if any(polar.window_change>0.02)
    warning('Some cases still drift in the averaging window. Extend endTime/avgStart before calling the polar steady.');
end
if any(abs(polar.tip_bending_m)>0.02*cfg.span | abs(polar.tip_twist_deg)>2)
    warning('Structural deformation exceeds the small-deflection teaching range. Increase EI/GJ or reduce the load.');
end
if cfg.makePlots, plot_lab_results(cfg,polar,cases,folder); end
fprintf('Results saved to %s\n',folder);
end
