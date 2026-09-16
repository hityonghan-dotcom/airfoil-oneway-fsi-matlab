function T = run_sensitivity
% alpha=4 deg: grid, time-step, penalty, and domain sensitivity study.
% It does not automatically claim convergence. Hold physical eta fixed when changing grid/time step.
base=teaching_config('demo'); base.makePlots=false; base.alphaDeg=4;
labels={'baseline','finer_grid','half_dt','half_eta','longer_time','larger_domain'};
data=zeros(numel(labels),8);
folder=fullfile(fileparts(mfilename('fullpath')),'results','sensitivity');
if ~exist(folder,'dir'), mkdir(folder); end
dt0=min(0.15*(base.Lx/base.nx)/base.U,0.20*(base.Lx/base.nx)^2/base.nu);
eta=base.etaRatio*base.endTime/ceil(base.endTime/dt0);
for k=1:numel(labels)
    cfg=base;
    switch labels{k}
        case 'finer_grid', cfg.nx=360; cfg.ny=240;
        case 'half_dt', cfg.dtScale=0.5;
        case 'half_eta'
        case 'longer_time', cfg.endTime=24; cfg.avgStart=16;
        case 'larger_domain'
            cfg.Lx=9; cfg.Ly=6; cfg.nx=360; cfg.ny=240;
            cfg.center=[3 3];
    end
    h=cfg.Lx/cfg.nx;
    targetdt=cfg.dtScale*min(0.15*h/cfg.U,0.20*h^2/cfg.nu);
    actualdt=cfg.endTime/ceil(cfg.endTime/targetdt);
    cfg.etaRatio=eta/actualdt;
    if strcmp(labels{k},'half_eta'), cfg.etaRatio=cfg.etaRatio/2; end
    old=fullfile(fileparts(mfilename('fullpath')),'results','demo','case_alpha_4.mat');
    if k==1 && exist(old,'file')
        saved=load(old,'r','cfg');
        expected=teaching_config('demo');
        physical={'rho','U','c','Re','nu','Lx','Ly','nx','ny','center', ...
            'thickness','subcells','dtScale','etaRatio','endTime','avgStart', ...
            'span','EI','GJ','elasticAxis'};
        for f=physical
            assert(isequal(saved.cfg.(f{1}),expected.(f{1})), ...
                'Baseline physical config changed; rerun demo.');
        end
        r=saved.r;
    else
        r=solve_airfoil_ns(cfg,4);
    end
    data(k,:)=[cfg.nx,cfg.ny,r.dt,r.eta,r.CL,r.CD,r.windowChange,r.maxDiv];
    save(fullfile(folder,[labels{k} '.mat']),'cfg','r');
    % Save after every case so long runs remain inspectable.
    T=array2table(data(1:k,:),'VariableNames',{'nx','ny','dt_s','eta_s','CL','CD','window_change','max_div_nondim'});
    T.case_name=labels(1:k)'; T=movevars(T,'case_name','Before',1);
    T.CL_change_pct=100*(T.CL-T.CL(1))/max(abs(T.CL(1)),eps);
    T.CD_change_pct=100*(T.CD-T.CD(1))/max(abs(T.CD(1)),eps);
    writetable(T,fullfile(folder,'sensitivity.csv')); disp(T);
end
end
