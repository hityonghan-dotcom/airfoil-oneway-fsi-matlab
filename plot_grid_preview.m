function filename = plot_grid_preview(cfg,alphaDeg,xb,yb,chiU,~,filename)
% Save a Cartesian MAC-grid, geometry, and Brinkman-mask preview.
dx=cfg.Lx/cfg.nx; dy=cfg.Ly/cfg.ny;
[Xu,Yu]=meshgrid(0:dx:cfg.Lx,dy/2:dy:cfg.Ly-dy/2);
[Xv,Yv]=meshgrid(dx/2:dx:cfg.Lx-dx/2,0:dy:cfg.Ly);
f=figure('Visible','off','Color','white','Position',[100 100 1200 620]);
if isprop(f,'Theme'), f.Theme='light'; end
cleaner=onCleanup(@() close_if_open(f));
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl); hold(ax1,'on');
gridStep=max(1,round(cfg.nx/30));
for x=0:gridStep*dx:cfg.Lx, plot(ax1,[x x],[0 cfg.Ly],'Color',[0.82 0.82 0.82]); end
for y=0:gridStep*dy:cfg.Ly, plot(ax1,[0 cfg.Lx],[y y],'Color',[0.82 0.82 0.82]); end
fill(ax1,xb,yb,[0.2 0.2 0.2],'EdgeColor','k','LineWidth',1.0);
plot(ax1,cfg.center(1),cfg.center(2),'r+','MarkerSize',8,'LineWidth',1.2);
axis(ax1,'equal'); xlim(ax1,[0 cfg.Lx]); ylim(ax1,[0 cfg.Ly]); box(ax1,'on');
xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)');
title(ax1,{sprintf('NACA 0012, alpha = %g deg',alphaDeg), ...
    sprintf('Whole domain: %d x %d pressure cells, Delta = %.3g m',cfg.nx,cfg.ny,dx)});
legend(ax1,{'Mesh lines','Airfoil','Section center'},'Location','southoutside','AutoUpdate','off');
ax2=nexttile(tl); hold(ax2,'on');
imagesc(ax2,0:dx:cfg.Lx,dy/2:dy:cfg.Ly-dy/2,chiU); axis(ax2,'xy');
contour(ax2,Xu,Yu,chiU,[0.5 0.5],'w-','LineWidth',1.1);
plot(ax2,xb,yb,'k-','LineWidth',1.1);
scatter(ax2,Xu(1:2:end,1:2:end),Yu(1:2:end,1:2:end),6,[1 1 1],'.');
scatter(ax2,Xv(1:2:end,1:2:end),Yv(1:2:end,1:2:end),6,[0 0 0],'.');
axis(ax2,'equal'); xlim(ax2,cfg.center(1)+[-0.75 0.75]*cfg.c); ylim(ax2,cfg.center(2)+[-0.42 0.42]*cfg.c);
colormap(ax2,parula); clim(ax2,[0 1]); cb=colorbar(ax2); cb.Label.String='u-face solid fraction chi_u';
xlabel(ax2,'x (m)'); ylabel(ax2,'y (m)');
title(ax2,'Zoom: Brinkman mask; white dots u faces, black dots v faces');
exportgraphics(f,filename,'Resolution',cfg.gridPreviewResolution);
end

function close_if_open(f)
if isgraphics(f), close(f); end
end
