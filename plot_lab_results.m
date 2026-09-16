function plot_lab_results(cfg,polar,cases,folder)
% Create lecture-ready PNG figures using only base MATLAB graphics.
C=palette();
f=styled_figure([100 100 1380 760]);
tl=tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
title(tl,sprintf('NACA 0012 teaching case | 2-D laminar flow | Re = %g',cfg.Re), ...
    'FontWeight','bold');
ax=nexttile(tl); hold(ax,'on');
errorbar(ax,polar.alpha_deg,polar.CL,polar.CL_std,'o-','LineWidth',1.8, ...
    'Color',C.blue,'MarkerFaceColor',C.blue,'CapSize',8);
yline(ax,0,'Color',C.gray,'LineWidth',0.8); style_axes(ax);
xlabel(ax,'Angle of attack, alpha (deg)'); ylabel(ax,'Lift coefficient, C_L');
title(ax,'Mean lift with temporal variation');
ax=nexttile(tl); hold(ax,'on');
errorbar(ax,polar.alpha_deg,polar.CD,polar.CD_std,'o-','LineWidth',1.8, ...
    'Color',C.orange,'MarkerFaceColor',C.orange,'CapSize',8);
style_axes(ax); xlabel(ax,'Angle of attack, alpha (deg)'); ylabel(ax,'Drag coefficient, C_D');
title(ax,'Mean drag with temporal variation');
ax=nexttile(tl); hold(ax,'on');
plot(ax,polar.alpha_deg,1e3*polar.tip_bending_m,'o-','LineWidth',1.8, ...
    'Color',C.teal,'MarkerFaceColor',C.teal); style_axes(ax);
xlabel(ax,'Angle of attack, alpha (deg)'); ylabel(ax,'Tip bending (mm)');
title(ax,'One-way structural response');
ax=nexttile(tl); hold(ax,'on');
plot(ax,polar.alpha_deg,polar.tip_twist_deg,'o-','LineWidth',1.8, ...
    'Color',C.purple,'MarkerFaceColor',C.purple); style_axes(ax);
xlabel(ax,'Prescribed angle of attack (deg)'); ylabel(ax,'Tip twist (deg)');
title(ax,'No structural feedback to CFD');
exportgraphics(f,fullfile(folder,'polar_and_structure.png'),'Resolution',200); close(f);

r=cases{end}; U=0.5*(r.U(:,1:end-1)+r.U(:,2:end)); V=0.5*(r.V(1:end-1,:)+r.V(2:end,:));
[X,Y]=meshgrid(r.x,r.y);
speed=hypot(U,V)/cfg.U;
Cp=r.p/(0.5*cfg.rho*cfg.U^2);
f=styled_figure([100 100 1480 700]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
title(tl,sprintf('Mean field at alpha = %g deg | averaging window: %.1f to %.1f s', ...
    r.alphaDeg,cfg.avgStart,cfg.endTime),'FontWeight','bold');
ax=nexttile(tl); contourf(ax,X,Y,speed,30,'LineColor','none'); hold(ax,'on');
q=quiver(ax,X(1:7:end,1:7:end),Y(1:7:end,1:7:end), ...
    U(1:7:end,1:7:end)/cfg.U,V(1:7:end,1:7:end),0.75,'Color',[0.12 0.12 0.12]);
q.LineWidth=0.6; fill(ax,r.xb,r.yb,C.airfoil,'EdgeColor',C.airfoil,'LineWidth',1.0);
axis(ax,'equal'); xlim(ax,[0 cfg.Lx]); ylim(ax,[0 cfg.Ly]); colormap(ax,turbo); clim(ax,[0 1.5]);
cb=colorbar(ax); cb.Label.String='Speed magnitude, |u| / U_inf';
style_axes(ax); xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); title(ax,'Velocity magnitude and direction');
ax=nexttile(tl); contourf(ax,X,Y,Cp,30,'LineColor','none'); hold(ax,'on');
fill(ax,r.xb,r.yb,C.airfoil,'EdgeColor',C.airfoil,'LineWidth',1.0);
axis(ax,'equal'); xlim(ax,[0 cfg.Lx]); ylim(ax,[0 cfg.Ly]); colormap(ax,flipud(turbo));
clim(ax,symmetric_limits(Cp)); cb=colorbar(ax); cb.Label.String='Pressure coefficient, C_p';
style_axes(ax); xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); title(ax,'Mean pressure coefficient');
exportgraphics(f,fullfile(folder,'flow_field.png'),'Resolution',200); close(f);

f=styled_figure([100 100 1380 720]);
tl=tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
title(tl,'Transient loads, numerical checks, and static structural shape','FontWeight','bold');
ax=nexttile(tl); hold(ax,'on');
plot(ax,r.history(:,1),r.history(:,2),'Color',C.blue,'LineWidth',1.4);
plot(ax,r.history(:,1),r.history(:,3),'Color',C.orange,'LineWidth',1.4);
xline(ax,cfg.avgStart,'--','Averaging starts','Color',C.gray,'LabelVerticalAlignment','bottom');
legend(ax,{'C_L','C_D'},'Location','best','Box','off'); style_axes(ax);
xlabel(ax,'Physical time (s)'); ylabel(ax,'Force coefficient'); title(ax,'Transient aerodynamic loads');
ax=nexttile(tl); hold(ax,'on');
semilogy(ax,r.history(:,1),max(r.history(:,5),1e-16),'Color',C.blue,'LineWidth',1.2);
semilogy(ax,r.history(:,1),max(r.history(:,6),1e-16),'Color',C.orange,'LineWidth',1.2);
semilogy(ax,r.history(:,1),max(r.history(:,7),1e-16),'Color',C.purple,'LineWidth',1.2);
legend(ax,{'Divergence c/U','Pressure residual','Solid speed/U'},'Location','best','Box','off'); style_axes(ax);
xlabel(ax,'Physical time (s)'); ylabel(ax,'Magnitude'); title(ax,'Numerical diagnostics');
ax=nexttile(tl); plot(ax,r.structure.z,1e3*r.structure.bending,'Color',C.teal,'LineWidth',1.8);
style_axes(ax); xlabel(ax,'Spanwise coordinate, z (m)'); ylabel(ax,'Bending (mm)'); title(ax,'Mean-load bending shape');
ax=nexttile(tl); plot(ax,r.structure.z,r.structure.twist*180/pi,'Color',C.purple,'LineWidth',1.8);
style_axes(ax); xlabel(ax,'Spanwise coordinate, z (m)'); ylabel(ax,'Twist (deg)'); title(ax,'Mean-load twist shape');
exportgraphics(f,fullfile(folder,'history_and_span.png'),'Resolution',200); close(f);
end

function f=styled_figure(position)
f=figure('Visible','off','Position',position,'Color','white');
if isprop(f,'Theme'), f.Theme='light'; end
end

function style_axes(ax)
grid(ax,'on'); box(ax,'on'); ax.FontName='Arial'; ax.FontSize=11;
ax.GridColor=[0.85 0.87 0.90]; ax.GridAlpha=0.7; ax.LineWidth=0.8;
end

function lim=symmetric_limits(values)
m=max(abs(values(~isnan(values)))); lim=[-m m];
end

function C=palette()
C.blue=[0.00 0.35 0.70]; C.orange=[0.90 0.38 0.05]; C.teal=[0.00 0.55 0.50];
C.purple=[0.47 0.27 0.65]; C.gray=[0.35 0.38 0.42]; C.airfoil=[0.08 0.12 0.18];
end
