function info = export_flow_video(r,cfg,filename,fps)
% Export a 16:9 teaching video from actual transient flow snapshots.
% The rendering uses base MATLAB VideoWriter and never substitutes mean flow for transient data.
if nargin<4, fps=15; end
assert(isfield(r,'snapshots') && numel(r.snapshots.time)>1, ...
    'Transient snapshots are required. Re-run with cfg.saveSnapshots = true.');
[folder,~,ext]=fileparts(filename);
if ~isempty(folder) && ~exist(folder,'dir'), mkdir(folder); end
if strcmpi(ext,'.mp4')
    writer=VideoWriter(filename,'MPEG-4');
elseif strcmpi(ext,'.avi')
    writer=VideoWriter(filename,'Motion JPEG AVI');
else
    error('Choose an .mp4 or .avi file extension.');
end
writer.FrameRate=fps; writer.Quality=95; open(writer);
f=figure('Visible','off','Color','white','Position',[50 50 1920 1080]);
if isprop(f,'Theme'), f.Theme='light'; end
cleaner=onCleanup(@() cleanup_video(writer,f));
C=palette(); [X,Y]=meshgrid(r.x,r.y);
tl=tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');
title(tl,'NACA 0012: transient incompressible flow on a Cartesian immersed-boundary grid', ...
    'FontWeight','bold');
axFlow=nexttile(tl,[2 2]);
field=imagesc(axFlow,r.x,r.y,zeros(size(X))); axis(axFlow,'xy'); hold(axFlow,'on');
% Opaque patch covers values inside the virtual solid without a transparent NaN halo.
fill(axFlow,r.xb,r.yb,C.airfoil,'EdgeColor',C.airfoil,'LineWidth',1.1);
arrows=quiver(axFlow,X(1:8:end,1:8:end),Y(1:8:end,1:8:end), ...
    zeros(size(X(1:8:end,1:8:end))),zeros(size(Y(1:8:end,1:8:end))),0.7, ...
    'Color',[0.08 0.08 0.08],'LineWidth',0.55);
axis(axFlow,'equal'); xlim(axFlow,[0 cfg.Lx]); ylim(axFlow,[0 cfg.Ly]);
colormap(axFlow,turbo); clim(axFlow,[0 1.5]); cb=colorbar(axFlow); cb.Label.String='Speed magnitude, |u| / U_inf';
style_axes(axFlow); xlabel(axFlow,'x (m)'); ylabel(axFlow,'y (m)');
axLoad=nexttile(tl); hold(axLoad,'on');
hL=plot(axLoad,NaN,NaN,'Color',C.blue,'LineWidth',1.6);
hD=plot(axLoad,NaN,NaN,'Color',C.orange,'LineWidth',1.6);
marker=plot(axLoad,NaN,NaN,'o','MarkerSize',5,'MarkerFaceColor',C.blue,'MarkerEdgeColor','w');
markerD=plot(axLoad,NaN,NaN,'o','MarkerSize',5,'MarkerFaceColor',C.orange,'MarkerEdgeColor','w');
legend(axLoad,[hL hD],{'C_L','C_D'},'Location','northeast','Box','off','AutoUpdate','off');
stable=r.history(:,1)>=min(0.1*cfg.c/cfg.U,cfg.endTime/2); bounds=r.history(stable,2:3);
pad=max(0.05,0.12*(max(bounds(:))-min(bounds(:)))); ylim(axLoad,[min(0,min(bounds(:))-pad) max(bounds(:))+pad]);
xlim(axLoad,[0 cfg.endTime]); style_axes(axLoad); xlabel(axLoad,'Physical time (s)'); ylabel(axLoad,'Force coefficient');
title(axLoad,'Instantaneous aerodynamic coefficients');
axInfo=nexttile(tl); axis(axInfo,'off');
infoText=text(axInfo,0.04,0.80,'','Units','normalized','FontName','Arial','FontSize',13, ...
    'VerticalAlignment','top','Interpreter','none');
temp=[tempname '.png']; tempCleaner=onCleanup(@() delete_if_exists(temp));
for k=1:numel(r.snapshots.time)
    speed=hypot(double(r.snapshots.U(:,:,k)),double(r.snapshots.V(:,:,k)))/cfg.U;
    field.CData=speed;
    arrows.UData=double(r.snapshots.U(1:8:end,1:8:end,k))/cfg.U;
    arrows.VData=double(r.snapshots.V(1:8:end,1:8:end,k))/cfg.U;
    t=r.snapshots.time(k); title(axFlow,sprintf('Re = %g | alpha = %g deg | flow snapshot at t = %.3f s',cfg.Re,r.alphaDeg,t));
    subset=r.history(:,1)<=t+1e-12;
    hL.XData=r.history(subset,1); hL.YData=r.history(subset,2);
    hD.XData=r.history(subset,1); hD.YData=r.history(subset,3);
    idx=find(subset,1,'last'); marker.XData=r.history(idx,1); marker.YData=r.history(idx,2);
    markerD.XData=r.history(idx,1); markerD.YData=r.history(idx,3);
    infoText.String=sprintf(['TIME\n%.3f s\n\nCOEFFICIENTS\nC_L = %+0.4f\nC_D = %0.4f\n\n' ...
        'GRID\n%d x %d pressure cells\nDelta = %.4f m'], ...
        t,r.history(idx,2),r.history(idx,3),cfg.nx,cfg.ny,cfg.Lx/cfg.nx);
    print(f,temp,'-dpng','-r110'); writeVideo(writer,imread(temp));
end
close(writer); close(f);
info=struct('filename',filename,'frames',numel(r.snapshots.time),'fps',fps, ...
    'playbackSeconds',numel(r.snapshots.time)/fps,'physicalTimeEnd',r.snapshots.time(end));
end

function style_axes(ax)
grid(ax,'on'); box(ax,'on'); ax.FontName='Arial'; ax.FontSize=11; ax.LineWidth=0.8;
ax.GridColor=[0.85 0.87 0.90]; ax.GridAlpha=0.7;
end

function C=palette()
C.blue=[0.00 0.35 0.70]; C.orange=[0.90 0.38 0.05];
C.gray=[0.30 0.34 0.39]; C.airfoil=[0.08 0.12 0.18];
end

function cleanup_video(writer,f)
try
    close(writer);
catch
end
if isgraphics(f), close(f); end
end

function delete_if_exists(filename)
if exist(filename,'file'), delete(filename); end
end
