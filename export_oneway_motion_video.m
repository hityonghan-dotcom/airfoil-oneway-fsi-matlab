function info = export_oneway_motion_video(r,cfg,filename)
% Export a four-panel one-way-FSI teaching video.
% Top left: frozen fixed-geometry CFD field. Top right: enlarged structural motion.
% Bottom panels: physical structural response and CFD force coefficients.
assert(isfield(r,'snapshots') && isfield(r,'dynamics'));
d=r.dynamics; fps=20; frameDt=0.025; times=(0:frameDt:cfg.endTime)';
writer=VideoWriter(filename,'MPEG-4'); writer.FrameRate=fps; writer.Quality=95; open(writer);
f=figure('Visible','off','Color','white','Position',[50 50 1920 1080]);
if isprop(f,'Theme'), f.Theme='light'; end
cleaner=onCleanup(@() close_resources(writer,f));
C=palette(); tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
title(tl,'One-way FSI: fixed-geometry CFD loads drive a visibly moving structural model', ...
    'FontWeight','bold');
[X,~]=meshgrid(r.x,r.y);
angle=-r.alphaDeg*pi/180;
xe=cfg.center(1)+(cfg.elasticAxis-0.5)*cfg.c*cos(angle);
ye=cfg.center(2)+(cfg.elasticAxis-0.5)*cfg.c*sin(angle);
x0=r.xb-xe; y0=r.yb-ye;

% Top-left: the original CFD field. The black body removes virtual-solid pixels from view.
axFlow=nexttile(tl); field=imagesc(axFlow,r.x,r.y,zeros(size(X))); axis(axFlow,'xy'); hold(axFlow,'on');
fill(axFlow,r.xb,r.yb,C.fixed,'EdgeColor',C.fixed);
flowReference=plot(axFlow,r.xb,r.yb,'--','Color',[0.85 0.85 0.85],'LineWidth',0.9);
flowMoving=fill(axFlow,r.xb,r.yb,C.moving,'FaceAlpha',0.70,'EdgeColor',[0.1 0.1 0.1],'LineWidth',0.8);
axis(axFlow,'equal'); xlim(axFlow,[0 cfg.Lx]); ylim(axFlow,[0 cfg.Ly]); colormap(axFlow,turbo); clim(axFlow,[0 1.5]);
cb=colorbar(axFlow); cb.Label.String='Speed magnitude, |u| / U_inf';
style_axes(axFlow); xlabel(axFlow,'x (m)'); ylabel(axFlow,'y (m)');
legend(axFlow,[flowReference flowMoving],{'Fixed CFD geometry','Moving structural overlay'}, ...
    'Location','southoutside','Box','on','AutoUpdate','off');

% Top-right: a local motion view. It uses a larger display scale than the field overlay.
axZoom=nexttile(tl); hold(axZoom,'on');
zoomReference=plot(axZoom,x0,y0,'--','Color',C.gray,'LineWidth',1.4);
zoomMoving=fill(axZoom,x0,y0,C.moving,'FaceAlpha',0.80,'EdgeColor',[0.08 0.08 0.08],'LineWidth',1.2);
plot(axZoom,0,0,'+','Color',C.gray,'MarkerSize',9,'LineWidth',1.2);
axis(axZoom,'equal'); xlim(axZoom,[-0.65 0.65]*cfg.c); ylim(axZoom,[-0.42 0.42]*cfg.c);
style_axes(axZoom); xlabel(axZoom,'x relative to elastic axis (m)'); ylabel(axZoom,'y relative to elastic axis (m)');
legend(axZoom,[zoomReference zoomMoving],{'Undeformed reference','Moving tip section'}, ...
    'Location','southoutside','Box','on','AutoUpdate','off');

% Bottom-left: actual, unamplified structural quantities.
axResponse=nexttile(tl); yyaxis(axResponse,'left');
hDisplacement=plot(axResponse,NaN,NaN,'Color',C.blue,'LineWidth',1.5); ylabel(axResponse,'Actual tip displacement (um)');
ylim(axResponse,range_pad(d.h*1e6));
yyaxis(axResponse,'right'); hTwist=plot(axResponse,NaN,NaN,'Color',C.orange,'LineWidth',1.5); ylabel(axResponse,'Actual tip twist (deg)');
ylim(axResponse,range_pad(d.theta*180/pi)); xlim(axResponse,[0 cfg.endTime]);
style_axes(axResponse); xlabel(axResponse,'Physical time (s)');
title(axResponse,sprintf('Structural response: bending %.2f Hz, torsion %.2f Hz',d.naturalHz));

% Bottom-right: unchanged CFD forcing used by the one-way structural solve.
axLoads=nexttile(tl); hold(axLoads,'on');
hLift=plot(axLoads,NaN,NaN,'Color',C.blue,'LineWidth',1.5);
hDrag=plot(axLoads,NaN,NaN,'Color',C.orange,'LineWidth',1.5);
steady=r.history(:,1)>=0.1*cfg.c/cfg.U;
ylim(axLoads,range_pad(r.history(steady,2:3))); xlim(axLoads,[0 cfg.endTime]);
legend(axLoads,[hLift hDrag],{'C_L','C_D'},'Location','northeast','Box','off','AutoUpdate','off');
style_axes(axLoads); xlabel(axLoads,'Physical time (s)'); ylabel(axLoads,'Force coefficient');
title(axLoads,'CFD loads passed to the structural model');

fieldScale=cfg.motionDisplayScale;
zoomScale=[300 60]; % Deliberate visual-only zoom: the plotted curves remain physical.
temp=[tempname '.png']; tempCleaner=onCleanup(@() delete_temp(temp));
motion=interp1(d.time,[d.h d.theta],times,'linear'); snapshotIndex=0;
for k=1:numel(times)
    t=times(k);
    j=find(r.snapshots.time<=t+1e-12,1,'last');
    if isempty(j)
        if snapshotIndex~=0, snapshotIndex=0; end
        if k==1, field.CData=ones(size(X)); end
        fluidTime=0;
    else
        if j~=snapshotIndex
            field.CData=hypot(double(r.snapshots.U(:,:,j)),double(r.snapshots.V(:,:,j)))/cfg.U;
            snapshotIndex=j;
        end
        fluidTime=r.snapshots.time(j);
    end
    h=motion(k,1); theta=motion(k,2);
    update_airfoil(flowMoving,xe,ye,x0,y0,h,theta,fieldScale,false);
    update_airfoil(zoomMoving,0,0,x0,y0,h,theta,zoomScale,true);
    title(axFlow,sprintf('Frozen CFD: alpha = %g deg, snapshot t = %.3f s | overlay t = %.3f s', ...
        r.alphaDeg,fluidTime,t));
    title(axZoom,sprintf('Enlarged structural motion: displacement x%g, twist x%g',zoomScale));
    selected=d.time<=t+1e-12;
    hDisplacement.XData=d.time(selected); hDisplacement.YData=1e6*d.h(selected);
    hTwist.XData=d.time(selected); hTwist.YData=d.theta(selected)*180/pi;
    selected=r.history(:,1)<=t+1e-12;
    hLift.XData=r.history(selected,1); hLift.YData=r.history(selected,2);
    hDrag.XData=r.history(selected,1); hDrag.YData=r.history(selected,3);
    print(f,temp,'-dpng','-r110'); writeVideo(writer,imread(temp));
    if mod(k,100)==0, fprintf('Motion video frame %d/%d\n',k,numel(times)); end
end
close(writer); close(f);
info=struct('filename',filename,'frames',numel(times),'fps',fps, ...
    'playbackSeconds',numel(times)/fps,'physicalTimeEnd',times(end), ...
    'fieldDisplayScale',fieldScale,'zoomDisplayScale',zoomScale);
end

function update_airfoil(handle,xe,ye,x0,y0,h,theta,scale,relative)
rotation=-scale(2)*theta;
x=xe+cos(rotation)*x0-sin(rotation)*y0;
y=ye+sin(rotation)*x0+cos(rotation)*y0+scale(1)*h;
if relative
    x=x-xe; y=y-ye;
end
handle.XData=x; handle.YData=y;
end

function style_axes(ax)
grid(ax,'on'); box(ax,'on'); ax.FontName='Arial'; ax.FontSize=10; ax.LineWidth=0.8;
ax.GridColor=[0.85 0.87 0.90]; ax.GridAlpha=0.7;
end

function lim=range_pad(values)
lo=min(values(:)); hi=max(values(:)); pad=max(1e-9,0.12*(hi-lo));
lim=[min(0,lo-pad) max(0,hi+pad)];
end

function C=palette()
C.blue=[0.00 0.35 0.70]; C.orange=[0.90 0.38 0.05];
C.gray=[0.35 0.38 0.42]; C.fixed=[0.08 0.12 0.18]; C.moving=[0.92 0.28 0.05];
end

function close_resources(writer,f)
try
    close(writer);
catch
end
if isgraphics(f), close(f); end
end

function delete_temp(filename)
if exist(filename,'file'), delete(filename); end
end
