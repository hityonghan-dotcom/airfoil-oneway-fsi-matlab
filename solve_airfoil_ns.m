function r = solve_airfoil_ns(cfg,alphaDeg)
%SOLVE_AIRFOIL_NS Fixed-airfoil CFD, load extraction, and static load transfer.
%
% Read this function in five passes:
%   1. Lines below "GRID AND FIXED GEOMETRY": make the MAC grid and chi masks.
%   2. "INITIAL CONDITION": choose dt, eta, the pressure matrix, and the inlet state.
%   3. "TIME MARCH": predict momentum, apply Brinkman penalty, then project.
%   4. "AERODYNAMIC LOADS": turn the penalty reaction into L', D', and M'.
%   5. "REPORT": average the late-time window and package the result.
%
% The CFD geometry is fixed throughout this function.  The final call to
% WING_STRUCTURE uses the mean loads only; it never changes u, v, p, or chi.
%
% Adapted staggered-grid momentum discretization from MathWorks courseware.
% Copyright (c) 2025, The MathWorks, Inc. All rights reserved.
% See vendor/mathworks-cfd/LICENSE (BSD-3-Clause).
% Extensions: airfoil mask, implicit Brinkman constraint, variable-coefficient
% pressure projection, and force/moment output.

%% 1. GRID AND FIXED GEOMETRY
nx=cfg.nx; ny=cfg.ny; dx=cfg.Lx/nx; dy=cfg.Ly/ny;
assert(abs(dx-dy)<1e-12,'This case requires dx = dy; scale nx and ny together.');
assert(cfg.avgStart<cfg.endTime && cfg.avgStart>=0);
dt0=cfg.dtScale*min([0.15*dx/cfg.U,0.20*dx^2/cfg.nu]);
nt=ceil(cfg.endTime/dt0); dt=cfg.endTime/nt; eta=cfg.etaRatio*dt;
[xb,yb]=airfoil_geometry(cfg,alphaDeg);
% u is ny-by-(nx+1), v is (ny+1)-by-nx, and p is ny-by-nx.
% The arrays u and v below also carry a one-cell ghost layer for boundary conditions.
[xu,yu]=meshgrid(0:dx:cfg.Lx,dy/2:dy:cfg.Ly-dy/2);
[xv,yv]=meshgrid(dx/2:dx:cfg.Lx-dx/2,0:dy:cfg.Ly);
% chi is the volume fraction occupied by the fixed solid at each velocity face.
chiU=fraction(xu,yu,xb,yb,dx,dy,cfg.subcells);
chiV=fraction(xv,yv,xb,yb,dx,dy,cfg.subcells);
if isfield(cfg,'makeGridPreview') && cfg.makeGridPreview
    if isfield(cfg,'gridPreviewFolder')
        previewFolder=cfg.gridPreviewFolder;
    else
        previewFolder=fullfile(fileparts(mfilename('fullpath')),'results','mesh_preview');
    end
    if ~exist(previewFolder,'dir'), mkdir(previewFolder); end
    previewFile=fullfile(previewFolder,sprintf('grid_preview_alpha_%g.png',alphaDeg));
    plot_grid_preview(cfg,alphaDeg,xb,yb,chiU,chiV,previewFile);
end

%% 2. INITIAL CONDITION AND PRESSURE PROJECTION
% The implicit factor a=1/(1+dt*chi/eta) is the Brinkman no-slip step.
aU=1./(1+dt*chiU/eta); aV=1./(1+dt*chiV/eta);
[A,solver]=pressure_matrix(aU,aV,dx,dy,nx,ny);
u=cfg.U*ones(ny+2,nx+2); v=zeros(ny+2,nx+2);
% Start with zero velocity inside the solid; the first projection builds pressure consistency.
u(2:ny+1,2:nx+2)=cfg.U*(1-chiU);
[u,v]=outer_bc(u,v,cfg.U,nx,ny);
hist=zeros(nt,10); rho=cfg.rho; dA=dx*dy;
a=-alphaDeg*pi/180;
xe=cfg.center(1)+(cfg.elasticAxis-0.5)*cfg.c*cos(a);
ye=cfg.center(2)+(cfg.elasticAxis-0.5)*cfg.c*sin(a);
sumU=zeros(ny,nx+1); sumV=zeros(ny+1,nx); sumP=zeros(ny,nx); count=0;
record=isfield(cfg,'saveSnapshots') && cfg.saveSnapshots;
if record
    assert(cfg.snapshotDt>=dt,'snapshotDt must not be smaller than the actual time step.');
    capacity=ceil(cfg.endTime/cfg.snapshotDt)+2;
    snapshots.time=zeros(capacity,1);
    snapshots.U=zeros(ny,nx,capacity,'single');
    snapshots.V=zeros(ny,nx,capacity,'single');
    ns=0; nextSnapshot=0;
end

%% 3. TIME MARCH: MOMENTUM PREDICTION -> PENALTY -> PRESSURE PROJECTION
for n=1:nt
    % (a) Apply outer boundary conditions and remember the previous face velocities.
    [u,v]=outer_bc(u,v,cfg.U,nx,ny);
    uprev=u(2:ny+1,2:nx+2); vprev=v(2:ny+2,2:nx+1);

    % (b) Momentum prediction: convection and diffusion, before pressure and penalty.
    % Retain the upstream central interpolation, conservative convection, and diffusion stencil.
    [uh,vh]=intermediateVelocity(u,v,u,v,rho,rho*cfg.nu,nx,ny,dx,dy,dt);

    % (c) Implicit Brinkman step.  In a solid face chi=1, this drives velocity to u_s=0.
    us=aU.*uh(2:ny+1,2:nx+2); vs=aV.*vh(2:ny+2,2:nx+1);

    % (d) Projection: solve pressure, then subtract its staggered gradient.
    rhs=(diff(us,1,2)/dx+diff(vs,1,1)/dy)*rho/dt;
    pvec=solver\(-rhs(:)); pc=reshape(pvec,ny,nx);
    gpU=zeros(ny,nx+1); gpV=zeros(ny+1,nx);
    gpU(:,2:nx)=diff(pc,1,2)/dx;
    gpU(:,end)=-2*pc(:,end)/dx; % Outlet p=0 at a half-cell distance
    gpV(2:ny,:)=diff(pc,1,1)/dy;
    un=us-dt/rho*aU.*gpU; vn=vs-dt/rho*aV.*gpV;
    u(2:ny+1,2:nx+2)=un; v(2:ny+2,2:nx+1)=vn;

    %% 4. AERODYNAMIC LOADS FROM THE IMMERSED-BOUNDARY REACTION
    % f is force applied to the virtual fluid by the fixed body.
    % q is the equal-and-opposite body load, including virtual-fluid momentum correction.
    fx=-rho/eta*chiU.*un; fy=-rho/eta*chiV.*vn;
    qx=-fx+rho*chiU.*(un-uprev)/dt;
    qy=-fy+rho*chiV.*(vn-vprev)/dt;
    Dprime=sum(qx(:))*dA; Lprime=sum(qy(:))*dA;
    % About the spanwise axis, positive moment raises the leading edge.
    mu=(yu-ye).*qx; mv=-(xv-xe).*qy;
    Mprime=(sum(mu(:))+sum(mv(:)))*dA;
    div=diff(un,1,2)/dx+diff(vn,1,1)/dy;
    pr=norm(A*pvec+rhs(:),inf)/max(1,norm(rhs(:),inf));
    deepU=chiU>0.99; deepV=chiV>0.99;
    slip=max([abs(un(deepU));abs(vn(deepV));0])/cfg.U;
    maxspeed=max([abs(un(:));abs(vn(:))]);
    cfl=dt*(max(abs(un(:)))/dx+max(abs(vn(:)))/dy);
    if any(~isfinite([un(:);vn(:)])) || cfl>0.8
        error('Flow instability or CFL limit exceeded; reduce cfg.dtScale. Step %d.',n);
    end
    CL=Lprime/(0.5*rho*cfg.U^2*cfg.c);
    CD=Dprime/(0.5*rho*cfg.U^2*cfg.c);
    CM=Mprime/(0.5*rho*cfg.U^2*cfg.c^2);
    hist(n,:)=[n*dt,CL,CD,CM,max(abs(div(:)))*cfg.c/cfg.U,pr,slip,cfl,maxspeed,Lprime];
    if record && (n*dt>=nextSnapshot-1e-12 || n==nt)
        ns=ns+1; snapshots.time(ns)=n*dt;
        snapshots.U(:,:,ns)=single(0.5*(un(:,1:end-1)+un(:,2:end)));
        snapshots.V(:,:,ns)=single(0.5*(vn(1:end-1,:)+vn(2:end,:)));
        nextSnapshot=nextSnapshot+cfg.snapshotDt;
    end
    if n*dt>=cfg.avgStart
        sumU=sumU+un; sumV=sumV+vn; sumP=sumP+pc; count=count+1;
    end
    if mod(n,max(1,floor(nt/6)))==0
        fprintf('alpha=%g: t=%.2f/%g, CL=%.4f CD=%.4f div=%.2e\n', ...
            alphaDeg,n*dt,cfg.endTime,CL,CD,hist(n,5));
    end
end

%% 5. REPORT: LATE-TIME MEAN FLOW, COEFFICIENTS, AND ONE-WAY STATIC RESPONSE
select=hist(:,1)>=cfg.avgStart;
r.alphaDeg=alphaDeg; r.history=hist; r.dt=dt; r.eta=eta;
r.CL=mean(hist(select,2)); r.CD=mean(hist(select,3)); r.CM=mean(hist(select,4));
r.CLstd=std(hist(select,2)); r.CDstd=std(hist(select,3));
k=find(select); mid=floor(numel(k)/2);
r.windowChange=max(abs(mean(hist(k(1:max(1,mid)),2:3),1)- ...
    mean(hist(k(max(1,mid+1):end),2:3),1)));
r.maxDiv=max(hist(select,5)); r.maxSlip=max(hist(select,7));
r.pressureResidual=max(hist(select,6)); r.maxCFL=max(hist(:,8));
r.U=sumU/count; r.V=sumV/count; r.p=sumP/count;
r.x=dx/2:dx:cfg.Lx-dx/2; r.y=dy/2:dy:cfg.Ly-dy/2;
r.xb=xb; r.yb=yb; r.chiU=chiU; r.chiV=chiV;
r.Lprime=r.CL*0.5*rho*cfg.U^2*cfg.c;
r.Dprime=r.CD*0.5*rho*cfg.U^2*cfg.c;
r.Mprime=r.CM*0.5*rho*cfg.U^2*cfg.c^2;
r.structure=wing_structure(cfg,r.Lprime,r.Dprime,r.Mprime);
if record
    snapshots.time=snapshots.time(1:ns);
    snapshots.U=snapshots.U(:,:,1:ns); snapshots.V=snapshots.V(:,:,1:ns);
    r.snapshots=snapshots;
end
end

function chi=fraction(x,y,xb,yb,dx,dy,ns)
% Estimate the solid fraction at a velocity face by ns-by-ns point sampling.
% This makes the stair-step Cartesian body boundary less abrupt than a binary mask.
chi=zeros(size(x)); offsets=((1:ns)-0.5)/ns-0.5;
for i=offsets
    for j=offsets
        chi=chi+inpolygon(x+i*dx,y+j*dy,xb,yb)/ns^2;
    end
end
end

function [A,solver]=pressure_matrix(aU,aV,dx,dy,nx,ny)
% -div(a grad) is symmetric positive definite: inlet/top/bottom Neumann, outlet Dirichlet.
% Cell identifiers let the finite-volume stencil be assembled as a sparse matrix.
id=reshape(1:nx*ny,ny,nx); N=nx*ny;
i=id(:,1:end-1); j=id(:,2:end); w=aU(:,2:nx)/dx^2;
I=[i(:);j(:);i(:);j(:)]; J=[i(:);j(:);j(:);i(:)]; W=[w(:);w(:);-w(:);-w(:)];
i=id(1:end-1,:); j=id(2:end,:); w=aV(2:ny,:)/dy^2;
I=[I;i(:);j(:);i(:);j(:)]; J=[J;i(:);j(:);j(:);i(:)]; W=[W;w(:);w(:);-w(:);-w(:)];
i=id(:,end); w=2*aU(:,end)/dx^2;
A=sparse([I;i(:)],[J;i(:)],[W;w(:)],N,N);
solver=decomposition(A,'chol');
end

function [u,v]=outer_bc(u,v,U,nx,ny)
% Uniform inlet, extrapolated outlet, and zero normal velocity at top/bottom.
u(:,2)=U; u(:,1)=U;
% Outlet extrapolation supports the next momentum prediction; projection updates outlet velocity.
u(:,nx+2)=u(:,nx+1);
u(1,:)=2*U-u(2,:); u(ny+2,:)=2*U-u(ny+1,:);
v(2,:)=0; v(ny+2,:)=0; v(1,:)=-v(3,:);
v(:,1)=-v(:,2); v(:,nx+2)=v(:,nx+1);
end

% The following function is retained line-for-line from the upstream teaching source.
function [u,v] = intermediateVelocity(u,v,u_old,v_old,rho,mu,nx,ny,dx,dy,dt)
    for i=3:nx+1
        for j=2:ny+1
            % Interpolating velocities
            u_e = 0.5*(u_old(j,i) + u_old(j,i+1));
            u_w = 0.5*(u_old(j,i-1) + u_old(j,i));
            u_n = 0.5*(u_old(j,i) + u_old(j+1,i));
            u_s = 0.5*(u_old(j,i) + u_old(j-1,i));
            
            v_n = 0.5*(v_old(j+1,i-1) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i-1) + v_old(j,i));
            
            % Solving div(rho*u*u) and div(tau) 
            convection = -(rho*u_e*u_e - rho*u_w*u_w)/dx -(rho*v_n*u_n - rho*v_s*u_s)/dy;
            diffusion = mu*(u_old(j,i-1) - 2.0*u_old(j,i) + u_old(j,i+1))/dx/dx + mu*(u_old(j+1,i) - 2.0*u_old(j,i) + u_old(j-1,i))/dy/dy;
            
            % Calculate intermediate u velocity
            u(j,i) = rho*u_old(j,i) + dt*(diffusion + convection);
            u(j,i) = u(j,i)/rho;
        end
    end

    for i=2:nx+1
        for j=3:ny+1
            % Interpolating velocities
            v_e = 0.5*(v_old(j,i) + v_old(j,i+1));
            v_w = 0.5*(v_old(j,i) + v_old(j,i-1));
            v_n = 0.5*(v_old(j,i) + v_old(j+1,i));
            v_s = 0.5*(v_old(j,i) + v_old(j-1,i));

            u_e = 0.5*(u_old(j-1,i+1) + u_old(j,i+1));
            u_w = 0.5*(u_old(j-1,i) + u_old(j,i));

            % Solving div(rho*u*u) and div(tau)
            convection = -(rho*v_e*u_e - rho*v_w*u_w)/dx -(rho*v_n*v_n - rho*v_s*v_s)/dy;
            diffusion = mu*(v_old(j,i-1) - 2*v_old(j,i) + v_old(j,i+1))/dx/dx + mu*(v_old(j+1,i) - 2*v_old(j,i) + v_old(j-1,i))/dy/dy;
            
            % Calculate intermediate v velocity
            v(j,i) = rho*v_old(j,i) + dt*(diffusion + convection);
            v(j,i) = v(j,i)/rho;
        end
    end
end
