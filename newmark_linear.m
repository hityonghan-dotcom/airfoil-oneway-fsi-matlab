function [q,v,a] = newmark_linear(t,F,M,C,K,q0,v0)
% Linear structural dynamics: M*qdd + C*qd + K*q = F(t).
% Standard Newmark average-acceleration scheme: beta=1/4, gamma=1/2.
% M, C, and K are vectors of uncoupled diagonal modal coefficients.
if nargin<6, q0=zeros(size(M)); end
if nargin<7, v0=zeros(size(M)); end
N=numel(t); nd=numel(M);
assert(all(diff(t)>0) && isequal(size(F),[N nd]));
assert(all(M>0) && all(K>0) && all(C>=0));
q=zeros(N,nd); v=q; a=q; q(1,:)=q0; v(1,:)=v0;
a(1,:)=(F(1,:)-C.*v0-K.*q0)./M;
beta=0.25; gamma=0.5;
for n=1:N-1
    dt=t(n+1)-t(n);
    qp=q(n,:)+dt*v(n,:)+dt^2*(0.5-beta)*a(n,:);
    vp=v(n,:)+dt*(1-gamma)*a(n,:);
    a(n+1,:)=(F(n+1,:)-C.*vp-K.*qp)./(M+gamma*dt*C+beta*dt^2*K);
    q(n+1,:)=qp+beta*dt^2*a(n+1,:);
    v(n+1,:)=vp+gamma*dt*a(n+1,:);
end
end
