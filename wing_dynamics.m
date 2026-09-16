function d = wing_dynamics(cfg,r)
% One-way flow-to-structure transfer: read fixed-geometry loads only.
% One bending mode and one torsion mode use the uniform-load static shapes.
b=cfg.span; rho=cfg.rho; U=cfg.U; c=cfg.c;
d.time=[0;r.history(:,1)];
d.Lprime=[0;r.history(:,2)]*(0.5*rho*U^2*c);
d.Mprime=[0;r.history(:,4)]*(0.5*rho*U^2*c^2);
% phi(s)=(6s^2-4s^3+s^4)/3; psi(s)=2s-s^2; s=z/b.
M=[cfg.massPerLength*b*104/405, cfg.polarInertiaPerLength*b*8/15];
K=[16*cfg.EI/(5*b^3),4*cfg.GJ/(3*b)];
C=2*cfg.dampingRatio.*sqrt(M.*K);
F=[(2*b/5)*d.Lprime,(2*b/3)*d.Mprime];
[q,v,a]=newmark_linear(d.time,F,M,C,K);
d.h=q(:,1); d.theta=q(:,2); d.velocity=v; d.acceleration=a;
d.M=M; d.C=C; d.K=K; d.generalizedForce=F;
d.naturalHz=sqrt(K./M)/(2*pi);
d.quasiStatic=F./K;
z=linspace(0,b,101)'; s=z/b;
d.z=z; d.phi=(6*s.^2-4*s.^3+s.^4)/3; d.psi=2*s-s.^2;
balance=a.*M+v.*C+q.*K-F;
d.maxBalanceResidual=max(abs(balance(:)))/max(1,max(abs(F(:))));
end
