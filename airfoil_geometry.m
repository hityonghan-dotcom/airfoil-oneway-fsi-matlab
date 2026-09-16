function [xb,yb] = airfoil_geometry(cfg,alphaDeg)
% Closed-trailing-edge NACA 00xx. Positive incidence raises the leading edge.
s = linspace(0,pi,201)'; x = (1-cos(s))/2;
yt = 5*cfg.thickness*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2 ...
    +0.2843*x.^3-0.1036*x.^4);
X = [flipud(x);x(2:end)]*cfg.c-cfg.c/2;
Y = [flipud(yt);-yt(2:end)]*cfg.c;
a = -alphaDeg*pi/180;
xb = cfg.center(1)+cos(a)*X-sin(a)*Y;
yb = cfg.center(2)+sin(a)*X+cos(a)*Y;
end
