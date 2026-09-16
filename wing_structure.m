function s = wing_structure(cfg,Lprime,Dprime,Mprime)
% One-way load transfer: repeat the same 2-D section load along the span.
% z is spanwise, x is streamwise, and y is vertical. Positive Mprime raises the leading edge.
z = linspace(0,cfg.span,101)'; b = cfg.span;
s.z = z;
s.bending = Lprime*z.^2.*(6*b^2-4*b*z+z.^2)/(24*cfg.EI);
s.streamwiseBending = Dprime*z.^2.*(6*b^2-4*b*z+z.^2)/(24*cfg.EI);
s.twist = Mprime*z.*(2*b-z)/(2*cfg.GJ);
s.tipBending = s.bending(end); s.tipTwistDeg = s.twist(end)*180/pi;
s.rootLift = Lprime*b; s.rootDrag = Dprime*b;
s.rootBendingMoment = Lprime*b^2/2; s.rootTorque = Mprime*b;
% Positive theta raises the leading edge. This deformation is never fed back to CFD.
end
