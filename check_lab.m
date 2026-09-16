function checks = check_lab(fullRun)
% Numerical and symmetry checks; these are not claims of experimental agreement.
if nargin==0, fullRun=false; end
cfg=teaching_config('smoke'); cfg.makePlots=false;
r=solve_airfoil_ns(cfg,0);
assert(all(isfinite(r.history(:))),'Nonfinite solution.');
assert(r.maxDiv<1e-7,'Projection failed to preserve mass.');
assert(r.pressureResidual<1e-9,'Pressure solve residual too large.');
assert(abs(r.CL)<1e-7,'Symmetric airfoil at zero incidence has nonzero lift.');
cfg2=cfg; cfg2.EI=2*cfg.EI; cfg2.GJ=2*cfg.GJ;
s1=wing_structure(cfg,1,0.2,0.1); s2=wing_structure(cfg2,1,0.2,0.1);
assert(abs(s2.tipBending/s1.tipBending-0.5)<1e-12);
assert(abs(s2.tipTwistDeg/s1.tipTwistDeg-0.5)<1e-12);
checks=struct('smokePassed',true,'fullSymmetryPassed',false);
if fullRun
    cfg=teaching_config('demo'); cfg.makePlots=false;
    p=solve_airfoil_ns(cfg,4); m=solve_airfoil_ns(cfg,-4);
    assert(abs(p.CL+m.CL)<1e-6,'Lift antisymmetry failed.');
    assert(abs(p.CD-m.CD)<1e-6,'Drag symmetry failed.');
    assert(p.CD>0 && p.CL>0,'Unexpected force signs.');
    checks.fullSymmetryPassed=true;
end
disp(checks);
end
