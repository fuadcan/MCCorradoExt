
new;
library tsm,optmum,pgraph;
TSMset;

load y = kal1a;

proc ml(theta);
  local Z,d,H,T,c,R,Q,a0,P0,Logl;
  Z = (theta[1]~0)|(0~theta[2]);
  d = theta[3]|0;
  H = (theta[4]~0)|(0~theta[5]);
  T = (theta[6]~theta[7])|(0~theta[8]);
  c = theta[9 10];
  R = 1|1;
  Q = theta[11];
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(y,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;

let sv = 1 1 10 2 1 0.5 0.3 0.2 1 4 1;

_tsm_Mcov = 1;
_tsm_optmum = 1;

output file = kalman1h.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

print; print; call pause(2);

print "       True values    Estimated values ";
call printfmt(sv~theta,1~1);

output off;

save kal1e = theta;
