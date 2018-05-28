
new;
library tsm,optmum,pgraph;
TSMset;

load y = kal1a;

/* Time-invariant SSM */

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */

load a0 = kal1c;
load P0 = kal1d;

call KFiltering(y,a0,P0);
Logl = KF_ml;

st = ftos(sumc(Logl),"Value of the Log-Likelihood function: %lf",3,4);

graphset;
  _pdate = ""; _pnum = 2;
  title(st);
  xlabel("Observation");
  ylabel("Log-likelihood");
  xy(seqa(1,1,100),Logl);

