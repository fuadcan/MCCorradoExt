
new;
library tsm,optmum,pgraph;
TSMset;

load y = kal1a;
load a = kal1b;

/* Time-invariant SSM */

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */

{a0,P0} = SSM_ic;                  /* Initial conditions */

call KFiltering(y,a0,P0);

at = KF_matrix(3);                 /* Read a(t) */
Pt = KF_matrix(4);                 /* Read P(t) */

alpha = at[.,1];

bound = alpha + (-1.96~1.96).*sqrt(Pt[.,1]);

graphset;
  fonts("simplex simgrma"); _pdate = ""; _pnum = 2;
  _pnumht = 0.20; _paxht = 0.20; _ptitlht = 0.25;
  begwind;
  window(2,1,0);
  setwind(1);
    title("\201True values of \202a\201]t[");
    _pltype = 6;
    xlabel("Observation");
    xy(seqa(1,1,100),a[.,1]);
  nextwind;
    title("\201Estimated values (with the confidence interval) of"\
          " \202a\201]t[");
    _pltype = 6|1|1;
    _pcolor = 10|11|11;
    xy(seqa(1,1,100),alpha~bound);
  endwind;


