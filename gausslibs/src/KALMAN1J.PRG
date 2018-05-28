
new;
library tsm,optmum,pgraph;
TSMset;

load y = kal1a;
load ESTIMATES = kal1e;

a0 = zeros(2,1);
P0 = zeros(2,2);


Z = (ESTIMATES[1]~0)|(0~ESTIMATES[2]);
d = ESTIMATES[3]|0;
H = (ESTIMATES[4]~0)|(0~ESTIMATES[5]);
T = (ESTIMATES[6]~ESTIMATES[7])|(0~ESTIMATES[8]);
c = ESTIMATES[9 10];
R = 1|1;
Q = ESTIMATES[11];

call SSM_build(Z,d,H,T,c,R,Q,0);
call KFiltering(y,a0,P0);

at = KF_matrix(3);
Pt = KF_matrix(4);

{as,Ps} = KSmoothing;

output file = kalman1j.out reset;

j = 1;
do until j>10;
  print chrs(vals("=")*ones(30,1));
  print ftos(j,"Observation nø %lf",1,0); print;
  print ftos(j,"Filter: a(%lf)",1,0); print at[j,.]';
  print ftos(j,"Smoothing: a(%lf|NOBS)",1,0); print as[j,.]';
  print ftos(j,"Filter: P(%lf)",1,0); print xpnd(Pt[j,.]');
  print ftos(j,"Smoothing: P(%lf|NOBS)",1,0); print xpnd(Ps[j,.]');

  print; call pause(2);
  
  j = j+1;
endo;

output off;


graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  _pnumht = 0.20; _paxht = 0.20; _ptitlht = 0.25;
  begwind;
  window(2,1,0);
  setwind(1);
    title("\201Estimates of \202a\201]t[");
    xlabel("Observations");
    _plegstr="Filter\0Smoothing"; _plegctl = { 2 5 1 1.25};
    xy(seqa(1,1,rows(y)),at[.,1]~as[.,1]);
  nextwind;
    title("\201Estimates of \202b\201]t[");
    xy(seqa(1,1,rows(y)),at[.,2]~as[.,2]);
  endwind;


