
new;
library tsm,optmum;
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
save kal1c = a0;
save kal1d = P0;


call KFiltering(y,a0,P0);

at = KF_matrix(3);                 /* Read a(t) */
Pt = KF_matrix(4);                 /* Read P(t) */

output file = kalman1d.out reset;

j = 1;
do until j>10;
  print chrs(vals("=")*ones(30,1));
  print ftos(j,"Observation nø %lf",1,0); print;
  print ftos(j,"a(%lf)",1,0); print at[j,.]';
  print ftos(j,"P(%lf)",1,0); print xpnd(Pt[j,.]');

  print; call pause(2);
 
  j = j+1;
endo;

output off;
