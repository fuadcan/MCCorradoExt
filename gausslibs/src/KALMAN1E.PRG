
new;
library tsm,optmum;
TSMset;

load y = kal1a;
load a = kal1b;

/* Time-variant SSM */

proc Z(t);
  local Zt;
  Zt = eye(2);
  retp(Zt);
endp;

proc d(t);
  local dt;
  dt = 10|0;
  retp(dt);
endp;

proc H(t);
  local Ht;
  Ht = {2 0, 0 1};
  retp(Ht);
endp;

proc T(t);
  local Tt;
  Tt = {0.5 0.3, 0 0.2};
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = 1|4;
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = 1|1;
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = 1;
  retp(Qt);
endp;

call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);   /* cn = 1 */

load a0 = kal1c;
load P0 = kal1d;

call KFiltering(y,a0,P0);

at = KF_matrix(3);                 /* Read a(t) */
Pt = KF_matrix(4);                 /* Read P(t) */

output file = kalman1e.out reset;

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

print;
print "Now, compare the results from kalman1d.out and kalman1e.out...";

