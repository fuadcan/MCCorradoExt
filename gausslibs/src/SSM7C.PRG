/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Forecast Error Variance Decomposition
*/

new;
library tsm,optmum,pgraph;

cls;

load beta = arma1;
load sigma = sigma1;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

Omega = SSM_fevd(8);

output file = ssm7c.out reset;

mask = 1~1~1;
let fmt[3,3]=  "*.*lf" 13 0 "*.*lf" 18 4 "*.*lf" 18 4;

print "FORECAST ERROR VARIANCE DECOMPOSITION in Investment (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment      Business Inventories  ";
print "         periods      innovations      innovations   ";
print chrs(45*ones(79,1));
i = 1;
do until i > 8;
  rsp = xpnd2(Omega,i);
  rsp = rsp[1,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in Business Inventories (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment      Business Inventories  ";
print "         periods      innovations      innovations   ";
print chrs(45*ones(79,1));
i = 1;
do until i > 8;
  rsp = xpnd2(Omega,i);
  rsp = rsp[2,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

output off;




