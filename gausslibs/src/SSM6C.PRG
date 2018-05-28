/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Forecast Error Variance Decomposition
**  See LUTKEPOHL, section 3.7.3
**
*/

new;
library tsm,optmum;

cls;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
data = dinv~dinc~dcons;

_print = 0;
{theta,stderr,Mcov,LOGL} = varx_LS(data,1,2);  /* Constant and 2 lags */

beta = theta[1:18];  /* Estimates of the AR part */
SIGMA = _varx_SIGMA;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

Omega = SSM_fevd(8);

i = seqa(1,1,8);

output file = ssm6c.out reset;

mask = 1~1~1~1;
let fmt[4,3]=  "*.*lf" 13 0 "*.*lf" 18 4 "*.*lf" 18 4 "*.*lf" 18 4;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in INVESTMENT (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
i = 1;
do until i > 8;
  rsp = xpnd2(Omega,i);
  rsp = rsp[1,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in INCOME (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
i = 1;
do until i > 8;
  rsp = xpnd2(omega,i);
  rsp = rsp[2,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in CONSUMPTION (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
i = 1;
do until i > 8;
  rsp = xpnd2(Omega,i);
  rsp = rsp[3,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

output off;



