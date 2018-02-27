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

call arma_fevd(beta,2,0,SIGMA,8);

E1 ={}; E2 = {}; E3 = {};

Nr = 8;

i = 1;
do until i>Nr;
  E = varget("ERROR"$+ftos(i,"%lf",1,0));
  E1 = E1|E[1,.];
  E2 = E2|E[2,.];
  E3 = E3|E[3,.];
  i = i+1;
endo;

i = seqa(1,1,Nr);

output file = varx1d.out reset;

mask = 1~1~1~1;
let fmt[4,3]=  "*.*lf" 13 0 "*.*lf" 18 4 "*.*lf" 18 4 "*.*lf" 18 4;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in INVESTMENT (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
call printfm(i~(E1*100),mask,fmt);

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in INCOME (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
call printfm(i~(E2*100),mask,fmt);

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in CONSUMPTION (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment           Income      Consumption";
print "         periods      innovations      innovations      innovations";
print chrs(45*ones(79,1));
call printfm(i~(E3*100),mask,fmt);

output off;

