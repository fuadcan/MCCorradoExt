/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Testing Granger Causality
**  HARVEY, pages 35-43 and 93-97
**
**  Does not INCOME/CONSUMPTION cause INVESTMENT in the Granger sense ?
*/

new;
library tsm,optmum;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUMPT = y[.,3];
dcons = CONSUMPT-lag1(CONSUMPT);
data = dinv~dinc~dcons;

Y_ = data[4:76,.];

Z = eye(3)~zeros(3,3); d = zeros(3,1); H = zeros(3,3);
R = eye(3)|zeros(3,3);

a0 = data[2,.]'|data[3,.]'; P0 = zeros(6,6);

proc ml(beta);
  local const,PHI1,PHI2,Pchol;
  local T,c,Q,LogL;
  const = beta[19:21];
  PHI1 = reshape(beta[1:9],3,3)';
  PHI2 = reshape(beta[10:18],3,3)';
  Pchol = (beta[22]~0~0)|(beta[23]~beta[25]~0)|(beta[24]~beta[26]~beta[27]);
  T = (PHI1~PHI2)|(eye(3)~zeros(3,3));
  c = const|zeros(3,1);
  Q = Pchol*Pchol';
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(Y_,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;

output file = kalman5b.out reset;


_print = 0;
load sv = kal5a;

_tsm_Mcov = 2;
_tsm_optmum = 1;
{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

/*
**  Testing Granger Causality is equivalent to test:
**  theta[4] = theta[7] = theta[13] = theta[16] = 0
*/


CC = design(4|7|13|16)~zeros(4,11);

proc restrict(theta);
  retp(CC*theta);
endp;

Dr_ = gradp(&restrict,theta);
r_ = restrict(theta);

wald = r_'*invpd(Dr_*Mcov*Dr_')*r_;
pvalue = cdfchic(wald,4);

print;
print ftos(wald,"Wald statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


output off;





