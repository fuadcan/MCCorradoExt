/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Estimation of a VAR process by Kalman filter
**  (Harvey, pages 70-73, 82-85 and chapter 13)
**
*/

new;
library tsm,optmum;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
data = dinv~dinc~dcons;

Y_ = data[2:76,.];

Z = eye(3)~zeros(3,3); d = zeros(3,1); H = zeros(3,3);
R = eye(3)|zeros(3,3);

proc ml(beta);
  local const,PHI1,PHI2,Pchol;
  local T,c,Q,a0,P0,LogL;
  const = beta[19:21];
  PHI1 = reshape(beta[1:9],3,3)';
  PHI2 = reshape(beta[10:18],3,3)';
  Pchol = (beta[22]~0~0)|(beta[23]~beta[25]~0)|(beta[24]~beta[26]~beta[27]);
  T = (PHI1~PHI2)|(eye(3)~zeros(3,3));
  c = const|zeros(3,1);
  Q = Pchol*Pchol';
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(Y_,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;

/*  Starting values */

_print = 0;

{theta1,stderr,Mcov,LogL} = varx_ML(data,1,2);

output file = kalman5a.out reset;

_print = 1;
_tsm_Mcov = 2;
_tsm_optmum = 1;
_tsm_parnm =
              "PHI1_11"|"PHI1_21"|"PHI1_31"|
              "PHI1_12"|"PHI1_22"|"PHI1_32"|
              "PHI1_13"|"PHI1_23"|"PHI1_33"|
              "PHI2_11"|"PHI2_21"|"PHI2_31"|
              "PHI2_12"|"PHI2_22"|"PHI2_32"|
              "PHI2_13"|"PHI2_23"|"PHI2_33"|
              "mu1"|"mu2"|"mu3"|
              "Pchol1"|"Pchol2"|"Pchol3"|"Pchol4"|"Pchol5"|"Pchol6";

sv = theta1;

{theta2,stderr,Mcov,Logl} = TD_ml(&ml,sv);
save kal5a = theta2;

output off;



