/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, section 13.7, pages 471-478
**
**  DAVIDSON and MACKINNON [1990], Specification tests based on artificial
**  regressions, Journal of the American Statistical Association, 85, 220-227
**
**  Compute the LM test with an OPG artificial regression
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


/*
**  Test the following restrictions:
**
**  Theta[2] = 0
**  Theta[3] = 0
**  Theta[4] = 0
**  Theta[5] = 0
**  Theta[7] = 0
**  Theta[10] = 0
**  Theta[11] = 0
**  Theta[12] = 0
**  Theta[13] = 0
**  Theta[14] = 0
**  Theta[16] = 0
**  Theta[17] = 0
**  Theta[18] = 0
**  Theta[19] = 0
**  Theta[23] = 0
*/


load theta = kal5b;

{J_cml,G_cml,H_cml,I_cml} = TDml_derivatives(&ml,theta);

Nobs = rows(J_cml);
ll = ones(Nobs,1);

beta = ll/J_cml;                           /* OPG artificial regression */
resid = ll-J_cml*beta;
R2 = 1 - (resid'resid)/(ll'll);            /* R2 statistic              */
LM = Nobs*R2;                              /* n * R2                    */
pvalue = cdfchic(LM,15);

output file = kalman5f.out reset;

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;








