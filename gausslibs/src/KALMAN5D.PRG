/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  LM and LR tests
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


/*
**  The Kal5a.fmt matrix corresponds to the unrestricted estimates
**  while the kal5b.fmt matrix corresponds to the restricted estimates
*/


load theta1 = kal5a;
load theta2 = kal5b;

w = 1|0|0|0|0|2|0|3|4|0|0|0|0|0|5|0|0|0|0|6|7|8|0|9|10|11|12;
RR_ = design(w); r_ = zeros(27,1);

_print = 0;
_tsm_Mcov = 1;
_tsm_optmum = 1;

sv = theta2[1 6 8 9 15 20 21 22 24 25 26 27];

call TD_cml(&ml,sv,RR_,r_);

output file = kalman5d.out reset;

/*
**  Likelihood ratio statistic
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.06)
*/

Logl1 = sumc(ml(theta1));
Logl2 = sumc(ml(theta2));

LR = 2*(logl1-logl2);      /* Likelihood ratio                              */
pvalue = cdfchic(LR,15);    /* Approximating the noncentral chi-squared CDF */

print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/

J_cml = vread(_ml_derivatives,"J_matrix");
H_cml = vread(_ml_derivatives,"H_matrix");
G_cml = vread(_ml_derivatives,"G_matrix");

I1 = -H_cml;       /* First approximation of the Information matrix */
I2 = J_cml'J_cml;      /* Second approximation of the Information matrix */

LM1 = G_cml'*inv(I1)*G_cml;
pvalue = cdfchic(LM1,15);

print;
print ftos(LM1,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


LM2 = G_cml'*inv(I2)*G_cml;
pvalue = cdfchic(LM2,15);

print;
print ftos(LM2,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;






