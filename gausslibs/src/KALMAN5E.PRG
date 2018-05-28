/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Compute the LM test with the procedure ML_derivatives
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

/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/

I1 = -H_cml;       /* First approximation of the Information matrix */
I2 = J_cml'J_cml;      /* Second approximation of the Information matrix */

output file = kalman5e.out reset;

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







