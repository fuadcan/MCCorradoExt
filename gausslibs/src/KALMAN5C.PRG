/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Estimation of a VAR process using
**  the Kalman filter under linear restrictions
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
**  We impose the following linear restrictions:
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
w = 1|0|0|0|0|2|0|3|4|0|0|0|0|0|5|0|0|0|0|6|7|8|0|9|10|11|12;
RR_ = design(w); r_ = zeros(27,1);

_print = 1;
_tsm_Mcov = 1;
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

load sv = kal5a;
sv = inv(RR_'RR_)*RR_'(sv-r_);

output file = kalman5c.out reset;

{theta2,stderr,Mcov,Logl2} = TD_cml(&ml,sv,RR_,r_);

output off;

save kal5b = theta2;




