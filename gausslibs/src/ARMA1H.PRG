/*
**  Using vread to obtain the gradient and the Hessian
**  of the likelihood function under linear restrictions
**
**  Construction of the LR, Wald and LM tests
*/

new;
library optmum,tsm;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

_print = 0;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

/*
**  Unrestricted model
*/

load sv = arma1;

{beta1,stderr1,Mcov1,LogL1} = arma_ML(y,1,1,sv);

/*
**  4 restrictions:
**
**  AR1_21 = 0
**  AR1_12 = 0
**  MA1_21 = 0
**  MA1_12 = 0
*/

w = 1|0|0|2|3|0|0|4;
RR = design(w); r = zeros(8,1);

load sv = arma2; sv = sv[1 4 5 8];

{beta2,stderr2,Mcov2,LogL2} = arma_CML(y,1,1,sv,RR,r);

H = vread(_ml_derivatives,"H_matrix");
G = vread(_ml_derivatives,"G_matrix");

I = -H;       /* Approximation of the Information matrix */

output file = arma1h.out reset;

/*
**  Likelihood ratio statistic
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.06)
*/


LR = 2*(logl1-logl2);      /* Likelihood ratio                             */
pvalue = cdfchic(LR,4);    /* Approximating the noncentral chi-squared CDF */

print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/

LM = G'*inv(I)*G;
pvalue = cdfchic(LM,4);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


/*
**  Wald test
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.05)
*/

proc restrict(beta);
  local r;
  r = beta[2]|beta[3]|beta[6]|beta[7];
  retp(r);
endp;

Dr = gradp(&restrict,beta1);
r = restrict(beta1);

wald = r'*invpd(Dr*Mcov1*Dr')*r;
pvalue = cdfchic(wald,4);

print;
print ftos(wald,"Wald statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;

