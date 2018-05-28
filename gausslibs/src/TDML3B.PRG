/*
**  Estimation of a linear model under restrictions
*/

new;
library tsm,optmum;

declare external epsilon_;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
INCOME = y[.,2];
CONSUMP = y[.,3];

X = lag1(income)~lag1(consump)~ones(76,1)~lag1(invest);
Y = INCOME;

{data,retcode} = MISSING(Y~X,0);

Y_ = data[.,1];
X_ = data[.,2:5];

proc ml(theta);
  local beta,sigma,epsilon,logl;

  beta = theta[1:4];
  sigma = theta[5];
  epsilon = Y_-X_*beta;
  epsilon_ = epsilon;

  logl = -0.5*ln(2*pi)-0.5*ln(sigma^2)-0.5*(epsilon^2)/(sigma^2);
  retp(logl);
endp;

output file = tdml3b.out reset;

_print = 1;
_tsm_optmum = 1;
_tsm_Mcov = 2;
_tsm_parnm  = "INC(t-1)"|"CON(t-1)"|"CONSTANT"|"INV(t-1)"|"sigma";
sv = rndu(4,1);

/*
** The INC(t-1) coefficient is set equal to 1
*/

r = 1|zeros(4,1);
RR = design(0|1|2|3|4);

{theta,stderr,Mcov,Logl} = TD_cml(&ml,sv,RR,r);

output off;



