/*
**  Estimation of a Vector ARMA(1,1) process under linear restrictions
**  Conditional Maximum Likelihood
*/

new;
library optmum,tsm;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

/*
**  AR1_11 = 1
*/

RR = design(0|seqa(1,1,7)); r = 1|zeros(7,1);

output file = arma1c.out reset;

{beta,stderr,Mcov,LogL} = arma_CML(y,1,1,0,RR,r);

/* You may use the unconstrained parameters to compute the starting values */

print;
print "You may use the unconstrained parameters to compute "\
         "the starting values."; print; call pause(2);

load arma1;
sv = inv(RR'RR)*RR'(arma1-r);

call arma_CML(y,1,1,sv,RR,r);

output off;


