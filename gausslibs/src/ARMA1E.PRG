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
**  AR1_12 = 0
**  MA1_11 = 0
**  MA1_21 = 0
*/

w = 1|2|0|3|0|0|4|5;
RR = design(w); r = zeros(8,1);

output file = arma1e.out reset;

{beta,stderr,Mcov,LogL} = arma_CML(y,1,1,0,RR,r);

output off;





