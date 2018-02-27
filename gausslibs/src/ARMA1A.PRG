/*
**  Estimation of a Vector ARMA(1,1) process
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

output file = arma1a.out reset;

t = hsec;
{beta,stderr,Mcov,LogL} = arma_ML(y,1,1,"HR2");
print;
print ftos((hsec-t)/100,"Computation time: %lf seconds",10,2);

output off;

save arma1 = beta;
save sigma1 = _arma_sigma;


