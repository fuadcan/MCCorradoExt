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

{beta,stderr,Mcov,LogL} = arma_ML(y,5,0,"hr1");
SIGMA = _arma_SIGMA;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,5,0,SIGMA);
call SSM(Z,d,H,T,c,R,Q,0);


