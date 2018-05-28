/*
** Estimation of an univariate ARMA process with
** the arma_Filter procedure
*/

new;
library tsm,optmum;
TSMset;

Nobs = 200;

y = RND_arma(0.6|0.3|0.25,2,1,1,miss(0,0),Nobs);

proc ml(beta);
  local e,sigma2,Logl;
  e = arma_filter(y,beta[1:3],2,1);
  e = missrv(e,0);
  sigma2 = beta[4]^2;
  logL = -0.5*ln(2*pi)-0.5*ln(sigma2)-0.5*(e^2)/sigma2;
  retp(logL);
endp;

output file = filter1a.out reset;

_tsm_Mcov = 1;
{beta1,stderr1,Mcov1,LogL1} = TD_ml(&ml,0.6|0.3|0.25|1);

_tsm_optmum = 0;
{beta2,stderr,Mcov,LogL2} = arma_ML(y,2,1,beta1[1:3]);

output off;
