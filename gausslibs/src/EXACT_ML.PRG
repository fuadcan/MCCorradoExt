/*
**  Estimation of a vector ARMA process
**  Exact Maximum Likelihood
**
*/

new;
library optmum,tsm;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);

K = 2;
p_ar = 1;
q_ma = 1;
Y = di~invent;

/*
** K is the dimension of the vector ARMA
** p_ar is the autoregressive order
** q_ma is the moving average order
** y is the matrix Nobs*K of data
*/


proc ml(theta);
  local beta,v,Pchol,SIGMA,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  /* AR and MA parameters */

  beta = theta[1:K^2*(p_ar+q_ma)];

  /* vech_(chol(SIGMA)) parameters */

  v = theta[1+K^2*(p_ar+q_ma):K^2*(p_ar+q_ma)+K*(K+1)/2];

  /* chol(SIGMA) */

  Pchol = xpnd_(v).*lowmat(ones(K,K));

  /* SIGMA */

  SIGMA = Pchol*Pchol';

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,p_ar,q_ma,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);

  {a0,P0} = SSM_ic;
  if ismiss(a0);
    cls;
    print "Not implemented for Non stationnary ARMA Models.";
    end;
  endif;

  call KFiltering(y,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 1;
_tsm_gtol = 0.001;

_tsm_Mcov = 0;

{sv,stderr,Mcov,Logl} = varx_ls(y,0,1);

sv = sv|zeros(4,1)|vech_(chol(_varx_SIGMA)');

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);






