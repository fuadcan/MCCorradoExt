/*
**  Estimation of a Vector ARMA(1,1) process
**  Exact Maximum Likelihood
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

/* First, perform a Conditional Maximum Likelihood estimation */

{beta,stderr,Mcov,LogL} = arma_ML(y,1,1,"hr1");
SIGMA = _arma_SIGMA;


/* Second, convert the ARMA process to a state-space model */

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);


/* Third, set the starting values for Exact Maximum Likelihood */

Pchol = chol(sigma');
coeff = beta|pchol[1,1]|pchol[2,1]|pchol[2,2];


/* Fourth, we use the CML estimates as initial conditions */

{a0,P0} = SSM_ic;


/* Finally, build a procedure for the exact likelihood */

proc ml(coeff);
  local beta,sigma,Pchol,Z,d,H,T,c,R,Q;
  local Logl;

  beta = coeff[1:8];
  Pchol = (coeff[9]~0)|(coeff[10]~coeff[11]);
  SIGMA = Pchol*Pchol';

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(y,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 0;      /* BHHH algorithm */
_tsm_gtol = 0.001;

output file = arma1b.out reset;

{beta,stderr,Mcov,Logl} = TD_ml(&ml,Coeff);

output off;



