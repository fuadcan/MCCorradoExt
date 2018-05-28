/*
**  Estimation of a Vector ARMA(1,1) process
**  Exact Maximum Likelihood
**
**  This programm is a modification of the example in the ARMA1B.prg file
**
*/

new;
library optmum,tsm;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;


/*  Build a procedure for the exact likelihood function */

proc ml(coeff);
  local beta,sigma,Pchol,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  beta = coeff[1:8];
  Pchol = (coeff[9]~0)|(coeff[10]~coeff[11]);
  SIGMA = Pchol*Pchol';

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);

  {a0,P0} = SSM_ic;
  if ismiss(a0);
    cls;
    print "Not implemented for non-stationnary ARMA models.";
    end;
  endif;

  call KFiltering(y,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

load ARMA1;
load SIGMA1;

sv = arma1|vech_(chol(SIGMA1));    /* Starting values */

output file = arma2e.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;

save EXACT_ML = theta;
