new;
library tsm,optmum,pgraph;
TSMset;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

{lambda,Iy} = PDGM2(y);

/*
** Procedure to compute the sgf of the SSM
*/

proc sgf(theta,lambda);
  local beta,Pchol,SIGMA,T,Q,H,Z,d,c,R;
  local Gy;

  beta = theta[1:8];
  Pchol = (theta[9]~0)|(theta[10]~theta[11]);
  SIGMA = Pchol*Pchol';

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  Gy = sgf_SSM(lambda);

  retp(Gy);
endp;


proc ml(theta);
  local Gy,Nstar,logl,j,Ij,Gj;

  Gy = sgf(theta,lambda);
  Nstar = rows(lambda);

  logl = zeros(Nstar,1);

  j = 1;
  do until j > Nstar;
   Ij = xpnd2(Iy,j);
   Gj = xpnd2(Gy,j);
   logl[j] = -0.5*ln(det(Gj)) -0.5*sumc(diag((inv(Gj)*Ij)));
   j = j + 1;
  endo;

  logl = real(logl);

  retp(logl);
endp;

load ARMA1;
load SIGMA1;

sv = arma1|vech_(chol(SIGMA1));    /* Starting values */

output file = fdml6.out reset;

_tsm_Mcov = 1;
{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;

