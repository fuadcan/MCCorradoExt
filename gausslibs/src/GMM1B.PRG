new;
library tsm,optmum;
TSMset;


Nobs = 1000;
x = 10*rndu(Nobs,4);
beta = rndn(4,1);
sigma = 2;
y = x*beta + rndn(Nobs,1)*sigma;

sv = beta|sigma; /* starting values */

proc ml(theta);
  local beta,sigma,u,logl;

  beta = theta[1:4];
  sigma = theta[5];

  u = y - x*beta;

  logl = -0.5*ln(2*pi)-0.5*ln(sigma^2)-0.5*u.*u/(sigma^2);

  retp(logl);
endp;


proc H(theta);
  local beta,sigma,h1,M,i;

  M = zeros(Nobs,6);

  beta = theta[1:4];
  sigma = theta[5];

  /* first moment */
  h1 = y - x*beta;
  M[.,1] = h1;
  /* second moment */
  M[.,2] = h1.*h1 - sigma^2;

  i = 1;
  do until i > 4;
    M[.,2+i] = h1.*x[.,i];
    i = i + 1;
  endo;

  retp(M);
endp;

proc Jproc(theta);
  local beta,sigma,D,i;

  D = zeros(5,6);

  beta = theta[1:4];
  sigma = theta[5];

  D[.,1] = meanc(-x)|0;
  D[.,2] = meanc(-2.*x.*(y-x*beta))|-2*sigma;

  i = 1;
  do until i > 4;
    D[.,2+i] = meanc(-x.*x[.,i])|0;
    i = i + 1;
  endo;

  D = D'; /* Jacobian matrix 6*5 */

  retp(D);
endp;

output file = gmm1b.out reset;

/*
** beta[1] = beta[2];
*/

RR = design(1|1|2|3|4); r = zeros(5,1);

_gmm_RR = RR;
_gmm_r = r;

sv = invpd(RR'RR)*RR'(sv-r);

_tsm_parnm = "beta1"|"beta2"|"beta3"|"beta4"|"sigma";

_gmm_lags = 2;
_gmm_Jacobian_proc = &Jproc;
{theta1,stderr,Mcov,Qmin} = gmm(&h,sv);

{theta2,stderr,Mcov,Logl} = TD_cml(&ml,sv,RR,r);

output off;


