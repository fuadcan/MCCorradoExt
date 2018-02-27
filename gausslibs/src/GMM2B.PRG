new;
library tsm,optmum;

Nobs = 250;
alpha0 = 0.5; alpha1 = 0.6;

u = zeros(Nobs,1);
h = zeros(Nobs,1);
i = 2;
do until i > Nobs;
  h[i] = sqrt(alpha0^2 + alpha1^2*u[i-1]^2);
  u[i] = rndn(1,1)*h[i];
  i = i+1;
endo;

x =  3*rndu(Nobs,1);
y = 2 + 3*x + u;


proc ml(theta);
  local u,u2,h2,logl;

  u = y-theta[1]-theta[2]*x;

  u2 = u.*u;
  h2 = theta[3]^2 + theta[4]^2*lag1(u2);
  h2[1] = theta[3]^2;

  logl = -0.5*ln(2*pi)-0.5*ln(h2)-0.5*u2./h2;

  retp(logl);
endp;


proc ARCH(theta);
  local M,u,h2,u2,i;

  M = zeros(Nobs,4);

  u = y-theta[1]-theta[2]*x;

  /* first moment */

  M[.,1] = u;

  /* computing h(t)^2 */

  u2 = u.*u;
  h2 = theta[3]^2 + theta[4]^2*lag1(u2);
  h2[1] = theta[3]^2;

  /* second moment */

  M[.,2] = u2-h2;

  /* u(t) uncorrelated with x(t) */

  M[.,3] = x.*u;

  /* u(t)^2-h(t)^2 uncorrelated with u(t-i)^2, i = 1 */

  M[.,4] = (u2-h2).*lagn(u2,1);

  retp(M);
endp;

output file = gmm2b.out reset;

sv = 2|3|0.5;

/* theta[4] = 0, no arch effect */

RR = design(1|2|3|0); r = zeros(4,1);
_tsm_parnm = "beta1"|"beta2"|"alpha0"|"alpha1";

{theta1,stderr,Mcov,Qmin} = TD_cml(&ml,sv,RR,r);

_gmm_RR = RR;
_gmm_r = r;

{theta2,stderr,Mcov,Qmin} = gmm(&arch,sv);

output off;


