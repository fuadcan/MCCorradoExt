/*
** Exact Maximum Likelihood of the Linear Model with AR(1) errors
**
** BEACH, C.M. and J.G. MacKinnon [1978], A maximum likelihood
** procedure for regression with autocorrelated errors,
** Econometrica, 46, pages 51-58
**
*/

new;
library tsm,optmum,CML;

rndseed 1234;

u = recserar(rndn(100,1),0,0.3);
x = 10+rndn(100,1);
y = 2+3*x+u;

output file = cml1.out reset;

{theta1,f,g,cov,ret} = CMLPrt(CML(y~ones(100,1)~x,0,&_AR1_ml,ones(4,1)));

/* Test rho = 0 */

_CML_A = {0 0 1 0}; _CML_B = {0};

{theta2,f,g,cov,ret} = CMLPrt(CML(y~ones(100,1)~x,0,&_AR1_ml,ones(4,1)));

output off;

/*
** Maximum Likelihood of the Linear Model with AR(1) errors
** under the linear restriction theta = RR*gamma + r
**
** Set RR equal to 0 for unconstrained ML
**
*/


proc _AR1_ml(theta,data);
  local k,y,x,beta,rho,sigma,T,u,du,logl;

  k = cols(data) - 1;
  y =  data[.,1]; x = data[.,2:k+1];

  beta = theta[1:k]; rho = theta[k+1]; sigma = theta[k+2];

  T = rows(y);
  u = zeros(T,1);

  u = y-x*beta;
  du = u[2:T]-rho*u[1:T-1];

  logl  = zeros(T,1);
  logl[2:T] = -0.5*ln(2*pi)-0.5*ln(sigma^2)-0.5*(du.*du)./(sigma^2);

  /* Log-likelihood for the first observation */

  logl[1] = -0.5*ln(2*pi)-0.5*ln(sigma^2)-0.5*(1-rho^2)*(u[1]^2)/(sigma^2);

  retp(logl);
endp;

