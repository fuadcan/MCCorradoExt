/*
** Exact Maximum Likelihood of the Linear Model with AR(1) errors
**
** BEACH, C.M. and J.G. MacKinnon [1978], A maximum likelihood
** procedure for regression with autocorrelated errors,
** Econometrica, 46, pages 51-58
**
*/

new;
library tsm,optmum;
declare external _ar1_data;

cls;

rndseed 1234;

Nobs = 5000;
u = recserar(rndn(Nobs,1),0,0.3);
x = 10+rndn(Nobs,1);
y = 2+3*x+u;
y[19] = miss(0,0);

_print = 1;
__output = 0;
_tsm_optmum = 1;
_tsm_Mcov = 3;

{theta1,stderr1,Mcov1,Logl1} = AR1_ml(y,ones(Nobs,1)~x,0,0);


/* Test rho = 0 */

RR = design(1|2|0|3); r = zeros(4,1);
{theta2,stderr2,Mcov2,Logl2} = AR1_ml(y,ones(Nobs,1)~x,RR,r);


H = vread(_ml_derivatives,"H_matrix");
G = vread(_ml_derivatives,"G_matrix");

I = -H;       /* Approximation of the Information matrix */

output file = tdml4a.out reset;

/*
**  Likelihood ratio statistic
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.06)
*/


LR = 2*(logl1-logl2);      /* Likelihood ratio                             */
pvalue = cdfchic(LR,1);    /* Approximating the noncentral chi-squared CDF */

print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/

LM = G'*inv(I)*G;
pvalue = cdfchic(LM,1);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


output off;


/*
** Maximum Likelihood of the Linear Model with AR(1) errors
** under the linear restriction theta = RR*gamma + r
**
** Set RR equal to 0 for unconstrained ML
**
*/


proc (4) = AR1_ml(y,x,RR,r);
  local k,data,retcode,y_,x_,sv,theta,stderr,Mcov,Logl;

  k = cols(x);
  if RR == 0;
    RR = eye(k+2);
    r = zeros(k+2,1);
  endif;

  _ar1_data = y~x;

  {data,retcode} = Missing(y~x,0);
  y_ = data[.,1];
  x_ = data[.,2:k+1];

  sv = y_/x_;
  sv = sv|0|stdc(y_-x_*sv);
  sv = inv(RR'RR)*RR'*(sv-r);

  __title = "Linear Model with AR(1) errors";

  _tsm_parnm = 0$+"X"$+ftocv(seqa(1,1,k),2,0);
  _tsm_parnm = _tsm_parnm|"rho"|"sigma";

  {theta,stderr,Mcov,Logl} = TD_cml(&_AR1_ml,sv,RR,r);

  retp(theta,stderr,Mcov,Logl);
endp;


proc _AR1_ml(theta);
  local data,k,y,x,beta,rho,sigma,T,u,du,logl;

  data = _ar1_data;
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





