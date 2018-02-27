new;
library tsm,optmum;
TSMset;

X0 = 10;
mu = 0.2;
sigma = 0.5;
h = 0.1;
Nobs = 500;
sv = mu|sigma;

/* Generate a geometric Brownian motion */

t = seqa(0,h,Nobs);
xt = x0*exp(   (mu-0.5*sigma^2)*t
               + sigma*cumsumc(0|rndn(Nobs-1,1)*sqrt(h))
          );


proc ml(theta);
  local e,mu,sigma,epsilon,logl;

  e = ln(xt[2:Nobs]./xt[1:Nobs-1]);

  mu = theta[1];
  sigma = theta[2];

  epsilon = e - (mu-0.5*sigma^2)*h;

  logl = -0.5*ln(2*pi) - 0.5*ln(sigma^2*h) - 0.5*epsilon^2/(sigma^2*h);

  retp(logl);
endp;

proc ml_jproc(theta);
  local e,mu,sigma,epsilon,J1,J2,J;

  e = ln(xt[2:Nobs]./xt[1:Nobs-1]);

  mu = theta[1];
  sigma = theta[2];

  epsilon = e - (mu-0.5*sigma^2)*h;

  J1 = epsilon/(sigma^2);
  J2 = -1/sigma + epsilon.*(epsilon-sigma^2*h)/(sigma^3*h);

  J = J1~J2;
  retp(J);
endp;


proc cm(theta);
  local e,mu,sigma,epsilon,M;

  e = ln(xt[2:Nobs]./xt[1:Nobs-1]);

  mu = theta[1];
  sigma = theta[2];

  epsilon = e - (mu-0.5*sigma^2)*h;

  M = epsilon~(epsilon^2-(sigma^2*h));

  retp(M);
endp;

proc gmm_jproc(theta);
  local e,mu,sigma,epsilon,J1,J2,J;

  e = ln(xt[2:Nobs]./xt[1:Nobs-1]);

  mu = theta[1];
  sigma = theta[2];

  epsilon = e - (mu-0.5*sigma^2)*h;

  J1 = -h~(sigma*h);
  J2 = 2*meanc(epsilon).*J1;
  J2[1,2] = J2[1,2] - 2*sigma*h;

  J = J1|J2;
  retp(J);
endp;


output file = gmm3a.out reset;

_tsm_parnm = "mu"|"sigma";

_ml_Jacobian_proc = 0;
__title = "Geometric Brownian Motion";

{theta1,stderr,Mcov,Logl} = TD_ml(&ml,sv);

_ml_Jacobian_proc = &ml_jproc;
{theta2,stderr,Mcov,Logl} = TD_ml(&ml,sv);

{theta3,stderr,Mcov,Qmin} = gmm(&cm,sv);

_gmm_Jacobian_proc = &gmm_jproc;
{theta4,stderr,Mcov,Qmin} = gmm(&cm,sv);

print; print;
print "                 True values        TDML            GMM";
call printfmt(_tsm_parnm~sv~theta1~theta3,0~1~1~1);

output off;

