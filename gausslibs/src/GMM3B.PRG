new;
library tsm,optmum;
TSMset;

X0 = 10;
a = 0.8;
b = 0.1;
sigma = 0.06;
h = 0.1;
Nobs = 1000;
sv = a|b|sigma;

/* Generate an Ornstein-Uhlenbeck process */

xt = zeros(Nobs,1);
xt[1] = x0;

k1 = exp(-a*h);
k2 = (1-exp(-2*a*h))/(2*a);

u = (sigma*sqrt(k2))*rndn(Nobs-1,1);
u = u+b*(1-k1);

i = 2;
do until i > Nobs;
  xt[i] = k1*xt[i-1]+u[i-1];
  i = i + 1;
endo;


proc ml(theta);
  local a,b,sigma,k1,k2,epsilon,logl;

  a = theta[1];
  b = theta[2];
  sigma = theta[3];

  k1 = exp(-a*h);
  k2 = sigma^2*(1-exp(-2*a*h))/(2*a);

  epsilon = xt[2:Nobs]-k1.*xt[1:Nobs-1]-b*(1-k1);

  logl = -0.5*ln(2*pi) - 0.5*ln(k2) - 0.5*epsilon^2/k2;

  retp(logl);
endp;

proc cm(theta);
  local a,b,sigma,k1,k2,epsilon,M;

  a = theta[1];
  b = theta[2];
  sigma = theta[3];

  k1 = exp(-a*h);
  k2 = sigma^2*(1-exp(-2*a*h))/(2*a);

  epsilon = xt[2:Nobs]-k1.*xt[1:Nobs-1]-b*(1-k1);

  M = epsilon~(epsilon^2-k2)~(epsilon.*xt[1:Nobs-1]);

  retp(M);
endp;


output file = gmm3b.out reset;

_tsm_parnm = "a"|"b"|"sigma";

_ml_Jacobian_proc = 0;
__title = "Ornstein-Ulhenbeck process";

{theta1,stderr,Mcov,Logl} = TD_ml(&ml,sv);

{theta2,stderr,Mcov,Qmin} = gmm(&cm,sv);

print; print;
print "                 True values        TDML            GMM";
call printfmt(_tsm_parnm~sv~theta1~theta2,0~1~1~1);

output off;



