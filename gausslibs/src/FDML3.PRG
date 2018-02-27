new;
library tsm,optmum;
TSMset;

rndseed 123;

sv = 0.5|1|0.5;
sv = 1.5|0.5|2;

Nobs = 500;

_tsm_parnm = "theta1"|"sigma_u"|"sigma_v";
_tsm_Mcov = 1;
_fourier = 0;

z = Process;

output file = fdml3.out reset;

{coeff,stderr,Mcov,Logl} = FD_ml(z-lag1(z),&sgf,sv);

output off;

proc sgf(coeff,lambda);
  local phi1,theta1,sigma_u,sigma_v;
  local w,g;
  theta1 = coeff[1];
  sigma_u = coeff[2];
  sigma_v = coeff[3];
  w = cos(lambda);
  g = (sigma_u^2).*(1-2*theta1*w+theta1^2) +
      (sigma_v^2).*(2*(1-w));
  retp(g);
endp;

proc process;
  local u,x,v,y,z;
  u = rndn(Nobs,1)*1;
  x = recserar(u - 0.5*(0|trimr(u,0,1)),0,1);
  v = rndn(Nobs,1)*0.5;
  y = v;
  z = x + y;
  retp(z);
endp;

