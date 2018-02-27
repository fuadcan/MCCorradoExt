/*
**  Chan, Karolyi, Longstaff and Sanders [1992], An empirical
**  comparison of alternative models of the short-term interest rate,
**  Journal of Finance, 47, 1209-1227
**
*/

new;
library tsm,optmum;
TSMset;

rndseed 123;

y0 = 10;
a = 0.8;
b = 0.1;
sigma = 0.06;
h = 0.1;
Nobs = 1000;

/* Generate an Ornstein-Ulhenbeck process */

y = zeros(Nobs,1);
y[1] = y0;

k1 = exp(-a*h);
k2 = (1-exp(-2*a*h))/(2*a);

u = (sigma*sqrt(k2))*rndn(Nobs-1,1);
u = u+b*(1-k1);

i = 2;
do until i > Nobs;
  y[i] = k1*y[i-1]+u[i-1];
  i = i + 1;
endo;


proc CKLS(theta);
  local alpha,beta,gamma_,sigma;
  local et,vt,M1,M2,M3,M4,M;

  alpha = theta[1];
  beta = theta[2];
  gamma_ = theta[3];
  sigma = theta[4];

  /* equations 2 and 3, page 1213 */

  et = y[2:Nobs] - y[1:Nobs-1] - (alpha+beta*y[1:Nobs-1])*h;
  vt = h*(sigma*abs(y[1:Nobs-1])^gamma_)^2;

  /* equation 4, page 1214 */

  M1 = et;
  M2 = et.*y[1:Nobs-1];
  M3 = et^2 - vt;
  M4 = (et^2 - vt).*y[1:Nobs-1];

  M = M1~M2~M3~M4;

  retp(M);
endp;


output file = gmm3c.out reset;

_tsm_parnm = "alpha"|"beta"|"gamma"|"sigma";

sv1 = a*b|-a|0|sigma;
{theta1,stderr,Mcov,Qmin} = gmm(&CKLS,sv1);

sv2 = a*b|-a|sigma;
_gmm_RR = design(1|2|0|3);
_gmm_r = zeros(4,1);
{theta2,stderr,Mcov,Qmin} = gmm(&CKLS,sv2);

print; print;
print "                 True values       GMM          const. GMM";
call printfmt(_tsm_parnm~sv1~theta1~theta2,0~1~1~1);

output off;




