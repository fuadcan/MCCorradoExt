/*
**  REINSEN [1994], Estimation of the fractional difference parameter
**  in the arima(p,d,q) model using the smoothed periodogram,
**  Journal of Time Series Analysis, 15, pages 335-350
**
**  See PRIESTLEY [1994], page 463 for the variances of the smoothed periodogram
*/

new;
library tsm,optmum;

T = 1000;

y = recserar(rndn(T,1),10,1); /* Simulation of a random walk process  */
_fourier = 1;                 /* Using the discrete Fourier transform */
{lambda,I} = PDGM(y);

_smoothing = 4|5|1|0.23;      /* Bandwidth = 4 and Parzen lag window  */
Is = smoothing(I);            /* Smoothed periodogram                 */

M = trunc(sqrt(T));

{d1,std} = GPH(I,lambda,M);   /* GPH estimator                            */
std1 = std*sqrt(pi^2/6);
{d2,std} = GPH(Is,lambda,M);  /* REINSEN estimator with Parzen lag window */
std2 = std*sqrt(0.53928);

d = d1|d2;
std = std1|std2;

/* Testing d = 1 */
tstat = (d-1)./std;
pvalue = 2*cdfnc(abs(tstat));

names = "GPH"|"Parzen";
omat = names~d~std~tstat~pvalue;

output file = gph2.out reset;

print "                     d           std. err."\
      "        test d=1          p-value ";
print;
call printfmt(omat,0~1~1~1~1);

output off;

/*
** Procedure for the GPH estimator
*/

proc (2) = gph(I,lambda,M);
  local Y,w,X,beta,d,XXinv,std_;

  Y = ln(I[2:M+1]);
  w = 4*sin(lambda[2:M+1]/2)^2;
  X = ones(M,1)~-ln(w);

  beta = Y/X;
  d = beta[2];

  XXinv = invpd(X'X);
  std_ = sqrt(XXinv[2,2]);
  retp(d,std_);
endp;
