/*
**  REINSEN [1994], Estimation of the fractional difference parameter
**  in the arima(p,d,q) model using the smoothed periodogram,
**  Journal of Time Series Analysis, 15, pages 335-350
**
**  See PRIESTLEY [1994], page 463 for the variances of the smoothed periodogram
*/

new;
library tsm,optmum;

T = 10000;

y = rndn(T,1);
_fourier = 0; {lambda,I} = PDGM(y);
M = 10;

d = zeros(7,1); std = zeros(7,1);
names = "GPH"|"Rectangular"|"Bartlett"|"Daniell"|"Tukey"|"Parzen"|"B-Priestley";
_smoothing = 2|0|1|0.23;

/* Variances of the smoothed periodogram */
var = zeros(7,1);
var[1] = pi^2/6;                /* periodogram              */
var[2] = 2;                     /* Truncated or Rectangular */
var[3] = 2/3;                   /* Bartlett                 */
var[4] = 1;                     /* Daniell                  */
var[5] = 2*(1-4*0.23+6*0.23^2); /* Tukey with a = 0.23; */
var[6] = 0.539285;              /* Parzen                   */
var[7] = 6/5;                   /* Bartlett-Priestley       */


{d_,std_} = GPH(I,lambda,M);
d[1] = d_; std[1] = std_;

j = 1;
do until j>6;
  _smoothing[2] = j;
  Is = smoothing(I);
  {d_,std_} = GPH(Is,lambda,M);
  d[j+1] = d_; std[j+1] = std_;
  j = j+1;
endo;

std = std.*sqrt(var);

/* Testing d = 0 */
tstat = d./std;
pvalue = 2*cdfnc(abs(tstat));

omat = names~d~std~tstat~pvalue;

output file = gph3.out reset;

print "                     d           std. err."\
      "        test d=0          p-value ";
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

