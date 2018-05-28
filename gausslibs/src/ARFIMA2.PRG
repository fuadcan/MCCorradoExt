/*
**  The problem of several maxima when estimating a fractional process
*/

new;
library tsm,optmum,pgraph;

rndseed 1234;

Nobs = 250;               /* Number of observations */

u = rndn(Nobs,1)*2;
u = u~(0|trimr(u,0,1));
u = u*(1|-0.5);
y = recserar(u,0,0.95);   /* Simulating an ARIMA(1,0,1) process */

sv = 0.95|0.5|0|2;
cnt = 0|0|0|0;            /* d is to be estimated */

_tsm_Mcov = 1;
_fourier = 0;
_print = 0;
__output = 1;

output file = arfima2.out reset;

_tsm_optmum = 0;   /* using the scoring algorithm */
{theta1,stderr,Mcov,Logl1} = arfima(y,1,1,sv,cnt);

_tsm_optmum = 1;   /* using OPTMUM library        */
{theta2,stderr,Mcov,Logl2} = arfima(y,1,1,sv,cnt);

tv = 0.95|0.5|0|2; /* True value */

print;
print "There is two maxima:";
print;
print "Parameter         True value      Scoring         OPTMUM";
call printfmt(_tsm_parnm~tv~theta1~theta2,0~1~1~1);

output off;

lambda = seqa(0,pi/100,100);
g1 = _arfima_sgf(theta1,lambda);
g2 = _arfima_sgf(theta2,lambda);

graphset;
  fonts("simplex simgrma");
  _pdate = ""; _pnum = 2;
  title("But the two spectral generating functions are identical !");
  ylabel("spectrum");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  logy(lambda,g1~g2);

