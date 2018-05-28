/*
**  Wornell and Oppenheim [1992], Estimation of fractal signals
**  from noisy measurements using wavelets,
**  IEEE transactions on Signals Processing, 40, 611-623
*/

new;
library tsm,optmum,pgraph;
#include fractal.src;
TSMset;

{H,G,Htilde,Gtilde} = Daubechies(4);

M = 12;
Nobs = 2^M;
{x,retcode} = RND_arfima(0.25,0,0,2,1000,Nobs,1);

w = wt(x,H,G,0);
_wcoeff = w;

sv = 2|0.25;

_print = 1;

output file = fractal2.out reset;

s = 2;
do until s > M;
  print;
  print ftos(s,"dilation indice: %lf",1,0);
  print "---------------------------------------";
  print;
  _scale = seqa(1,1,s);
  {theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);
  s = s + 1;
endo;

output off;


