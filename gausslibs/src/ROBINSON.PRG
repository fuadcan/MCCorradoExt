/*
**  ROBINSON [1995], Gaussian semiparametric estimation of long range
**  dependence, Annals of Statistics, 23, 1630-1661
*/

new;
library tsm,optmum;
TSMset;

cls;

#include robinson.src;

output file = robinson.out reset;

_fourier = 0;
__output = 0;

d = -0.5;
do until d > 0.5;

  {x,retcode} = RND_arfima(d,0,0,2,1000,1000,1);

  H = Robinson(x,0,1,0);
  print ftos(d,"Fractal process simulation with d = %lf",1,2);
  print;
  print ftos(h,"Hurst exponent: %lf",4,3);
  d_ = H - 0.5;
  print; print ftos(d_,"Fractional difference: %lf",4,3);
  print; print; print;

  d = d + 0.05;

endo;

output off;
