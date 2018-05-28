/*
**  Monte Carlo investigations of the GPH estimator's power
*/

new;
library tsm,optmum,pgraph;

T = 1000;
Ns = 250; /* Number of simulations */

_fourier = 0; M = trunc(sqrt(T)); d = zeros(Ns,3);
_smoothing = 2|0|1|0.20; /* a = .20 */

j = 1;
do until j>Ns;

  y = rndn(T,1);
  {lambda,I} = PDGM(y);

  d1 = GPH(I,lambda,M);      /* Geweke-Porter-Hudak     */

  _smoothing[2] = 2;         /* Bartlett                */
  Is1 = smoothing(I);
  d2 = GPH(Is1,lambda,M);

  _smoothing[2] = 4;         /* Tukey with a = .20      */
  Is2 = smoothing(I);
  d3 = GPH(Is2,lambda,M);

  d[j,.] = d1~d2~d3;
  j = j+1;
endo;


let _Kernel = {-0.75,0.75,128,2,.,.}; /* Epanechnikov Kernel */

dens = {};

j = 1;
do until j>3;
  {x,dens_,F,retcode} = Kernel(d[.,j]);
  dens = dens~dens_;
  j = j + 1;
endo;

graphset;
title("Monte Carlo investigations"\
      "\LDensity of the fractional root d estimators");
_pdate = ""; _pnum = 2; _pltype = 6|1|3;
_plegstr = "GPH\0Bartlett\0Tukey"; _plegctl = 1;
xy(x,dens);


/*
** Procedure for the GPH estimator
*/

proc (1) = gph(I,lambda,M);
  local Y,w,X,beta,d;

  Y = ln(I[2:M+1]);
  w = 4*sin(lambda[2:M+1]/2)^2;
  X = ones(M,1)~-ln(w);

  beta = Y/X;
  d = beta[2];

  retp(d);
endp;


