/*
**  Moulin [1994], Wavelet thresholding techniques for
**  power spectrum estimation, IEEE Transactions on Signal
**  Processing, 42, 3126-3136
**
**  Gao, H-Y. [1996], Choice of thresholds for wavelet
**  shrinkage estimate of the spectrum, Research report,
**  StatSci Division of MathSoft
**
*/

new;
library tsm,pgraph;
TSMset;

rndseed 123456;

t = seqa(1,1,1024);
beta = -0.9|-0.5;
sigma = 0.5;
x = recserar(rndn(1024,1)*sigma,0|0,beta);
x = x + rndn(1024,1)*0.25;


_fourier = 0;
{lambda,I} = PDGM(x);

_smoothing = 25|5|0|0.23;   /* Parzen lag window with bandwidth = 25 */
I2 = smoothing(I);
sgf = arma_sgf(beta|sigma,lambda,2,0);

{H,G,Htilde,Gtilde} = Daubechies(8);
_wcenter = 0;

w = wt(I,H,G,0);
w = extract(w,11|10|9|8|7|6);
y = iwt(w,Htilde,Gtilde,0);

q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _ptitlht = 0.25;
  _pltype = 6|1|3|5;

  begwind;
  window(2,1,0);
  setwind(1);
    title("Periodogram");
    xy(lambda[1:q],I[1:q]);
  nextwind;
    title("Power Spectrum Estimation");
    ytics(0,2,1,0);
    _plegstr = "spectral generating function\0Parzen window\0Wavelet";
    _plegctl = {2 5 1 3};
    xy(lambda[1:q],sgf[1:q]~I2[1:q]~y[1:q]);
  endwind;



/*
** spectral generating function for ARMA process
*/

proc arma_sgf(beta,lambda,p,q);
  local Np,M,phi,theta,d,sigma;
  local w1,w2,w3,w4,w5,w6,g;

  Np = p+q+1;

  if p /= 0;
    phi = beta[1:p];
  endif;

  if q /= 0;
    theta = beta[p+1:p+q];
  endif;

  sigma = beta[p+q+1];

  M = maxc(p|q);

  w1 = lambda.*seqa(1,1,M)';
  w2 = cos(w1);
  w3 = sin(w1);

  if q /= 0;
    w4 = (1 - w2[.,1:q]*theta)^2 + (w3[.,1:q]*theta)^2;
  else;
    w4 = 1;
  endif;

  if p /= 0;
    w5 = (1 - w2[.,1:p]*phi)^2 + (w3[.,1:p]*phi)^2;
  else;
   w5 = 1;
  endif;

  w6 = w4./w5;
  g = (sigma^2)*w6;

  retp(g);
endp;

