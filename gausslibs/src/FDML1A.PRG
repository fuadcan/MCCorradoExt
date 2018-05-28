new;
library tsm,optmum,pgraph;
TSMset;

declare external r;

rndseed 123;
y = recserar(rndn(500,1),0|0,0.5|0.4);

_fourier = 1;
_tsm_optmum = 0;
_tsm_Mcov = 1;
r = 4;
sv = ones(r+1,1);

_tsm_parnm = "gamma1"|"gamma2"|"gamma3"|"gamma4"|"sigma";

output file = fdml1a.out reset;

{theta,stderr,Mcov,Logl} = FD_ml(y,&sgf,sv);

output off;

{lambda,I} = PDGM(y);
_smoothing = 5|5|0|0.23;   /* Parzen lag window with bandwidth = 2 */
I = smoothing(I);
g = sgf(theta,lambda);
q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  xlabel("frequency");
  _plegstr = "Periodogram of the data"\
           "\000Estimated spectral generating function";
  _plegctl = {2 5 3 4};
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda[1:q],I[1:q]~g[1:q]);


proc sgf(theta,lambda);       /* Dzhaparidze, page 125 */
  local N,gamma_,sigma,j,w,g;
  N = rows(lambda);
  gamma_ = theta[1:r];
  sigma = theta[r+1];
  j = seqa(1,1,r);
  w = lambda.*j';
  w = cos(w);
  g = (sigma^2)*exp(2*w*gamma_);
  retp(g);
endp;
