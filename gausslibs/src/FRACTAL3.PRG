/*
**  Wornell and Oppenheim [1992], Estimation of fractal signals
**  from noisy measurements using wavelets,
**  IEEE transactions on Signals Processing, 40, 611-623
**
**
**  GEWEKE and PORTER-HUDAK [1983], The estimation and application
**  of long memory time series model,
**  Journal of Time Series Analysis, 4, pages 221-238
*/

new;
library tsm,optmum,pgraph;
#include fractal.src;
TSMset;

M = 8;
Nobs = 2^M;
Nr = 500;

data = rndn(Nobs,Nr)*4;

Estimate = zeros(Nr,2);

Nbsimul = 1;
do until Nbsimul > Nr;

  locate 5,5; print ftos(Nbsimul,"Iteration nø %lf",2,0);
  x = data[.,Nbsimul];

  /* GPH Method */

  d1 = GPH(x);

  /* Wavelet Method */

  d2 = WaveletMethod(x);

  Estimate[Nbsimul,.] = d1~d2;

  Nbsimul = Nbsimul + 1;

endo;


x = {}; dens ={};

j = 1;
do until j > 2;

  {x_,dens_,F,retcode} = Kernel(estimate[.,j]);

  x = x~x_; dens = dens~dens_;

  j = j + 1;

endo;


graphset;
   title("Monte Carlo Experiments"\
         "\LEmpirical density of the fractional "\
         "parameter for white noise process");
  _pdate = ""; _pnum = 2;
  _plegstr = "GPH method\0Wavelet Method";
  _plegctl = {2 6 1 4};
  _pltype = 6|1;
  xy(x,dens);




/*
** Procedure for the GPH method
*/

proc (1) = gph(x);
  local lambda,I,Nobs,M,Y,w,XX,beta,d;

  _fourier = 0;
  {lambda,I} = PDGM(x);

  Nobs = rows(I);
  M = trunc(sqrt(Nobs));

  Y = ln(I[2:M+1]);
  w = 4*sin(lambda[2:M+1]/2)^2;
  XX = ones(M,1)~-ln(w);

  beta = Y/XX;
  d = beta[2];

  retp(d);
endp;


/*
** Procedure for the Wavelet Method
*/

proc WaveletMethod(x);
  local H,G,Htilde,Gtilde,w;
  local sv,theta,stderr,Mcov,Logl;

  M = ln(rows(x))/ln(2);

  {H,G,Htilde,Gtilde} = Daubechies(4);

  w = wt(x,H,G,0);
  _wcoeff = w;
  _scale = seqa(1,1,M);
  sv = 4|0;

  __output = 0;
  _print = 0;
  _tsm_Mcov = 0;

  {theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

  retp(theta[2]);
endp;


