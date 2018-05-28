/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 86-89
**
**  Cycle Model
*/

new;
library tsm,optmum,pgraph;

load rainfall[] = rainfall.asc;

y = rainfall[1:131];
y = y - 142.2;

_cycle_prmt = 0;
_tsm_parnm = "rho"|"lambda_c"|"sig_kappa"|"sig_epsilon";
_tsm_optmum = 1;
_print = 1;

/*
**  Using Harvey's starting values, page 198
*/

rho = 0.5; lambda_c = 0.78; sig_epsilon = 1; sig_kappa = 1;
sv = rho|lambda_c|sig_epsilon|sig_kappa;

{theta,stderr,Mcov,Logl} =  FD_ml(y,&sgf,sv);




{lambda,I} = PDGM(y);
_smoothing = 6|5|0|0.23;   /* Parzen lag window with bandwidth = 6 */
I = smoothing(I);
g = sgf(theta,lambda);
q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  _plegstr = "Periodogram of the Rainfall data"\
             "\000Estimated spectral generating function";
  _plegctl = {2 5 3 4};
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda[1:q],I[1:q]~g[1:q]);

/*
**  PAWITAN and O'SULLIVAN [1994], Nonparametric spectral density estimation
**  using penalized whittle likelihood,
**  Journal of the American Statistical Association, pages 600-610
*/

/*
**  I: periodogram
**  f: spectral density function
**
**  The autors use the fact that 2I/f is asymptotically
**  distributed as a chi-squared function with 2 degrees of freedom
**  to plot the observed quantiles of 2I/f against
**  the quantiles of the chi(2) distribution.
*/

z = 2*I./g;

z = sortc(z,1); Nobs = rows(z);
CDF = seqa(1,1,Nobs)/Nobs;
x = seqa(0,10/Nobs,Nobs);
CDF_ = 1-cdfchic(x,2);

graphset;
  _pdate = ""; _pnum = 2;
  fonts("simplex simgrma");
  _plegstr = "\202h\201[2](2)\000empirical distribution of 2I/f";
  _plegctl = {2 6 2.5 2};
  let w = {.,.};
  scale(0|10,w);
  xy(x~z,CDF_~CDF);


/*
**  Spectral generating function for a stochastic cycle plus noise model
*/

proc sgf(theta,lambda);
  local g_cycle,g_noise,g;
  g_cycle = _cycle_sgf(theta[1:3],lambda);
  g_noise = theta[4]^2;
  g = g_cycle + g_noise;
  retp(g);
endp;





