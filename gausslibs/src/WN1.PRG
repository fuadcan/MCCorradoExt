/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press
*/


new;
library tsm,optmum,pgraph;

rndseed 123456;

y = rndn(1000,1)*2;

proc sgf(beta,lambda);       /* Harvey, page 194, formula 4.3.15 */
  local sigma,N,g;
  sigma = beta; N = rows(lambda);
  g = ones(N,1)*(sigma^2);
  retp(g);
endp;

proc Jproc(beta,lambda);
  local sigma,N,J;
  sigma = beta; N = rows(lambda);
  J = 2*ones(N,1)*sigma;
  retp(J);
endp;


_fourier = 0;
_tsm_optmum = 0;
_sgf_Jacobian_proc = &Jproc;
{sigma,stderr,Mcov,Logl} = FD_ml(y,&sgf,1);

/*
**  PAWITAN and O'SULLIVAN [1994], Nonparametric spectral density estimation
**  using penalized whittle likelihood,
**  Journal of the American Statistical Association, 89, pages 600-610
*/

/*
**  I: periodogram
**  f: spectral density
**
**  The autors use the fact that 2I/f is asymptotically
**  distributed as a chi-squared function with 2 degrees of freedom
**  to plot the observed quantiles of 2I/f against
**  the quantiles of the chi(2) distribution
*/

{lambda,I} = PDGM(y);
f = sgf(sigma,lambda);

q = trunc(rows(lambda)/2);
z = 2*I[1:q]./f[1:q];

z = sortc(z,1); Nobs = rows(z);
CDF = seqa(1,1,Nobs)/Nobs;
x = seqa(0,maxc(z)/Nobs,Nobs);
CDF_ = 1-cdfchic(x,2);

graphset;
  title("Is that really a white noise process?");
  _pdate = ""; _pnum = 2;
  fonts("simplex simgrma");
  _plegstr = "\202h\201[2](2)\000empirical distribution of 2I/f";
  _plegctl = {2 6 2.5 2};
  let w = {.,.};
  scale(0|10,w);
  xy(x~z,CDF_~CDF);

