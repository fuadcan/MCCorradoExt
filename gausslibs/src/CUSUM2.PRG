/*
**  CUSUM and CUSUMsq tests
**
**  SPANOS [1986], Statistical Foundations of Econometric Modelling,
**  Cambridge University Press, page 477
**
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 256-257
*/

new;
library optmum,tsm,pgraph;

rndseed 123456;
Nobs = 100;
u = rndn(100,1)*2;
u = u~(0|trimr(u,0,1))~(0|0|trimr(u,0,2));
data = u*(1|-0.5|0.25);                      /* MA(2) process */

_print = 1;
_tsm_optmum = 1;
_tsm_gtol = 0.001;

{beta_,stderr_,Mcov_,LogL_} = arma_ML(data,0,2,"HR2");
SIGMA_ = _arma_SIGMA;

/*  Convert the ARMA process to a state space model */

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta_,0,2,SIGMA_);
call SSM_build(Z,d,H,T,c,R,Q,0);

/*  Starting values for the exact Maximum Likelihood estimation */

coeff = beta_|sqrt(SIGMA_);


/*  Build a procedure for the exact likelihood function */

proc ml(coeff);
  local beta,sigma,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  beta = coeff[1:2];
  SIGMA = coeff[3]^2;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,0,2,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(data,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

{beta,stderr,Mcov,Logl} = TD_ml(&ml,Coeff);

SIGMA = beta[3]^2; beta = beta[1:2];
{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,0,2,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);
{a0,P0} = SSM_ic;
call KFiltering(data,a0,P0);

v = KF_matrix(2);
F = KF_matrix(7);


w = v./sqrt(F);   /* Standardized innovations */
t_ = seqa(1,1,Nobs);
sw = stdc(w);


/* CUSUM statistic */
cusum = cumsumc(w)/stdc(w);         /* Wt statistic, Spanos, formula 21.168  */
                                    /* Harvey, formula 5.4.3a                */

Bound = (-1.96~1.96).*sqrt(t_);     /* distribution approximated by N(0,t_)  */

graphset;
title("CUSUM test --- significance level: 5%");
_pdate = "";
xy(t_,cusum~Bound);

a = 0.948; /* 5% */

Bound = (-a~a).*(sqrt(Nobs)+2*t_/sqrt(Nobs));   /* Spanos, formula 21.169   */
                                                /* Harvey, page 257         */

graphset;
title("CUSUM test --- significance level: 5%");
_pdate = "";
xy(t_,cusum~Bound);


/* CUSUMsq statistic */
cusumsq = cumsumc(w^2)/sumc(w^2);   /* Vt statistic, Spanos, formula 21.170 */
center = t_/Nobs;

graphset;
title("CUSUMsq test");
_pdate = "";
xy(t_,cusumsq~center);


