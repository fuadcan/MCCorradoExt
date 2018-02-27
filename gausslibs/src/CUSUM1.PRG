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
library tsm,optmum,pgraph;

load purse[] =  purse.asc;
Nobs = rows(purse);

{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

Z = 1; d = 0; H = beta[2]^2;
T = 1; c = 0; R = 1; Q = beta[1]^2;

call SSM_build(Z,d,H,T,c,R,Q,0);

a_0 = 0;
P_0 = 0;

call KFiltering(purse,a_0,P_0);

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


