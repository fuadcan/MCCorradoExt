/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 89-90
**
**  Estimation in the time domain
*/

new;
library tsm,optmum,pgraph;
declare external a_0,P_0;

load purse[] = purse.asc;
Nobs = rows(purse);

Z = 1; d = 0;
T = 1; c = 0; R = 1;

proc ml(theta);
  local H,Q,LogL;
  H = theta[2]^2;
  Q = theta[1]^2;
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(purse,a_0,P_0);
  LogL = KF_ml;
  retp(LogL);
endp;


_print = 1;
_tsm_optmum = 1;

output file = ll3.out reset;

/* First estimation */
a_0 = 0; P_0 = 0;
{beta,stderr,Mcov,Logl} = TD_ml(&ml,1|1);

/* Second estimation */
a_0 = purse[1]; P_0 = 0;
{beta,stderr,Mcov,Logl} = TD_ml(&ml,1|1);

output off;


