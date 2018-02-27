/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 90-93
**
**  Kalman filtering and the Local Linear Trend Model
*/

new;
library tsm,optmum,pgraph;
load gnp[] = gnp.asc;

Z = 1~0; d = 0;
T = {1 1,0 1}; c = 0|0; R = {1 0,0 1};
a0 = gnp[1]|0; P0 = zeros(2,2);

proc ml(theta);
  local H,Q,LogL;
  theta = theta^2;  /* sigma^2 */
  H = theta[3];
  Q = (theta[1]~0)|(0~theta[2]);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(gnp,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;


_tsm_optmum = 0;    /* BHHH algorithm  */
_tsm_Mcov = 3;      /* Heteroskedasticity covariance */

output file = llt3.out reset;

{beta,stderr,Mcov,Logl} = TD_ml(&ml,5|5|5);

output off;






