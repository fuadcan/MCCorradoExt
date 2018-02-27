/*
**  Consider the univariate MA(1) process:
**
**   y(t) = e(t) + 0.3*e(t-1)
**   with Var[e(t)] = 2
**
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, Cambridge
**
**  HARVEY, section 8.1.2
*/

new;
library tsm,optmum;
TSMset;

beta = -0.3;  /* phi1 = -0.3 */
sigma = 2;

call arma_autocov(beta,0,1,sigma,3);
acv0 =  varget("acv0");

/*
** acv0 corresponds to E[y(t)y(t)]
*/

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,0,1,sigma);
call SSM_build(Z,d,H,T,c,R,Q,0);

{a0,P0} = SSM_ic;

/*
** Exact Maximum Likelihood needs the computation of
**    E[y(t)y(t)]      (acv0)
**    E[y(t)e(t)] = 2  (acv1)
**    E[e(t)y(t)] = 2  (acv2)
**    E[e(t)e(t)] = 2  (acv3)
**
** P0 corresponds in fact to the following (2*2) matrix
**
**      acv0  acv1
**      acv2  acv3
*/

print P0;
print (acv0~2)|(2~2);

