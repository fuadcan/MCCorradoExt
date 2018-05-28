/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, Cambridge
**
**  HARVEY, section 8.1.2
*/

new;
library tsm,optmum;
TSMset;

load beta = ARMA1;
load SIGMA = SIGMA1;

call arma_autocov(beta,1,1,SIGMA,3);
acv0 =  varget("acv0");

/*
** acv0 corresponds to E[y(t)y(t)]
*/

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,sigma);
call SSM_build(Z,d,H,T,c,R,Q,0);

{a0,P0} = SSM_ic;

/*
** Exact Maximum Likelihood needs the computation of
**    E[y(t)y(t)]      (acv0)
**    E[y(t)e(t)] = SIGMA  (acv1)
**    E[e(t)y(t)] = SIGMA  (acv2)
**    E[e(t)e(t)] = SIGMA  (acv3)
*/

print P0;
print (acv0~SIGMA)|(SIGMA~SIGMA);


