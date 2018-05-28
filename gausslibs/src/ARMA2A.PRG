/*
**  Consider the univariate AR(1) process:
**
**   y(t) = 0.5*y(t-1) + e(t)
**   with Var[e(t)] = 1
*/

new;
library tsm,optmum;
TSMset;

beta = 0.5;
sigma = 1;

call arma_autocov(beta,1,0,sigma,3);
acv0 =  varget("acv0");

/*
** acv0 corresponds to E[y(t)y(t)]
*/

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,0,sigma);
call SSM_build(Z,d,H,T,c,R,Q,0);

{a0,P0} = SSM_ic;

/*
** Exact Maximum Likelihood needs the computation of
**    E[y(t)y(t)]      (acv0)
**    E[y(t)e(t)] = 1  (acv1)
**    E[e(t)y(t)] = 1  (acv2)
**    E[e(t)e(t)] = 1  (acv3)
**
** P0 corresponds in fact to the following (2*2) matrix
**
**      acv0  acv1
**      acv2  acv3
*/

print P0;
print (acv0~1)|(1~1);
