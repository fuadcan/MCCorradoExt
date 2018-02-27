/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 89-90
**
**  Estimation in the frequency domain
*/


new;
library tsm,optmum;

load purse[] = purse.asc;

_tsm_optmum = 0;   /* Method of scoring */

output file = LL1.out reset;

{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

print; print;

_tsm_optmum = 1;   /* OPTMUM library    */

{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

output off;

