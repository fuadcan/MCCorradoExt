/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 90-93
**
**  Estimation in the frequency domain
*/

new;
library tsm,optmum;
load gnp[] = gnp.asc;


output file = LLT1.out reset;

_tsm_optmum = 0;   /* Method of scoring */

{beta,stderr,Mcov,Logl} =  sm_LLT(gnp,1|1|1);

_tsm_optmum = 1;   /* OPTMUM library    */
_tsm_Mcov = 0;
{beta,stderr,Mcov,Logl} =  sm_LLT(gnp,1|1|1);
_tsm_Mcov = 1;
{beta,stderr,Mcov,Logl} =  sm_LLT(gnp,abs(beta));

output off;




