/*
**  The problem of several maxima when estimating a fractional process
**
**  How to overcome this difficulty ?
**
*/

new;
library tsm,optmum,pgraph;

rndseed 1234;

Nobs = 250;

u = rndn(Nobs,1)*2;
u = u~(0|trimr(u,0,1));
u = u*(1|-0.5);
y = recserar(u,0,0.95);   /* Simulating an ARIMA(1,0,1) process */

sv = 0.95|0.5|0|2;
cnt = 0|0|1|0;            /* d is fixed */

_tsm_Mcov = 1;
_fourier = 0;
_print = 0;
__output = 1;

output file = arfima3.out reset;

_tsm_optmum = 0;   /* using the scoring algorithm */
{theta1,stderr,Mcov,Logl1} = arfima(y,1,1,sv,cnt);

_tsm_optmum = 1;   /* using OPTMUM library        */
{theta2,stderr,Mcov,Logl2} = arfima(y,1,1,sv,cnt);

print;
print "There is only one maximum when the fractional root d is fixed:";
print;
print "Parameter         Scoring         OPTMUM";
call printfmt(_tsm_parnm~theta1~theta2,0~1~1);
print;
print "So, d can be estimated before computing the MLE of the ARFIMA process.";
print "For example, we can use the GPH (Geweke Porter-Hudak) "\
      "estimate to fix d.";

output off;

