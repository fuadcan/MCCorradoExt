new;
library tsm,optmum;

rndseed 123456;

beta = 0.8|0.4|0.25;
{y,retcode} = RND_arfima(beta,1,1,2,1000,100,1);

sv = 0.8|0.4|0.25|2;
cnt = 0|0|0|0;       /* No fixed parameters           */

_tsm_Mcov = 1;
_fourier = 0;
_print = 0;

{beta1,stderr,Mcov,Logl1} = arfima(y,1,1,sv,cnt);

sv = 0.8|0.4|0.25|2;
cnt = 0|0|1|0;       /* d is fixed                    */

{beta2,stderr,Mcov,Logl2} = arfima(y,1,1,sv,cnt);

sv = 0.8|0.4|0.25|2;
cnt = 0|1|1|0;       /* d and the MA parameter are fixed */

{beta3,stderr,Mcov,Logl3} = arfima(y,1,1,sv,cnt);

/*
**  Likelihood ratio statistic
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.06)
*/

output file = arfima4.out reset;

LR = 2*(logl1-logl2);      /* Likelihood ratio                             */
pvalue = cdfchic(LR,1);    /* Approximating the noncentral chi-squared CDF */

print "Testing d = 0.25"; print chrs(45*ones(40,1));
print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);
print; print;


LR = 2*(logl1-logl3);      /* Likelihood ratio                             */
pvalue = cdfchic(LR,2);    /* Approximating the noncentral chi-squared CDF */

print "Testing d = 0.25 and theta1 = 0.4"; print chrs(45*ones(40,1));
print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;





