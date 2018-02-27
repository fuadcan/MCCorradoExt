new;
library tsm,optmum;

rndseed 123;

Nobs = 1000;             /* Number of observations */

u = rndn(Nobs,1)*2;
u = u~(0|trimr(u,0,1));
u = u*(1|-0.5);
y = recserar(u,0,0.95);   /* Simulating an ARIMA(1,0,1) process */

sv = 0.95|0.5|0|2;
cnt = 0|0|0|0;

_tsm_Mcov = 1;
_fourier = 0;
_print = 1;

output file = arfima1.out reset;

_tsm_optmum = 0;   /* using the scoring algorithm */

{theta1,stderr,Mcov,Logl} = arfima(y,1,1,sv,cnt);

output off;


