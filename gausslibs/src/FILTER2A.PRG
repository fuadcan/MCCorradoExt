/*
** Estimation of a fractional ARMA(2,1) process
*/

new;
library tsm,optmum;
TSMset;

Nobs = 200;
rndseed 123456;
y = RND_arma(0.6|0.3|0.25,2,1,1,miss(0,0),Nobs);

proc ml(beta);
  local e,sigma2,logL;

  e = fractional_Filter(y,beta[4]);      /* d = beta[4] */
  e = arma_filter(e,beta[1:3],2,1);
  e = missrv(e,0);
  sigma2 = beta[5];
  logL = -0.5*ln(2*pi)-0.5*ln(sigma2)-0.5*(e^2)./sigma2;

  retp(logL);
endp;


output file = filter2a.out reset;

_tsm_parnm  = "Phi1"|"Phi2"|"Theta1"|"d"|"sigma";
_tsm_Mcov = 1;
sv = rndu(5,1);
_print = 1;

Print "Time-Domain Estimation (CSS method):";
print "------------------------------------";

{beta,stderr,Mcov,LogL} = TD_ml(&ml,sv);

print; print;
Print "Frequency-Domain Estimation (Whittle method):";
print "---------------------------------------------";

_fourier = 0;
{beta,stderr,Mcov,Logl} = arfima(y,2,1,sv,zeros(5,1));

output off;
