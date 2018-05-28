/*
** Estimation of an univariate ARMA-GARCH process with
** the arma_Filter and garch_Filter procedures
**
** The ARMA(2,1)-GARCH(1,1) model is
**
**     y(t) = phi1*y(t-1) + phi2*y(t_2) + epsilon(t) - theta1*epsilon(t-1)
**
**     epsilon(t) --> N(0,h(t)^2)
**
**     h(t)^2 = alpha0^2 + alpha1^2*epsilon(t-1)^2 + alpha2^2*h(t-1)^2
**
*/

new;
library tsm,optmum,pgraph;
TSMset;

Nobs = 200;
rndseed 123456;
y = RND_arma(0.6|0.3|0.25,2,1,1,miss(0,0),Nobs);

proc ml(beta);
  local e,h2,logL;

  e = arma_filter(y,beta[1:3],2,1);
  e = missrv(e,0);
  h2 = garch_filter(e,beta[4:6]^2,1,1);
  logL = -0.5*ln(2*pi)-0.5*ln(h2)-0.5*(e^2)./h2;

  retp(logL);
endp;

output file = filter1b.out reset;

_tsm_parnm  = "Phi1"|"Phi2"|"Theta1"|"alpha0"|"alpha1"|"alpha2";
_tsm_Mcov = 2;
sv = 0.6|0.3|0.25|1|0.5|0.5;
{beta,stderr,Mcov,LogL} = TD_ml(&ml,sv);

output off;

e = arma_filter(y,beta[1:3],2,1);
e = missrv(e,0);
h2 = garch_filter(e,beta[4:6]^2,1,1);

graphset;
  _pdate = ""; _pnum = 2;
  title("GARCH component h(t)[2]");
  xlabel("observation");
  xy(seqa(1,1,Nobs),h2);
