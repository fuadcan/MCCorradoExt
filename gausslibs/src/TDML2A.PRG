/*
**  OZAKI [1982], The statistical analysis of perturbed limit cycle process
**  using nonlinear time series models, Journal of Times Series Analysis,
**  3, pages 29-41
**
**  See also:
**  TONG [1990], Non-linear Time Series, Oxford University Press
**
**  Estimate an EXPAR model
**  Note: Tong suspects that there has been a misprint or a computing error.
**        That is why we don't find OZAKI's estimates
**
*/

new;
library tsm,optmum;

load y[] = lynx.asc;
y = log(y);
x = y - meanc(y);

Nobs = rows(y);


proc EXPAR(theta);
  local beta,delta,sigma,resid,i,logL;

  beta = theta[1:6]; delta = theta[7]; sigma = theta[8];
  resid = miss(zeros(Nobs,1),0);
  LogL = miss(zeros(Nobs,1),0);

  i = 3;
  do until i>Nobs;

    resid[i] = x[i]
       -(beta[1] +(beta[2]+beta[3]*x[i-1])*exp(-delta*(x[i-1]^2)))*x[i-1]
       +(beta[4] +(beta[5]+beta[6]*x[i-1])*exp(-delta*(x[i-2]^2)))*x[i-2];
    i=i+1;
  endo;

  LogL[3:Nobs] = -0.5*ln(2*pi) -0.5*ln(sigma^2)
                 -0.5*(resid[3:Nobs]^2)/(sigma^2);

  retp(Logl);
endp;

output file = tdml2a.out reset;

_print = 1;
_tsm_optmum = 1;
_tsm_parnm  = "beta1"|"beta2"|"beta3"|"beta4"|"beta5"|
             "beta6"|"delta"|"sigma";
__title = "EXPAR model";

sv = ones(8,1);

{theta,stderr,Mcov,Logl} = TD_ml(&EXPAR,sv);

output off;




