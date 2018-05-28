/*
**  TERASVIRTA [1994], Specification, Estimation, and Evaluation of
**  Smooth Transition Autoregressive Models, Journal
**  of the American Statistical Association, 89, pages 208-218
**
**  Estimate a LSTAR model
**  See TERASVIRTA, page 216
*/

new;
library tsm,optmum;

load y[] = lynx.asc;
y = log(y);
Nobs = rows(y);


proc LSTAR(theta);
  local beta,sigma,resid,i,logL;

  beta = theta[1:8]; sigma = theta[9];
  resid = miss(zeros(Nobs,1),0);
  LogL = miss(zeros(Nobs,1),0);

  i = 12;
  do until i>Nobs;
    resid[i] = y[i] - beta[1]*y[i-1] -
                 ( beta[2]*y[i-2] + beta[3]*y[i-3]  + beta[4]*y[i-4]
                 + beta[5]*y[i-9] + beta[6]*y[i-11])/
                 ( 1 + exp(-beta[7]*1.8*(y[i-3]-beta[8])) );
    i=i+1;
  endo;

  LogL[12:Nobs] = -0.5*ln(2*pi) -0.5*ln(sigma^2)
                  -0.5*(resid[12:Nobs]^2)/(sigma^2);

  retp(Logl);
endp;


_tsm_parnm  = "y(t-1)"|"y(t-2)"|"y(t-3)"|"y(t-4)"|"y(t-9)"|
             "y(t-11)"|"gamma"|"c"|"sigma";

output file = tdml1a.out reset;

_print = 1;
_tsm_optmum = 1;
__title = "LSTAR model";

sv = ones(9,1);

{theta,stderr,Mcov,Logl} = TD_ml(&LSTAR,sv);

output off;

save lstar1 = theta;



