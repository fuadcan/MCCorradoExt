/*
**  The Kolmogorov-Smirnov test
**
**  BROCKWELL and DAVIS [1991], Times Series: Theory and Methods,
**  Springer-Verlag, New York, pages 339-342
*/

new;
library tsm,optmum,pgraph;

load purse[] = purse.asc;

_ml_optmum = 0;
{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

{x,Yi,C_bound} = KS_test(purse-lag1(purse),&f,beta);

graphset;
title("Kolmogorov-Smirnov test using the spectral analysis"\
      "\LIs purse a random walk plus noise process?");
_pdate = "";
_pltype = 6|1|1;
scale(1|36,0|1);
xy(x,Yi~C_bound);


proc f(theta,lambda);
  local sigma1,sigma2,g;
  sigma1  = theta[1]; sigma2 = theta[2];
  g = (sigma1^2) + 2*(1-cos(lambda))*(sigma2^2)/(2*pi);
  retp(g);
endp;


proc (3) = KS_test(y,f,theta);
  local f:proc;          /* This is the spectral density */
  local y_,lambda,I,N,q,I_,Yi,x,k_alpha,C_bound;

  y_ = packr(y);
  {lambda,I} = PDGM(y_);
  N = rows(I);
  q = N/2;

  I_ = I[1:q]./f(theta,lambda[1:q]);

  Yi = cumsumc(I_)/sumc(I_);
  x = seqa(1,1,q);

  k_alpha = -1.36~1.36;                          /* 5% level */
  C_bound = (x-1)/(q-1)+k_alpha*(q-1)^(-0.5);
  retp(x,Yi,C_bound);
endp;

