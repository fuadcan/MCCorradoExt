/*
**  Compute the autocovariances, autocorrelations and partial autocorrelations
**  of an ARFIMA process
**
**  Suppose that the process is of the form:
**  (1-0.5L+0.25L^2)(1-L)^d * y(t) = u(t) - 0.3*u(t-1);
**
**  where d takes the different values: -0.5 -0.25 0 0.25 0.5
*/

new;
library tsm,optmum,pgraph;
cls;

d = -0.5|-0.25|0|0.25|0.5;

order = 10;

output file = canon7.out reset;

i=1;
do until i>rows(d);
  beta = 0.5|-0.25|0.3|d[i];
  {phi,theta,retcode} = canonical_arfima(beta,2,1,1000,1000);
  {k,acv,acr,pacr} = AUTOCOV(theta,order);

  locate 1,1; print;
  print ftos(d[i],"                          ARFIMA process with d = %lf",5,3);
  print chrs(45*ones(79,1));
  print "             Order      ACF (cov)        ACF (cor)       PACF (cor)" ;
  print chrs(45*ones(79,1));
  call printfmt(k~acv~acr~pacr,1~1~1~1);

  print; print; call pause(2);
 
  i = i+1;
endo;


output off;

/*
**  Given the Infinite MA representation of a univariate process,
**  AUTOCOV computes the autocovariance, the autocorrelation
**  and the partial autocorrelation functions.
*/

proc (4) = AUTOCOV(theta,order);
  local k,theta_mat,acv,acr,pacr,i,RR,r;

  k = seqa(0,1,order+1);

  theta_mat = shiftr(theta'.*ones(order+1,1),k,zeros(order+1,1));

  acv = theta_mat*theta;  /* Autocovariances  */

  acr = acv/acv[1];       /* Autocorrelations */

  pacr = zeros(order,1);

  pacr[1] = acr[2];

  i = 2;
  do until i > order;
    RR = toeplitz(acr[1:i]);
    r = acr[2:i+1]/RR;
    pacr[i] = r[i];       /* Partial autocorrelations */
    i = i+1;
  endo;

  pacr = miss(0,0)|pacr;

  retp(k,acv,acr,pacr);
endp;

