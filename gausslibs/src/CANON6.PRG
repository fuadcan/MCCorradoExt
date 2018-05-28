/*
**  Compute the Impulse Responses and the Accumulated Impulse Responses
**  of an ARFIMA process
**
**  Suppose that the process is of the form:
**
**  (1-0.5L+0.25L^2)(1-L)^d * y(t) = u(t) - 0.3*u(t-1);
**
**  where d takes the different values: -0.5 -0.25 0 0.25 0.5
*/

new;
library tsm,optmum,pgraph;

d = -0.5|-0.25|0|0.25|0.5;

Nresp = 20;
Period = seqa(1,1,Nresp);

phi = {}; theta = {};

i=1;
do until i>5;
  beta = 0.5|-0.25|0.3|d[i];
  {phi_,theta_,retcode} = canonical_arfima(beta,2,1,1000,Nresp-1);
  phi = phi~phi_; theta = theta~theta_;
  i = i+1;
endo;

graphset;
  _pdate = ""; _pnum = 2;
  str = "\L(1-0.5L+0.25L[2])(1-L)[d]y]t[ = u]t[-0.3u]t-1[";
  title("Impulse Responses of the ARFIMA process\L" $+ str );
  _plegstr = "d=-0.5\0d=-0.25\0d=0\0d=0.25\0d=0.50";
  _plegctl = {2 4 5 3.5};
  xlabel("periods");
  xtics(1,Nresp,1,0);
  xy(Period,theta[1:Nresp,.]);

  title("Accumulated Impulse Responses of the ARFIMA process\L" $+ str );
  ylabel("interim multipliers");
  _plegctl = {2 4 1.75 3.75};
  xy(Period,cumsumc(theta[1:Nresp,.]));

