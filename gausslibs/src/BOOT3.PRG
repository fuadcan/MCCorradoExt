/*
** Monte Carlo Bootstrap for state space models
** Example of Stoffer and Wall [1991]
*/

new;
library tsm,optmum,pgraph;
TSMset;

declare external sigma;
declare external f12;
declare external f22;

declare external data;

Nobs = 50;
sigma = sqrt(0.01|0.0025);
f12 = -.85;
f22 = 1.40;
xt = -5 + 10*rndu(Nobs,1);

/* STATE SPACE MODEL */

proc Z(i);
  retp(0~1);
endp;

proc d(i);
  retp(0);
endp;

proc H(i);
  local w;
  w = sigma[1]^2;
  retp(w);
endp;

proc T(i);
  local w;
  w = (0~f12)|(1~f22);
  retp(w);
endp;

proc c(i);
  local w;
  w = 0|(.3*xt[i]);
  retp(w);
endp;

proc R(i);
  retp(0|1);
endp;

proc Q(i);
  local w;
  w = sigma[2]^2;
  retp(w);
endp;

/* STATE SPACE MODEL simulation */

proc simul;
  local i,ai,Qi,nu,Hi,epsilon,yi,y;

  ai = 0|0;
  y = zeros(Nobs,1);

  i = 1;
  do until i > Nobs;

    Qi = Q(i);
    nu = sqrt(Qi)*rndn(1,1);
    ai = T(i)*ai + c(i) + R(i)*nu;
    Hi = H(i);
    epsilon = sqrt(Hi)*rndn(1,1);
    yi = Z(i)*ai + d(i) + epsilon;

    y[i] = yi;
    i = i + 1;
  endo;

  retp(y);
endp;

/* STATE SPACE MODEL log-likelihood */

proc ml(theta);
  local a0,P0,LogL;
  f12 = theta[1];
  f22 = theta[2];
  sigma = theta[3 4];
  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  a0 = 0|0; P0 = zeros(2,2);
  call KFiltering(data,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;

yt = simul;

_tsm_optmum = 1;   /* BFGS algorithm */
_tsm_Mcov = 0;

sv = -.85|1.40|sqrt(0.01|0.0025);
data = yt;
{theta,stderr,Mcov,LogL} = TD_ml(&ml,sv);

f12 = theta[1]; f22 = theta[2]; sigma = theta[3 4];
call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
a0 = 0|0; P0 = zeros(2,2);
call KFiltering(yt,a0,P0);

Nr = 100;  /* Number of the Monte Carlo replications */

ys = zeros(Nobs,Nr);

i = 1;
do until i > Nr;
  {ys[.,i],as} = bootstrap_SSM(a0);
  i = i + 1;
endo;

boot = zeros(Nr,4);

i = 1;
do until i > Nr;
  data = ys[.,i];
  __title = ftos(i,"Monte Carlo Bootstrap -- Iteration nø %lf",2,0);
  {thetas,stderr,Mcov,LogL} = TD_ml(&ml,sv);
  boot[i,.] = thetas';
  i = i + 1;
endo;

graphset;
  title("Bootstrap Histogram for f12");
  _pdate = ""; _pnum = 2; _pgrid = 1|0;
  call hist(boot[.,1],15);

  title("Bootstrap Histogram for f22");
  _pdate = ""; _pnum = 2; _pgrid = 1|0;
  call hist(boot[.,2],15);

