/*
**  Time-Variant Model
**
**  y(t) = b0*x0(t) + b1*x1(t) + u(t)
**  b0(t) = b0(t-1) + v0(t)
**  b1(t) = b1(t-1) + v1(t)
**
**  Maximum Likelihood Estimation
*/


new;
library tsm,optmum,pgraph;

declare external sigma;

rndseed 123;

s = seqa(1,1,100);

b0 = recserar(rndn(100,1)*1.5,10,1);
b1 = recserar(rndn(100,1)*0.2,4,1);

x0 = ones(100,1);
x1 = rndu(100,1)*25;

u=rndn(100,1)*2;
X = x0~x1;
Y = X0.*b0 + X1.*b1 + u;

proc Z(i);
  local w;
  w = X[i,.];
  retp(w);
endp;

proc d(i);
  local w;
  w = 0;
  retp(w);
endp;

proc T(i);
  local w;
  w = eye(2);
  retp(w);
endp;

proc c(i);
  local w;
  w = 0|0;
  retp(w);
endp;

proc R(i);
  local w;
  w = eye(2);
  retp(w);
endp;


proc H(i);
  local w;
  w = sigma[1]^2;
  retp(w);
endp;

proc Q(i);
  local w;
  w = eye(2).*sigma[2:3]';
  w = w^2;
  retp(w);
endp;

proc ml(theta);
  local a0,P0,LogL;
  sigma = theta[1:3];
  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  a0 = b0[1]|b1[1];
  P0 = zeros(2,2);
  call KFiltering(Y,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;

_tsm_optmum = 0;   /* BHHH algorithm */
_tsm_parnm = "sig_u"|"sigma0"|"sigma1";

output file = kalman6b.out reset;

{theta,stderr,Mcov,LogL} = TD_ml(&ml,2|1.5|0.2);

output off;

beta = invpd(X'X)*X'Y;
sigma = theta[1:3];
a0 = b0[1]|b1[1];
P0 = zeros(2,2);

call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
call KFiltering(Y,a0,P0);

yc = KF_matrix(1);     /* y[t|t-1]      */
v = KF_matrix(2);      /* v[t]          */
a = KF_matrix(3);      /* a[t]          */
P = KF_matrix(4);      /* P[t]          */
ac = KF_matrix(5);     /* a[t|t-1]      */
Pc = KF_matrix(6);     /* P[t|t-1]      */
F = KF_matrix(7);      /* F[t|t-1]      */
invF = KF_matrix(8);   /* F[t|t-1]^(-1) */

{as,Ps} = Ksmoothing;

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  xlabel("Observation");
  _plegstr = "\202b\201]0[\000Kalman filter\000Kalman smoothing\000OLS";
  _plegctl = {2 6 5 1};
  _pltype = 6|1|3|2;
  xy(s,b0~a[.,1]~as[.,1]~(beta[1]*ones(100,1)));

  _plegstr = "\202b\201]1[\000Kalman filter\000Kalman smoothing\000OLS";
  xy(s,b1~a[.,2]~as[.,2]~(beta[2]*ones(100,1)));

