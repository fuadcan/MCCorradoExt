/*
**  Time-Variant Model
**
**  y(t) = b0*x0(t) + b1*x1(t) + u(t)
**  b0(t) = b0(t-1) + v0(t)
**  b1(t) = b1(t-1) + v1(t)
**
**  Kalman filtering and Kalman smoothing
*/


new;
library tsm,optmum,pgraph;

rndseed 123;

s = seqa(1,1,100);

b0 = recserar(rndn(100,1)*1.5,10,1);
b1 = recserar(rndn(100,1)*0.2,4,1);

x0 = ones(100,1);
x1 = rndu(100,1)*25;

u = rndn(100,1)*2;
X = x0~x1;
Y = X0.*b0 + X1.*b1 + u;

beta = invpd(X'X)*X'Y;

sigma = 2|1.5|0.2;

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


call SSM(&Z,&d,&H,&T,&c,&R,&Q,1);

