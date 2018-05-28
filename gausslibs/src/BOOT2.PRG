/*
** Bootstrapping a time-variant state space model
** The SSM is a linear model with random coefficients
** (see kalman6a.prg)
*/

new;
library optmum,tsm,pgraph;
TSMset;

Nobs = 100;

b0 = recserar(rndn(Nobs,1)*1.5,10,1);
b1 = recserar(rndn(Nobs,1)*0.2,4,1);

x0 = ones(Nobs,1);
x1 = rndu(Nobs,1)*25;

u = rndn(Nobs,1)*2;
X = x0~x1;
Y = X0.*b0 + X1.*b1 + u;

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

call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
a0 = 10|4;
P0 = zeros(2,2);
call KFiltering(Y,a0,P0);

{ys,as} = bootstrap_SSM(a0);

graphset;
  title("Bootstrapping state space models");
  _pdate = ""; _pnum = 2;
  _plegstr = "Original data\0bootstrap data set";
  _plegctl = 1;
  s = seqa(1,1,Nobs);
  xy(s,y~ys);


