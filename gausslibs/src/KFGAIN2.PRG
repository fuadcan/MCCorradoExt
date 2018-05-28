new;
library optmum,tsm,pgraph;
TSMset;

rndseed 123;

s = seqa(1,1,100);

b0 = recserar(rndn(100,1)*1.5,10,1);
b1 = recserar(rndn(100,1)*0.2,4,1);

b0[50] = miss(0,0);

x0 = ones(100,1);
x1 = rndu(100,1)*25;

u = rndn(100,1)*2;
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

K = KF_gain;

at_1 = KF_matrix(5);
vt = KF_matrix(2);

i = 1;
do until i > 99;

  Kt = reshape(K[i,.],1,2)';
  a = T(i+1)*at_1[i,.]' + c(i+1) + Kt*vt[i,.]';

  print ftos(i,"iteration nø %lf",1,0); print; print a~at_1[i+1,.]'; print; call pause(2);
 
  i = i + 1;

endo;



