new;
library pgraph,tsm;

rndseed 123;

Nobs = 500;

x = floor(rndu(Nobs,3)*100);

b1 = ones(Nobs,1);
b2 = sin(seqa(0,2*pi/Nobs,Nobs)) + rndn(Nobs,1)*0.25;
b3 = recserar(rndn(Nobs,1),10,0.9);
B = b1~b2~b3;

y = sumc(x'.*B') + rndn(Nobs,1);
y[50] = miss(0,0);

B = eye(3);
C = eye(3);
D = eye(3);
F = eye(3);

_print = 1;
{beta,u,SSE} = GFLS2(y,x,2,1|0,B,C,D,F);

print; print; call pause(2);


{Bfls,u,r2_M,r2_D} = FLS(y,x,1);

print; print; call pause(2);

graphset;
  _pdate = "";
  title("Difference between GFLS and FLS estimates");
  _pltype = 6;
  xy(seqa(1,1,Nobs),beta-Bfls);




