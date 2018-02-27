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

_print = 0;

mu = sqrt((seqm(100,5,15)));

r2D = {};
r2M = {};

/* Residual efficiency frontier */

i = 1;
do until i>rows(mu);
  {Bfls,u,_r2M,_r2D} = FLS(y,x,mu[i]);
  r2D = r2D|_r2D;
  r2M = r2M|_r2M;
  i = i+1;
endo;

graphset;
  _pdate = "";
  _pnum = 2;
  title("residual efficiency frontier");
  xlabel("r]D[[2]");
  ylabel("r]M[[2]");
  xy(r2D,r2M);



