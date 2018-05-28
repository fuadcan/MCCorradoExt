new;
library pgraph,tsm;

rndseed 123;

Nobs = 500;
t = seqa(1,1,Nobs);
x = floor(rndu(Nobs,1)*100);
b = recserar(rndn(Nobs,1),0,1)+ sin(pi*t^2);
y = x.*b + rndn(Nobs,1);

_print = 0;

Bfls = zeros(Nobs,6);

mu = 10^4|10^6|10^7|10^8|4*10^8|5*10^9;

i = 1;
do until i > 6;
  {Bfls[.,i],u,r2M,r2D} = FLS(y,x,mu[i]);
  i = i + 1;
endo;

Bols = invpd(x'x)*x'y;

graphset;
  _pdate = ""; _pnum = 2;
  title("Convergence of FLS estimates to OLS estimates"\
         "\Las the mu parameter increases.");
  _pltype = ones(6,1)*3|6;
  _plegctl = 1;
  _plegstr = "mu = 10[4]\0mu = 10[6]\0mu = 10[7]\0mu = 10[8]"\
             "\0mu = 4*10[8]\0mu = 5*10[9]\0ols";
  xy(t,Bfls~Bols*ones(Nobs,1));

