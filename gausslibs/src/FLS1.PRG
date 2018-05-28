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


/* Flexible Least Squares with mu = 1 */

mu = 10000|1|1;
{Bfls,u,r2_M,r2_D} = FLS(y,x,mu);

print; print; call pause(2);

{data,retcode} = MISSING(y~x,0);
y_ = data[.,1]; x_ = data[.,2:4];
beta = invpd(x_'x_)*x_'y_;
Bols = beta'.*ones(Nobs,1);

s = seqa(1,1,Nobs);
graphset;
  _pltype = 6|1|3;
  _plegstr = "true\0FLS\0OLS";
  _plegctl = 1;
  _pdate = ""; _pnum = 2; xlabel("Observation");
  title("Flexible Least Squares\LFirst coefficient");
  xy(s,b1~Bfls[.,1]~Bols[.,1]);
  title("Flexible Least Squares\Lsecond coefficient");
  xy(s,b2~Bfls[.,2]~Bols[.,2]);
  title("Flexible Least Squares\LThird coefficient");
  xy(s,b3~Bfls[.,3]~Bols[.,3]);



