new;
library pgraph,tsm;

rndseed 123;

Nobs = 500;

x = floor(rndu(Nobs,3)*100);

b1 = seqa(0,2*pi/Nobs,Nobs);
b2 = sin(seqa(0,2*pi/Nobs,Nobs));
b3 = cos(seqa(0,2*pi/Nobs,Nobs));
B = b1~b2~b3;

y = sumc(x'.*B');

mu = 100;
{Bfls,u,r2M,r2D} = FLS(y,x,mu);

s = seqa(1,1,Nobs);

graphset;
  _pltype = 6|1|3;
  _plegstr = "true\0FLS";
  _plegctl = 1;
  _pdate = ""; _pnum = 2; xlabel("Observation");
  title("Flexible Least Squares\LFirst coefficient");
  xy(s,b1~Bfls[.,1]);
  title("Flexible Least Squares\Lsecond coefficient");
  xy(s,b2~Bfls[.,2]);
  title("Flexible Least Squares\LThird coefficient");
  xy(s,b3~Bfls[.,3]);




