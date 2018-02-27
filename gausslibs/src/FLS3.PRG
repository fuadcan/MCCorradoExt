new;
library pgraph,tsm;

rndseed 123;

Nobs = 500;
S = 250;
z = 2; w = 5;

x = floor(rndu(Nobs,1)*100);
b = ones(S,1)*z|ones(Nobs-S,1)*w;
y = x.*b + rndn(Nobs,1);

_print = 0;
{Bfls1,u,r2_M,r2_D} = FLS(y,x,10^0);
{Bfls2,u,r2_M,r2_D} = FLS(y,x,10^7);
{Brls,stderr,u,rss} = RLS(y,x);
Bols = invpd(x'x)*x'y;

s = seqa(1,1,Nobs);

graphset;
  _pltype = 6;
  fonts("simplex simgrma");
  _plegstr = "true\0FLS  \202m\201 = 1\0FLS \202m\201 = 10[7]\0RLS\0OLS";
  _plegctl = 1;
  _pdate = ""; _pnum = 2; xlabel("Observation");
  title("Unanticipated shift from z to w");
  xy(s,b~Bfls1~Bfls2~(Brls~Bols*ones(Nobs,1)));
