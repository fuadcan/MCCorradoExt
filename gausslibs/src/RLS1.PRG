new;
library pgraph,tsm;

rndseed 156;

Nobs = 100;
Nreg = 2;

B = rndn(Nreg,1)*10;
exog = rndu(Nobs,Nreg)*100;
endog = exog*B + rndn(Nobs,1);


{Beta,stderr,u,RSS} = RLS(endog,exog);

s = seqa(1,1,Nobs);

graphset;
  _pdate = ""; _pnum = 2;

  bound = (-1.96~0~1.96).*stderr[.,1];

  title("Recursive Least Squares");
  ylabel("first coefficient");
  xlabel("observation");
  _plegstr = "True\0Lower CL\0RLS\0Upper CL";
  _plegctl = {2 6 6 4};
  _pltype = 6|3|1|3;
  xy(s,B[1]*ones(Nobs,1)~(beta[.,1]+bound));

