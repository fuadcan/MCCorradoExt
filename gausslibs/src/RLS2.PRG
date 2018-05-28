new;
library pgraph,tsm;

rndseed 156;

Nobs = 200;

B1 = ones(Nobs,1)*2;
B2 = recserar(rndn(Nobs,1)*0.5,10,1);
exog = rndu(Nobs,2)*100;
endog = exog[.,1].*B1 + exog[.,2].*B2 + rndn(Nobs,1);


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
  xy(s,B1~(beta[.,1]+bound));


  bound = (-1.96~0~1.96).*stderr[.,2];

  title("Recursive Least Squares");
  ylabel("second coefficient");
  xlabel("observation");
  _plegstr = "True\0Lower CL\0RLS\0Upper CL";
  _plegctl = {2 6 1 4};
  _pltype = 6|3|1|3;
  xy(s,B2~(beta[.,2]+bound));


