new;
library tsm,optmum,pgraph;
TSMset;

_Kernel = {0,.,.,6,0,.};

y = rndchi2(3,1000,1);

{x,dens,F,retcode} = Kernel(y);

graphset;
  _pdate = ""; _pnum = 2;
  title("Probability density function");
  _plegstr = "estimated\0true";
  _plegctl = 1;
  xy(x,dens~pdfchi2(x,3));

  title("Cumulative distribution function");
  xy(x,F~(1-cdfchic(x,3)));


proc (1) = rndchi2(nu,K,L);
  local y,i,u;

  y = zeros(K,L);

  i = 1;
  do until i > L;
    u = rndn(nu,K);
    y[.,i] = sumc(u^2);
    i = i+1;
  endo;

  retp(y);
endp;


proc pdfchi2(x,nu) ;
  retp(exp(-x/2).*(x^(nu/2-1))./((2^(nu/2))*gamma(nu/2))) ;
endp ;


