new;
library tsm,optmum,pgraph;
TSMset;

_Kernel = {0,1,.,6,0,1};

y = rndu(1000,1);

{x,dens,F,retcode} = Kernel(y);

graphset;
  _pdate = ""; _pnum = 2;
  title("Probability density function");
  _plegstr = "estimated\0true";
  _plegctl = 1;
  xy(x,dens~ones(rows(x),1));

  title("Cumulative distribution function");
  xy(x,F~x);

