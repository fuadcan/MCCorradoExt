new;
library tsm,optmum,pgraph;
TSMset;

rndseed 123;
y = rndn(1000,1);

{x,dens,F,retcode} = Kernel(y);

graphset;
  _pdate = ""; _pnum = 2;
  title("Probability density function");
  _plegstr = "estimated\0true";
  _plegctl = 1;
  xy(x,dens~pdfn(x));

  title("Cumulative distribution function");
  xy(x,F~cdfn(x));




