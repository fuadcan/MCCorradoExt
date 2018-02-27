new;
library tsm,pgraph,optmum;

Nobs = 200;
s = seqa(1,1,Nobs);

x = recserar(rndn(Nobs,1)*2,0,0.5);
y = surrogate(x);

graphset;
  _pdate = ""; _pnum = 2;
  title("Surrogate data technique");
  xlabel("Observation");
  _plegstr = "Original data\0Surrogate data";
  _plegctl = 1;
  xy(s,x~y);

{lambda,Ix} = PDGM(x);
{lambda,Iy} = PDGM(y);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  title("Surrogate data technique:\LThe two periodogram are the same");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda,Ix~Iy);

Rx = autoc(x,10);
Ry = autoc(y,10);

graphset;
  _pdate = ""; _pnum = 2;
  title("Surrogate data technique");
  _plegstr = "covariance function of the original data"\
             "\0covariance function of the surrogate data";
  _plegctl = {2 5 3 4};
  xlabel("lag");
  xy(seqa(0,1,11),Rx~Ry);


proc autoc(x,k);
  local t,rho;
  x=packr(x); t=rows(x); x=x-meanc(x);
  rho=rev(conv(x,rev(x),t-k,t));
  retp(rho/t);
endp;


