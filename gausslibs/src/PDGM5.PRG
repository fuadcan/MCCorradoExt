new;
library tsm,optmum,pgraph;
TSMset;

rndseed 123456;

_fourier = 0;

Z = eye(2); d = 0|0;
let H[2,2] = 0.2 0 0 0.1;
let T[2,2] = 0.5 0.3 0 -0.5;
c = 0|0; R = 1|1; Q = 0.25;

call SSM_build(Z,d,H,T,c,R,Q,0);
{y,a} = RND_SSM(0|0,100);

{lambda,I} = PDGM2(y);

_smoothing = 15|4|1|0.23;
I1 = smoothing(I[.,1]);
I2 = smoothing(I[.,4]);

G = sgf_SSM(lambda);

q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2;
  _ptitlht = 0.30; _pnumht = 0.30; _paxht = 0.30;
  fonts("simplex simgrma");
  begwind;
  window(1,2,0);
  setwind(1);
    title("first time series");
    xlabel("frequency");
    xtics(0,pi,pi/4,0);
    lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
    asclabel(lab,0);
    _plegstr = "PDGM\0spectral density"; _plegctl = {2 5 3 5};
    xy(lambda[1:q],I1[1:q]~G[1:q,1]);
  nextwind;
    title("second time series");
    xy(lambda[1:q],I2[1:q]~G[1:q,4]);
  endwind;




