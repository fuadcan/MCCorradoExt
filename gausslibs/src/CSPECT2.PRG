new;
library tsm,optmum,pgraph;
TSMset;

let SIGMA[2,2] = 0.5 0.25 0.25 0.5;

data = RND_arma(0.9|0.3|0.15|0.2,1,0,SIGMA,100|100,250);

x= data[.,1]; y = data[.,2];

_fourier = 1;
{lambda,I1} = PDGM(x);
{lambda,I2} = PDGM(y);
{lambda,I12} = CPDGM(x,y);

_smoothing = {2,2,1,0.25};
I1 = smoothing(I1);
I2 = smoothing(I2);
I12 = smoothing(I12);

{w,alpha,phi} = Cspectrum(I1,I2,I12);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  title("absolute coherency or cross-amplitude spectra");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda,alpha);
  title("phase spectra");
  xy(lambda,phi);


