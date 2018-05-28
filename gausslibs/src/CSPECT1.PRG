new;
library tsm,pgraph;
TSMset;

data = rndn(100,2);

x = data[.,1];
y = data[.,2];

_fourier = 1;
{lambda,I1} = PDGM(x);
{lambda,I2} = PDGM(y);
{lambda,I12} = CPDGM(x,y);

_smoothing = {2,2,1,0.25};
I1 = smoothing(I1);
I2 = smoothing(I2);
I12 = smoothing(I12);

{w,alpha,phi} = Cspectrum(I1,I2,I12);
q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  title("absolute coherency or cross-amplitude spectra");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda[1:q],alpha[1:q]);
  title("phase spectra");
  xy(lambda[1:q],phi[1:q]);





