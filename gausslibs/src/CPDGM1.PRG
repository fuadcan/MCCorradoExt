new;
library tsm,optmum,pgraph;

y1 = rndn(1000,1);
y2 = recserar(y1,10,0.5);


{lambda,I1} = PDGM(y1);
{lambda,I2} = PDGM(y2);
{lambda,I21} = CPDGM(y2,y1);

_smoothing = 20|2|1|0.23;  /* Bandwidth = 20 and Bartlett lag window */
I1 = Smoothing(I1);
I2 = Smoothing(I2);
I21 = Smoothing(I21);

{w,alpha,phi} = CSpectrum(I1,I2,I21);

q = trunc(rows(lambda)/2);

graphset;
  _pdate=""; _pnum = 2; fonts("simplex simgrma");
  title("Cross-periodogram --- real part");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda[1:q],real(I21[1:q]));

  title("Cross-periodogram --- imaginary part");
  xy(lambda[1:q],imag(I21[1:q]));

  title("Cross-spectrum --- cross-amplitude");
  xy(lambda[1:q],alpha[1:q]);

  title("Cross-spectrum --- phase spectra");
  xy(lambda[1:q],phi[1:q]);




