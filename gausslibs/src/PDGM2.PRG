new;
library tsm,optmum,pgraph;

y = recserar(rndn(500,1),10|10,.5|.3);

_fourier = 0;    /* Use the fast Fourier transform */
{lambda,I} = PDGM(y);
_smoothing = 10|4|1|0.23;  /* Bandwidth = 10 and Tukey lag window */
Is = Smoothing(I);

q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2;
  fonts("simplex simgrma");
  begwind;
  makewind(9,6.855,0,0,0);
  makewind(5,4,3.5,1.5,1);
  setwind(1);
    title("Smoothed periodogram of the AR(2) process"\
          "\LTukey lag window");
    xlabel("frequency");
    xtics(0,pi,pi/4,0);
    lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
    asclabel(lab,0);
    xy(lambda[1:q],Is[1:q]);
  nextwind;
    _ptitlht = 0.20; _pnumht = 0.20; _paxht = 0.20;
    title("Periodogram");
    xy(lambda[1:q],I[1:q]);
  endwind;



