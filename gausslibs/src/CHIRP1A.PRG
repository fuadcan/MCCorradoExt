new;
library tsm,pgraph;

M = 11;
N = 2^M;
t = seqa(0,2/N,N);
chirp = sin(100*pi*t^2);

{H,G,Htilde,Gtilde} = Coiflet(2);

_wcenter = 0;
w = wt(chirp,H,G,0);

graphset;
  _pnum = 2; _pdate = ""; _ptitlht = 0.20;
  fonts("simplex simgrma");
  title("Linear chirp --- sin(100\202p\201t[2])");
  xy(t,chirp);

graphset;
    call wplot(w,0,1);
    _pnum = 2; _pdate = ""; _ptitlht = 0.20; _paxht = 0.20;
    fonts("simplex simgrma");
    title("Linear chirp --- sin(100\202p\201t[2])\LWavelet Coefficients");
    ylabel("scale");
    xtics(1,N/2,64,0);
    draw;

y = iwt(w,Htilde,Gtilde,0);

graphset;
  _pnum = 2; _pdate = ""; _ptitlht = 0.15;
  fonts("simplex simgrma");
  title("Difference between the original series and the reconstructed series");
  xy(t,chirp-y);

