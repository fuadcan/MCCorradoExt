new;
library tsm,pgraph;
TSMset;

M = 11;
N = 2^M;
t = seqa(0,2/N,N);
chirp = sin(100*pi*t^2);

{H,G,Htilde,Gtilde} = Coiflet(2);

_wcenter = 0;
pkt = wpkt(chirp,H,G,6);

_wgrid = 0;
_wcolor = seqa(1,1,6);

graphset;
    call wpkPlot(pkt,1);
    _pnum = 2; _pdate = ""; _ptitlht = 0.20; _paxht = 0.20;
    fonts("simplex simgrma");
    title("Linear chirp --- sin(100\202p\201t[2])"\
          "\LWavelet Packet Coefficients");
    xtics(0,N-1,128,0);
    draw;

_wgrid = 1;
_wcolor = 0;

graphset;
    call wpkPlot(pkt,1);
    _pnum = 2; _pdate = ""; _ptitlht = 0.20; _paxht = 0.20;
    fonts("simplex simgrma");
    title("Linear chirp --- sin(100\202p\201t[2])"\
          "\LWavelet Packet Coefficients");
    xtics(0,N/2-1,128,0);
    draw;

graphset;
    call wpkPlot(pkt,1);
    _pnum = 2; _pdate = ""; _ptitlht = 0.20; _paxht = 0.20;
    fonts("simplex simgrma");
    title("Linear chirp --- sin(100\202p\201t[2])"\
          "\LWavelet Packet Coefficients");
    xtics(N/2,N,128,0);
    draw;

B = 1|2|3|3;
wpk = Basis(pkt,B);
cn = wpk.==pkt;

_wgrid = 1;
_wcolor = seqa(1,1,6);

graphset;
    call wpkPlot(pkt,cn);
    _pnum = 2; _pdate = ""; _ptitlht = 0.20; _paxht = 0.20;
    fonts("simplex simgrma");
    title("Linear chirp --- sin(100\202p\201t[2])"\
          "\LWavelet Packet Basis Coefficients");
    xtics(0,N-1,128,0);
    draw;

y = iwpkt(wpk,Htilde,Gtilde,B);

graphset;
  _pnum = 2; _pdate = ""; _ptitlht = 0.15;
  fonts("simplex simgrma");
  title("Difference between the original series and the reconstructed series");
  xy(t,chirp-y);

