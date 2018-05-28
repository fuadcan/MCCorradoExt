new;
library tsm,optmum,pgraph;
TSMset;

cls;

{data,retcode} = RND_arfima(0.25,0,0,2,1000,512,1);

{H,G,Htilde,Gtilde} = Haar;

w = wt(data,H,G,0);

pkt = wpkt(data,H,G,0);

_wgrid = 0;

graphset;
    call wpkPlot(pkt,1);
    _pnum = 2; _pdate = ""; _ptitlht = 0.15; _paxht = 0.15;
    title("The wavelet packet table of the fractional process");
    xtics(0,512,50,0);
    draw;

Nobs = rows(w);
M = ln(Nobs)/ln(2);
B = M|seqa(M,-1,M);

wpk = Basis(pkt,B);


graphset;
  _pnum = 2; _pdate = ""; _ptitlht = 0.25; _paxht = 0.25;
  begwind;
  window(2,1,0);
  setwind(1);
    v = ones(M+1,1);
    y1 = maxc(w)*v;
    y2 = minc(w)*v;
    x1 = Nobs*2^(-B);
    x2 = x1;
    _pline = (ones(M+1,1).*(1~3))~x1~y1~x2~y2~(v.*(1~7~0));
    title("The wavelet transform");
    xy(seqa(1,1,Nobs),w);
  nextwind;
    title("The wavelet packet transform with the wavelet basis");
    xy(seqa(1,1,Nobs),wpk);
  endwind;

