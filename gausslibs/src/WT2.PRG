new;
library pgraph,tsm,optmum;

rndseed 1234;

Nobs = 2^9;

s = seqa(1,1,Nobs);
sigma = 1;
x0 = miss(0,0);
x = RND_arma(0.62|0.35,2,0,sigma,x0,Nobs);

{H,G,Htilde,Gtilde} = Daubechies(6);

_wcenter = 2;
w = wt(x,H,G,0);

p = 0|0.25|0.75|0.90|0.95|0.99|0.995;

i = 1;
do until i > rows(p);
  w_ = Thresholding(w,p[i]);
  y = iwt(w_,Htilde,Gtilde,0);

  begwind;
  window(2,1,0);
  setwind(1);
    graphset;
    _pdate = ""; _pnum = 2; _pnumht = 0.20; _ptitlht = 0.25; _paxht = 0.20;
    str = ftos(100*p[i],"Wavelet coefficients (Thresholding %lf",4,1);
    str = str $+ "%)";
    title(str);
    _wcolor = seqa(1,1,9); _wline = 0;
    call wplot(w,0,w_);
    xtics(0,Nobs/2,32,0);
    draw;
  setwind(2);
    graphset;
    _pdate = ""; _pnum = 2; _pnumht = 0.20; _ptitlht = 0.25; _paxht = 0.20;
    title("original/reconstructed series");
    xy(s,x~y);
  endwind;

  i = i + 1;

endo;



