new;
library pgraph,tsm;

rndseed 123;

Nobs = 2^10;

t = seqa(0,2*pi/Nobs,Nobs);
x_ = sin(t) + sin(2*t);
x = x_ + rndn(Nobs,1)*0.4;

{H,G,Htilde,Gtilde} = Daubechies(4);

_wcenter = 0;

w_ = wt(x_,H,G,0);
w = wt(x,H,G,0);

lambda = 2.5;
w1 = WaveShrink(w,lambda,"Hard");
w2 = Waveshrink(w,lambda,"Soft");

y1 = iwt(w1,Htilde,Gtilde,0);
y2 = iwt(w2,Htilde,Gtilde,0);

y = y1~y2;


graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.20; _ptitlht = 0.25;
  _wcolor = seqa(1,1,10); _wline = 0; _paxht = 0.25;
  begwind;
  window(2,2,0);
  setwind(1);
    title("Wavelet coefficients (data without noise)");
    call wplot(w_,0,1);
    xtics(0,Nobs/2,32,0);
    draw;
  setwind(2);
    title("Wavelet coefficients (data with noise)");
    call wplot(w,0,1);
    xtics(0,Nobs/2,32,0);
    draw;
  setwind(3);
    title("Wavelet coefficients (Hard Shrinkage)");
    call wplot(w,0,w1);
    xtics(0,Nobs/2,32,0);
    draw;
  setwind(4);
    title("Wavelet coefficients (Soft Shrinkage)");
    call wplot(w,0,w2);
    xtics(0,Nobs/2,32,0);
    draw;
endwind;


graphset;
  fonts("simplex simgrma");
  _pdate = ""; _pnum = 2;
  begwind;
  window(1,2,0);
  setwind(1);
    _pnumht = 0.20;
    _paxht = 0.25;
    _ptitlht = 0.25;
    title("x(t) = sin(t) + sin(2t) + e(t)");
    xlabel("t");
    xy(t,x);
  setwind(2);
    _pltype = 6|1|3|4;
    _plwidth = 5;
    _plegstr = "Hard Shrinkage with \202l\201 = 2.5"\
               "\0Soft Shrinkage with \202l\201 = 2.5"\
               "\0data without noise";
    _plegctl = {2 6 1.8 5.4};
    title("WaveShrink");
    xy(t,y~x_);
  endwind;



