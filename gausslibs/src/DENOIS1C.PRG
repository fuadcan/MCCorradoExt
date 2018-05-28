new;
library pgraph,tsm;

rndseed 123;

Nobs = 2^10;

t = seqa(0,2*pi/Nobs,Nobs);
x_ = sin(t) + sin(2*t);
x = x_ + rndn(Nobs,1)*0.4;

{H,G,Htilde,Gtilde} = Daubechies(4);

_wcenter = 0;
w = wt(x,H,G,0);

w1 = VisuShrink(w,"Hard");
w2 = VisuShrink(w,"Soft");

y1 = iwt(w1,Htilde,Gtilde,0);
y2 = iwt(w2,Htilde,Gtilde,0);

y = y1~y2;

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
    _plegstr = "Visu Hard Shrinkage"\
               "\0Visu Soft Shrinkage"\
               "\0data without noise";
    _plegctl = {2 6 1.8 5.4};
    title("WaveShrink");
    xy(t,y~x_);
  endwind;





