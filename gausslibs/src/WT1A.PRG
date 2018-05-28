new;
library pgraph,tsm;

w_ = zeros(1024,1);
t = seqa(1,1,1024);

w = w_; w[5] = 1;
{H,G,Htilde,Gtilde} = Daubechies(4);
_wcenter = 0;
y1a = iwt(w,Htilde,Gtilde,0);
_wcenter = 1;
y1b = iwt(w,Htilde,Gtilde,0);
_wcenter = 2;
y1c = iwt(w,Htilde,Gtilde,0);
y1 = y1a~y1b~y1c;


w = w_; w[30] = 1;
{H,G,Htilde,Gtilde} = Daubechies(20);
_wcenter = 0;
y2a = iwt(w,Htilde,Gtilde,0);
_wcenter = 1;
y2b = iwt(w,Htilde,Gtilde,0);
_wcenter = 2;
y2c = iwt(w,Htilde,Gtilde,0);
y2 = y2a~y2b~y2c;


graphset;
begwind;
window(2,1,0);
xtics(1,1024,128,0);
_ptitlht = 0.30; _paxht = 0.20; _pnumht = 0.20;
setwind(1);
  title("DAUBECHIES #4 e]5[ with different filter delays");
  xy(t,y1);
nextwind;
  title("DAUBECHIES #20 e]30[ with different filter delays");
  xy(t,y2);
endwind;


