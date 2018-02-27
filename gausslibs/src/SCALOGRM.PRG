/*
**  ARINO and VIDAKOVIC [1995], On wavelet scalograms and their applications
**  in economic time series, DP 95-21, ISDS, Duke University
**
**  ARINO [1995], Time series forecasts via wavelets: an application
**  to car sales in the spanish market, DP 94-30, ISDS, Duke University
**
*/

new;
library pgraph,tsm;

cls;

C1 = intquad1(&gfunc,(8*pi)|0)/(8*pi);
C2 = intquad1(&hfunc,(8*pi)|0)/(8*pi);

t = seqa(1,1,2^12);
yt = gfunc(8*pi*t/4096) - C1;
zt = hfunc(8*pi*t/4096) - C2;
xt = yt+zt;

{H,G,Htilde,Gtilde} = Daubechies(8);
w = wt(xt,H,G,0);

E = scalogram(w);

bnd = split(w,6);

x1 = iwt(bnd[.,1],Htilde,Gtilde,0);
x2 = iwt(bnd[.,2],Htilde,Gtilde,0);

graphset;
begwind;
window(2,2,0);
_ptitlht = 0.25; _pnumht = 0.20;
setwind(1);
  title("data");
  xy(t,xt);
nextwind;
  _pnumht = 0.16;
  lab = "]c]0[[2][ ]R(0)[ ]E(1)[ ]E(2)[ ]E(3)[ ]E(4)[ ]E(5)[ ]E(6)[ "\
        "]E(7)[ ]E(8)[ ]E(9)[ ]E(10)[ ]E(11)[";
  xtics(0,13,1,0);
  asclabel(lab,0);
  title("scalogram");
  bar(0,E);
nextwind;
  graphset; _ptitlht = 0.25; _pnumht = 0.20;
  title("trend (wavelet/true)");
  xy(t,x2~yt);
nextwind;
  title("cycle (wavelet/true)");
  xy(t,x1~zt);
endwind;


proc (1) = gfunc(x);
  local y;
  y = 2*sin(x);
  retp(y);
endp;

proc (1) = hfunc(x);
  local y;
  y = abs(arcsin(sin(16*x/pi)));
  retp(y);
endp;



