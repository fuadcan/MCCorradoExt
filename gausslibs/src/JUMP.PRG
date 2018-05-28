/*
** The following example is an illustration of jump detection
** by wavelets. The program reproduces the figure 2 of Wang
** [1995], page 392.
**
** Reference:
**
** Wang, Y. [1995], Jump and sharp cusp detection by wavelets,
** Biometrika, 82, 385-397
**
*/


new;
library tsm,optmum,pgraph;
TSMset;

n = 1024;
sigma = 0.2;

proc sign(x);
  local s;
  s = 2*(x.>0)-1;
  retp(s);
endp;

proc f(x);
  local y;
  y = 2*sin(4*pi*x)
      -6*abs(x-0.4)^(3/10)
      -0.5*sign(0.7-x);
  retp(y);
endp;

x = seqa(0,1,n)/n;
t = f(x);                  /* (a) true curve  */
y = t + sigma*rndn(n,1);   /* (b) noisy curve */

{H,G,Htilde,Gtilde} = Daubechies(6);
_wcenter = 0;

w = wt(y,H,G,0);     /* wavelet coefficients */

/* Estimation of sigma */

bnd = select(w,1);     /* finest level coefficients */
bnd = abs(bnd-meanc(bnd));
sigma = median(bnd)/0.6745;

Cn = sigma*sqrt(2*ln(n));

/* selection of a wavelet subband */

bnd = select(w,4);  /* level j = 6, that is 10 - 4 */
bnd = abs(bnd);

graphset;
  begwind;
  window(3,1,0);
  _pdate = ""; _pnum = 2; _paxht = 0.40; _pnumht = 0.30; _ptitlht = 0.30;
  setwind(1);
    title("(a)");
    ylabel("f");
    xtics(0,1,0.2,0);
    ytics(-6.5,-2,1.5,0);
    xy(x,t);
  setwind(2);
    title("(a)");
    ytics(-6.5,-1.5,1,0);
    ylabel("y");
    xy(x,y);
  setwind(3);
    title("(c)");
    ylabel("|TY|");
    _pline = 1~3~0~Cn~1~Cn~0~7~1;
    Nobs = rows(bnd);
    ytics(0,1.2,0.4,0);
    bar(seqa(0,1,Nobs)/Nobs,bnd);
   endwind;

