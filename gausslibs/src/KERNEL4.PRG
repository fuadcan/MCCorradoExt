/*
**  Vannucci [1995], Nonparametric density estimation using wavelets
**  Discussion Paper 95-26, ISDS, Duke University
*/

new;
library tsm,optmum,pgraph;
TSMset;

load data[] = frfdem.asc;

data = packr(data);

Nobs = rows(data);

dt = 1|2|5|10|30;

x1 = {}; d1 = {};

i = 1;
do until i > 5;
  y = data - lagn(data,dt[i]);
  {x_,dens_,F_,retcode_} = Kernel(y);
  x1 = x1~x_;
  d1 = d1~dens_;
  i = i + 1;
endo;

{H,G,Htilde,Gtilde} = Daubechies(4);

data = padding(data,2);

_wcenter = 2;
w = wt(data,H,G,0);
w = Thresholding(w,0.75);
y = iwt(w,Htilde,Gtilde,0);

data2 = y[1:Nobs];

x2 = {}; d2 = {};

i = 1;
do until i > 5;
  y = data2 - lagn(data2,dt[i]);
  {x_,dens_,F_,retcode_} = Kernel(y);
  x2 = x2~x_;
  d2 = d2~dens_;
  i = i + 1;
endo;

graphset;
  _pdate = ""; _pnum = 2;
  _pnumht = 0.20; _ptitlht = 0.25;
  fonts("simplex simgrma");
  begwind;
  window(2,1,0);
  setwind(1);
    title("Probability density of the return");
    xtics(-0.01,0.01,0.0025,0);
    ytics(0,500,100,0);
    _plegstr = "\202D\201t = 1\0\202D\201t = 2\0\202D\201t = 5"\
               "\0\202D\201t = 10\0\202D\201t = 30";
    _plegctl = {2 6 1 3};
    xy(x1,d1);
  nextwind;
    title("Probability density of the return from the denoised series");
    xtics(-0.01,0.01,0.0025,0);
    ytics(0,500,100,0);
    xy(x2,d2);
  endwind;
