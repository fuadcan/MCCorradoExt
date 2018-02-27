new;
library tsm,pgraph;

x = recserar(rndn(100,1),0,0.9);
y = recserar(x,0,1)+rndn(100,1);

_fourier = 1;
d1 = fourier(x); d2 = fourier(y);
{dx,dy} = fourier2(x,y);

diff = dy-d2;

graphset;
  _pdate = ""; _pnum = 2;
  begwind;
  window(2,1,0);
  setwind(1);
    title("Difference: real part");
    _ptitlht = 0.30; _paxht = 0.25; _pnumht = 0.25;
    xy(seqa(1,1,100),real(diff));
  nextwind;
    title("Difference: imaginary part");
    xy(seqa(1,1,100),imag(diff));
  endwind;

Nobs = 100;
do until Nobs > 1000;
  x = recserar(rndn(Nobs,1),0,0.9);
  y = recserar(x,0,1)+rndn(Nobs,1);

  _fourier = 1;
  t = hsec;
  d1 = fourier(x); d2 = fourier(y);
  tt = hsec;
  {dx,dy} = fourier2(x,y);
  ttt = hsec;
  print ftos(Nobs,"Number of observations: %lf",5,0);
  print ftos((tt-t)/100,"Computation time (fourier): %lf secondes",10,2);
  print ftos((ttt-tt)/100,"Computation time (fourier2): %lf secondes",10,2);
  print;
  Nobs = Nobs + 200;
endo;

