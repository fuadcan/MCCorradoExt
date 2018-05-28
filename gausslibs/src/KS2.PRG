/*
**  The Kolmogorov-Smirnov test for a unit root process
**
**  BROCKWELL and DAVIS [1991], Times Series: Theory and Methods,
**  Springer-Verlag, New York, pages 339-342
*/

new;
library tsm,optmum,pgraph;

load y[] = frfdem.asc;
T = rows(y);

/*  Scatterplot */

y1 = y[1:T-1]; y2 = y[2:T];

graphset;
  _pdate=""; _pnum = 2;
  title("Scatterplot");
  xlabel("y]t[");
  ylabel("y]t-1[");
  w = ones(T-1,1);
  _pline = 4*w~6*w~y2~y1~0*w~(2*pi)*w~w~15*w~1.5*w;
  m1 = minc(y);
  m2 = maxc(y);
  scale(m1|m2,m1|m2);
  m = seqa(3,0.01,50);
  xy(m,m);


/* If y is integrated of order one, (1-L)y is a white noise process */

dy = packr(y-lag1(y));    /* (1-L)y */

_fourier = 0;
{lambda,I} = PDGM(dy);
N = rows(I);
q = N/2;

I_ = I[1:q];
Yi = cumsumc(I_)/sumc(I_);
x = seqa(1,1,q);


k_alpha = -1.63~1.63;                          /* 1% level */
C_bound = (x-1)/(q-1)+k_alpha*(q-1)^(-0.5);


graphset;
  title("Kolmogorov-Smirnov test using the spectral analysis"\
      "\LTesting for a unit root process");
  _pdate = ""; _pnum = 2;
  _pltype = 6|1|1|2|2;
  scale(0|1000,0|1);
  xy(x,Yi~C_bound);


