/*
**  Range over standard deviation (R/S statistic)
**
**  Lo [1991], Long-term memory in stock market prices,
**  Econometrica, 59, pages 1279-1313
*/

new;
library tsm,optmum,pgraph;
TSMset;
#include hurst.src;

rndseed 123;

/* fractional process with d = 0.25 */

{data,retcode} = RND_arfima(0.25,0,0,1,1000,250,1);

{H,x} = HURST(data,1,5);

N = x[.,1];
Q = x[.,4];
V = x[.,5];
fit = x[.,6];

print ftos(h,"Hurst exponent: %lf",4,3);
d = H - 0.5;
print; print ftos(d,"Fractional difference: %lf",4,3);


graphset;
  _pdate = ""; _pnum = 2;
  title("R/S statistic");
  xlabel("n"); ylabel("log(R/S)");
  LOGX(N,log(Q)~fit);

  title("V]n[ statistic");
  xlabel("n"); ylabel("V]n[");
  Bound = ones(rows(V),2).*(0.809~1.862);
  _pltype = 6~1~1;
  xy(N,V~Bound);


