/*
** Detection of Changes of Variance of the FRF/DEM exchange rate
*/

new;
library tsm,optmum,pgraph;
TSMset;

#include icss.src;

load x[] = frfdem.asc;

y = Linear_Filter(x,-1|1,-1,0);  /* y(t) = x(t) - x(t-1) */
y = packr(y);

cp = icss(y);
NT = rows(cp);

graphset;
  title("Detection of Changes of Variance of the FRF/DEM exchange rate");
  _pdate = ""; _pline = {}; _pnum = 2;
  ytics(-0.0075,0.0075,0.005,0);
  i = 1;
  do until i > NT;
    _pline = _pline | (1~6~cp[i]~-100~cp[i]~100~1~11~1);
    i = i + 1;
  endo;
  xy(seqa(1,1,rows(y)),y);



