new;
library pgraph;

#include icss.src;

x = rndn(100,1)|2.5*rndn(100,1)|rndn(100,1)|0.5*rndn(100,1);
cp = icss(x);
NT = rows(cp);

graphset;
  title("ICSS algorithm");
  _pdate = ""; _pline = {}; _pnum = 2;
  i = 1;
  do until i > NT;
    _pline = _pline | (1~6~cp[i]~-100~cp[i]~100~1~11~1);
    i = i + 1;
  endo;
  xy(seqa(1,1,400),x);


