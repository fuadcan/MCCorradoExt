new;
library tsm,optmum;
TSMset;

#include kpss.src;

cls;

y = recserar(rndn(300,1),5,1);
{test1,test2} = KPSS(y,8);

print;

y = y - lag1(y);
{test1,test2} = KPSS(y,8);

