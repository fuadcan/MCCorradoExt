new;
library tsm,optmum;
TSMset;

#include nw.src;


x = rndn(1000,1);
y = recserar(x,0,1);

sig2 = zeros(16,2);

i = 0;
do until i>15;
  s2 = _nw_sigma2(x,i);
  sig2[i+1,1] = s2;
  s2 = _nw_sigma2(y,i);
  sig2[i+1,2] = s2;
  i = i+1;
endo;

print "            --------------------------------------";
print "             LAGS     white noise          AR1    ";
print "            --------------------------------------";
call printfmt(seqa(0,1,16)~sig2,1~1~1);

