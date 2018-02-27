new;
library tsm,optmum;
TSMset;

M1 = rndn(4,4);
M2 = rndn(4,4);
M3 = rndn(4,4);
M4 = rndn(4,4);
M5 = rndn(4,4);
M6 = rndn(4,4);
M7 = rndn(4,4);
M8 = rndn(4,4);
M9 = rndn(4,4);

M = zeros(9,16);

i = 1;
do until i > 9;
  s = ftos(i,"M%lf",1,0);
  M[i,.] = vec(varget(s))';
  i = i + 1;
endo;

print "M4";
print M4;
print;
print "xpnd2(M,4)";
print xpnd2(M,4);


