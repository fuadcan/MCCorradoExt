new;
library tsm,optmum;
TSMset;

M1 = complex(rndn(2,2),rndn(2,2));
M2 = complex(rndn(2,2),rndn(2,2));
M3 = complex(rndn(2,2),rndn(2,2));
M4 = complex(rndn(2,2),rndn(2,2));
M5 = complex(rndn(2,2),rndn(2,2));
M6 = complex(rndn(2,2),rndn(2,2));
M7 = complex(rndn(2,2),rndn(2,2));
M8 = complex(rndn(2,2),rndn(2,2));
M9 = complex(rndn(2,2),rndn(2,2));

M = zeros(9,4);

/*
** WARNING: The bookkeeping transpose operator has to be used
**          for complex matrices
*/

i = 1;
do until i > 9;
  s = ftos(i,"M%lf",1,0);
  M[i,.] = vec(varget(s)).';
  i = i + 1;
endo;

print "M4";
print M4;
print;
print "xpnd2(M,4)";
print xpnd2(M,4);



