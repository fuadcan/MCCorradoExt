new;
library tsm;
TSMset;

output file = band3.out reset;

format 10,2;

w = seqa(1,1,16);

L = 0;
do until L>6;
  bnd = select(w,L);
  st = ftos(L,"select procedure illustration with L = %lf",1,0);
  print st;
  print bnd;
  print;
  call pause(2);
  L = L + 1;
endo;

output off;
