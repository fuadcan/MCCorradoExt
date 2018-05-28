new;
library tsm;
TSMset;

output file = band1.out reset;

format 10,2;

w = seqa(1,1,16);

L = 1;
do until L>5;
  bnd = split(w,L);
  st = ftos(L,"split procedure illustration with L = %lf",1,0);
  print st;
  print bnd;
  print;
  call pause(2);
  L = L + 1;
endo;

L = 1|3|4;
bnd = split(w,L);
st = "split procedure illustration with L = 1, 3 and 4";
print st;
print bnd;
print;
call pause(2);

L = 2|4;
bnd = split(w,L);
st = "split procedure illustration with L = 2 and 4";
print st;
print bnd;
print;
call pause(2);

L = 0;
bnd = split(w,L);
st = "split procedure illustration with L = 0";
print st;
print bnd;
print;
call pause(2);

output off;

