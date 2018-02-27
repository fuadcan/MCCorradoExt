new;
library tsm;
TSMset;

output file = band2.out reset;

format 10,2;

w = seqa(1,1,16);

L = 0;
do until L>6;
  bnd = extract(w,L);
  st = ftos(L,"extract procedure illustration with L = %lf",1,0);
  print st;
  print bnd;
  print;
  call pause(2);
  L = L + 1;
endo;

L = 1|3|5;
bnd = extract(w,L);
st = "extract procedure illustration with L = 1, 3 and 5";
print st;
print bnd;
print;
call pause(2);

L = 2|4|5;
bnd = extract(w,L);
st = "extract procedure illustration with L = 2, 4 and 5";
print st;
print bnd;
print;
call pause(2);

L = 1|2;
bnd = extract(w,L);
st = "extract procedure illustration with L = 1 and 2";
print st;
print bnd;
print;

output off;
