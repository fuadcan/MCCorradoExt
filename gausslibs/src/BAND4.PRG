new;
library tsm;
TSMset;

output file = band4.out reset;

format 10,2;

w = seqa(1,1,16);

L = 5;
bnd1 = zeros(1,1);
bnd2 = insert(w,bnd1,L);
print bnd2;
call pause(2);

L = 4;
bnd1 = zeros(1,1);
bnd2 = insert(w,bnd1,L);
print bnd2;
call pause(2);

L = 3;
bnd1 = zeros(2,1);
bnd2 = insert(w,bnd1,L);
print bnd2;
call pause(2);

L = 1;
bnd1 = zeros(8,1);
bnd2 = insert(w,bnd1,L);
print bnd2;
call pause(2);

print;
L = 1;
bnd1 = zeros(9,1);
bnd2 = insert(w,bnd1,L);
print bnd2;

output off;
