new;
library tsm,optmum;

cls;

y = { 2,
      3,
      1,
      7,
      5 };

x = { 1 3 2,
      2 3 1,
      7 1 7,
      5 3 1,
      3 5 5 };

output file = varx3a.out reset;

call ols(0,y,x);
call varx_ls(y,ones(5,1)~x,0);

output off;
