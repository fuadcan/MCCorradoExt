/*
**     Min   f(x,y,z) = x * y * z
**
**     u.c.  x + 2*y -z = 1
*/

new;
library tsm,optmum;
TSMset;

proc fct(c);
  local x,y,z,u;
  x = c[1];
  y = c[2];
  z = c[3];
  u = x*y*z;
  retp(u);
endp;

CC = 1~2~-1;
c = 1;

sv = zeros(3,1);

{x,fmin,g,retcode} = optprt(optmum2(&fct,sv,CC,c));

print;
print "published solution";
print " 1/3      1/6       -1/3";
print;
print "computed solution";
print x';


