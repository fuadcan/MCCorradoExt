/*
**     Min   f(x,y,z) = x^2 * y * z
**
**     u.c.  2*x + 4*y + z = 64
*/

new;
library tsm,optmum;
TSMset;

proc fct(c);
  local x,y,z,u;
  x = c[1];
  y = c[2];
  z = c[3];
  u = x^2*y*z;
  retp(-u);
endp;

CC = 2~4~1;
c = 64;

sv = rndn(3,1);

{x0,f0,g0,retc0} = optprt(optmum2(&fct,sv,CC,c));

print;
print "published solution";
print " 16      4       16";
print;
print "computed solution";
print x0';

