/*
**     Min   f(x,y) = x^2 + 5*y^2
**
**     u.c.  x + 2*y = 1
*/

new;
library tsm,optmum;
TSMset;

proc f(c);
  local x,y,val;
  x = c[1];
  y = c[2];
  val = x^2 + 5*y^2;
  retp(val);
endp;

proc g(c);
  local x,y,g;
  x = c[1];
  y = c[2];

  g = 2*x~
      10*y;
  retp(g);
endp;

_op2gdprc = &g;

CC_ = 1~2;
c_ = 1;

{x,fmin,d,retcode} = optprt(optmum2(&f,150|150,CC_,c_));

