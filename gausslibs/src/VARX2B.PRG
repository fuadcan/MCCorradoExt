/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Estimation of a System of Dynamic Simultaneous Equations
**  See LUTKEPOHL, chapter 10
**
*/

new;
library tsm,optmum,pgraph;

cls;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
INCOME = y[.,2];
CONSUMP = y[.,3];

X = ones(76,1)~lag1(invest);
Y = INCOME~CONSUMP;

/*
**  Estimation of the reduced form (10.2.4)
**  See LUTKEPOHL, page 326
*/

/*
**  We impose the following restrictions:
**
**  The INC(t-1) coefficient in the INC(t) equation is 1
**  The CONS(t-1) coefficient in the INC(t) equation is 0
**  The constant in the INC(t) equation is 0
**  The INV(t-1) coefficient in the CONS(t) equation is 0
*/

w = 0|1|0|2|0|3|4|0;
RR = design(w); r = 1|zeros(7,1);

output file = varx2b.out reset;

{theta,stderr,Mcov,LOGL1} = varx_CLS(Y,X,1,RR,r);

BB = reshape(theta,4,2)';

BB = ("INC(t)"|"CONS(t)")~BB;
print; print;
print "             INC(t-1)   CONS(t-1)    Constant    Inv(t-1)";
let fmt[5,3] = "-*.*s" 8 8
               "*.*lf" 12 4
               "*.*lf" 12 4
               "*.*lf" 12 4
               "*.*lf" 12 4;
call printfm(BB,0~1~1~1~1,fmt);
print;


output off;


