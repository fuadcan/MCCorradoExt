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

output file = varx2a.out reset;

{theta,stderr,Mcov,LOGL1} = varx_LS(Y,X,1);

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


