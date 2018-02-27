/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Estimation of a restricted VAR(2) process
**  See LUTKEPOHL, chapter 5
**
*/

new;
library tsm,optmum,pgraph;

cls;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
data = dinv~dinc~dcons;

/*
**  Use the results of the VAR(2) estimation
**  and set non significant estimates equal to 0
*/


w = 1|0|0|0|0|2|0|3|4|0|0|0|0|0|5|0|0|0|0|6|7;
RR = design(w);
r = zeros(21,1);

output file = varx1h.out reset;

{theta1,stderr,Mcov,LOGL1} = varx_CLS(data,1,2,RR,r);
{theta2,stderr,Mcov,LOGL1} = varx_CML(data,1,2,RR,r);

BB1 = reshape(theta1,7,3)';
BB2 = reshape(theta2,7,3)';

print; print "B(ls) = ";
let fmt[1,3] = "*.*lf" 9 4;
call printfm(BB1,1,fmt);
print;

print; print "B(ml) = ";
call printfm(BB2,1,fmt);
print;


output off;

