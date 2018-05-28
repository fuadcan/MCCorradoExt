/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Estimation of a VAR(2) process and stability analysis
**  See LUTKEPOHL, section 3.2.3
**
**  Note: the constant is the FIRST column of B in LUTKEPOHL
**        and the LAST column of B in this program
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

output file = varx1a.out reset;

{theta,stderr,Mcov,LOGL1} = varx_LS(data,1,2);  /* Constant and 2 lags */
epsilon = _varx_epsilon;
print;
print "SIGMA:";
print _varx_sigma;


beta = theta[1:18];   /* Estimates of the AR part */
roots = arma_roots(beta,2,0);

print; print chrs(45*ones(79,1));
print "                 STABILITY analysis of the VAR(2) process  ";
print chrs(45*ones(79,1));
print "        Roots of the reverse            Modulus";
print "      characteristic polynomial";
print chrs(45*ones(79,1));
omat = roots~abs(roots);
mask=1~1;
let fmt[2,3]=
    "*.*lf" 14 5
    "*.*lf" 14 5;
call printfm(omat,mask,fmt);

output off;

graphset;

r = real(roots); i = imag(roots);
Nr = rows(roots);
_psym = r~i~ones(Nr,5).*(11~5~2~1~1);
t = seqa(0,2*pi/100,101);
x=cos(t); y=sin(t);

_pdate = "";
title("Stability of the VAR(2) process"\
      "\LRoots of the reverse characteristic polynomial");
scale(-3|3,-3|3);
_pcross = 1;
_pframe = {0,0};
xy(x,y);


