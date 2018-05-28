/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Stability analysis of an ARMA process
*/

new;
library tsm,optmum,pgraph;

load beta = arma1;
load sigma = sigma1;

output file = arma1i.out reset;

roots = arma_roots(beta,1,1);  /* ARMA(1,1) process */

print; print chrs(45*ones(79,1));
print "              STABILITY analysis of the Vector ARMA(1,1) process  ";
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


r = real(roots); i = imag(roots);
Nr = rows(roots);
t = seqa(0,2*pi/100,101);
x=cos(t); y=sin(t);

graphset;
  _pdate = ""; _pnum = 2; _pcross = 1; _pframe = {0,0};
  title("Stability of the Vector ARMA(1,1) process"\
      "\LRoots of the reverse characteristic polynomial");
  scale(-3|3,-3|3);
  _psym = r~i~ones(Nr,5).*(11~5~2~1~1);
  xy(x,y);



