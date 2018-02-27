/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Impulse Response Analysis
**  See LUTKEPOHL, section 3.7.3
**
*/

new;
library tsm,optmum;

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

_print = 0;
{theta,stderr,Mcov,LOGL} = varx_LS(data,1,2);  /* Constant and 2 lags */

beta = theta[1:18];  /* Estimates of the AR part */
SIGMA = _varx_SIGMA;

z = miss(zeros(3,1),0);

/* Responses to Forecast Errors */

output file = varx1e.out reset;

call arma_impulse(beta,2,0,8);

mask = ones(1,7);
let fmt[7,3]=
    "*.*lf" 10 3
    "*.*lf" 10 3
    "*.*lf" 10 3
    "*.*lf" 10 0
    "*.*lf" 10 3
    "*.*lf" 10 3
    "*.*lf" 10 3;

i = 0;
do until i>8;
  x = varget("IMPULSE"$+ftos(i,"%lf",1,0));
  y = varget("_IMPULSE"$+ftos(i,"%lf",1,0));
  w = x~z~y;

  print chrs(45*ones(79,1));
  print ftos(i,"Periods: %lf",1,0);
  print chrs(45*ones(79,1));
  print "        Estimated responses                   Accumulated responses";
  print "            PHI matrix                             PSI matrix      ";
  print chrs(45*ones(79,1));
  call printfm(w,mask,fmt);
  print;
  i = i+1;
endo;


/* Responses to Orthogonal Impulses */

let fmt[7,3]=
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*lf" 10 0
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*le" 10 2;

call arma_orthogonal(beta,2,0,SIGMA,8);

print; print;
i = 0;
do until i>8;
  x = varget("IMPULSE"$+ftos(i,"%lf",1,0));
  y = varget("_IMPULSE"$+ftos(i,"%lf",1,0));
  w = x~z~y;

  print chrs(45*ones(79,1));
  print ftos(i,"Periods: %lf",1,0);
  print chrs(45*ones(79,1));
  print "        Estimated responses                   Accumulated responses";
  print "           THETA matrix                            KSI matrix      ";
  print chrs(45*ones(79,1));
  call printfm(w,mask,fmt);
  print;
  i = i+1;
endo;


output off;

