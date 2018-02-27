/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Impulse Responses
*/

new;
library tsm,optmum,pgraph;

load beta = arma1;
load sigma = sigma1;

call arma_impulse(beta,1,1,8);

z = miss(zeros(2,1),0);

/* Responses to Forecast Errors */

output file = arma1k.out reset;

mask = ones(1,5);
let fmt[5,3]=
    "*.*lf" 10 3
    "*.*lf" 10 3
    "*.*lf" 25 0
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

print; print;

/* Responses to Orthogonal Impulses */

let fmt[5,3]=
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*lf" 25 0
    "*.*le" 10 2
    "*.*le" 10 2;

call arma_orthogonal(beta,1,1,SIGMA,8);

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



