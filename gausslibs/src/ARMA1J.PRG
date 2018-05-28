/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Forecast Error Variance Decomposition
*/

new;
library tsm,optmum,pgraph;

load beta = arma1;
load sigma = sigma1;


call arma_fevd(beta,1,1,SIGMA,8);

E1 ={}; E2 = {};

Nr = 8;  /* Number of responses */

i = 1;
do until i>Nr;
  E = varget("ERROR"$+ftos(i,"%lf",1,0));
  E1 = E1|E[1,.];
  E2 = E2|E[2,.];
  i = i+1;
endo;

i = seqa(1,1,Nr);

output file = arma1j.out reset;

mask = 1~1~1;
let fmt[3,3]=  "*.*lf" 13 0 "*.*lf" 18 4 "*.*lf" 18 4;

print "FORECAST ERROR VARIANCE DECOMPOSITION in Investment (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment      Business Inventories  ";
print "         periods      innovations      innovations   ";
print chrs(45*ones(79,1));
call printfm(i~(E1*100),mask,fmt);

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in Business Inventories (percent)";
print;
print chrs(45*ones(79,1));
print "                       Investment      Business Inventories  ";
print "         periods      innovations      innovations   ";
print chrs(45*ones(79,1));
call printfm(i~(E2*100),mask,fmt);

output off;



