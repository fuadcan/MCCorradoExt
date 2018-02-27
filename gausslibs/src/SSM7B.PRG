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

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

{Delta1,Ksi1,zeta1} = SSM_orthogonal(1|0,8);
{Delta2,Ksi2,zeta2} = SSM_orthogonal(0|1,8);

Delta = Delta1~Delta2;
Ksi = Ksi1~Ksi2;
zeta = zeta1~zeta2;

z = miss(zeros(2,1),0);

/* Responses to Orthogonal Impulses */

output file = ssm7b.out reset;

mask = ones(1,5);

let fmt[5,3]=
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*lf" 25 0
    "*.*le" 10 2
    "*.*le" 10 2;

i = 0;
do until i > 8;
  w = xpnd2(Delta,1+i)~z~xpnd2(Ksi,1+i);

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

print chrs(45*ones(79,1));
print "Long-term Accumulated responses matrix";
print chrs(45*ones(79,1));
print zeta;

output off;




