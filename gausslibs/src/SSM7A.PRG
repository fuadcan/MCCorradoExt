/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Impulse Responses
*/

new;
library tsm,optmum;

load beta = arma1;
load sigma = sigma1;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

{Lambda1,Psi1,zeta1} = SSM_impulse(1|0,8);
{Lambda2,Psi2,zeta2} = SSM_impulse(0|1,8);

Lambda = Lambda1~Lambda2;
Psi = Psi1~Psi2;
zeta = zeta1~zeta2;

z = miss(zeros(2,1),0);

/* Responses to Forecast Errors */

output file = ssm7a.out reset;

mask = ones(1,5);
let fmt[5,3]=
    "*.*lf" 10 3
    "*.*lf" 10 3
    "*.*lf" 25 0
    "*.*lf" 10 3
    "*.*lf" 10 3;

i = 0;
do until i > 8;
  w = xpnd2(Lambda,1+i)~z~xpnd2(Psi,1+i);

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

print chrs(45*ones(79,1));
print "Long-term Accumulated responses matrix";
print chrs(45*ones(79,1));
print zeta;

output off;



