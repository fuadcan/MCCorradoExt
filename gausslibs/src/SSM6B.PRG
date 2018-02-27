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

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

{Delta1,Ksi1,zeta1} = SSM_orthogonal(1|0|0,8);
{Delta2,Ksi2,zeta2} = SSM_orthogonal(0|1|0,8);
{Delta3,Ksi3,zeta3} = SSM_orthogonal(0|0|1,8);

Delta = Delta1~Delta2~Delta3;
Ksi = Ksi1~Ksi2~Ksi3;
zeta = zeta1~zeta2~zeta3;

z = miss(zeros(3,1),0);

/* Responses to Orthogonal Impulses */

output file = ssm6b.out reset;

mask = ones(1,7);
let fmt[7,3]=
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*le" 10 2
    "*.*lf" 10 0
    "*.*le" 10 2
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



