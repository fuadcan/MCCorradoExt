/*
** Theoretical and empirical Hankel matrices
*/

new;
library tsm,optmum;

load beta = arma1;
load SIGMA = sigma1;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);
H1 = SSM_Hankel(5);

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

H2 = Hankel(y,5);

format 7,3;

output file = hankel6.out reset;

print "Theoretical Hankel matrix:";
print H1;
print "Empirical Hankel matrix:";
print H2;

output off;

