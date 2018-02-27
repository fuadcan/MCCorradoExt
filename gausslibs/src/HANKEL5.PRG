/*
** Theoretical Hankel matrix of a VAR process
** McMillan order
*/

new;
library tsm,optmum;

let PHI1 = {.5 .1, .4 .5};
let PHI2 = {0 0, .25 0};
beta = vec(PHI1~PHI2);
let SIGMA = {.09 0, 0 .04};

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);
H = SSM_Hankel(10);

M = rank(H);
print;
print "The McMillan order of a VAR(p) process (K-dimensional process) is p*K.";
print;
print ftos(M,"McMillan order: %lf",1,0);





