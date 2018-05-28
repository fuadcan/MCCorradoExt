/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press
**
**  HARVEY, exercise 4.8 page 233
*/


new;
library tsm,optmum;

Z = 1~0; d = 0;
T = {1 1,0 0}; c = 0|0; R = {1,1};
H = 1; Q = 1;

call SSM_build(Z,d,H,T,c,R,Q,0);

P_bar = ARE(0.00001,100);

a0 = 0|0; P0 = eye(2);
call KFiltering(rndn(100,1),a0,P0);
P = KF_matrix(6);

print "The solution (with ARE procedure) is:";
print P_bar;
print;
print "The solution (with KFiltering procedure) is:";
print xpnd(P[100,.]');

print;
P = KF_matrix(4);
print "Pbar* is:";
print xpnd(P[100,.]');



