/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 172-173
**
**  Basic structural model (s = 4)
**
**  Estimation in the frequency domain
*/

new;
library tsm,optmum,pgraph;

rndseed 123;

variance =  1|0.25|1.25|3;
Nobs = 125;
t_ = seqa(1,1,Nobs);

T1 = {1 1,0 1};
T2 = zeros(2,3);
T3 = zeros(3,2);
T4 = {-1 -1 -1,1 0 0,0 1 0};

Z = {1 0 1 0 0}; d = 0; H = variance[4];
T = (T1~T2)|(T3~T4); c = zeros(5,1);
R = eye(3)|zeros(2,3); Q = diagrv(zeros(3,3),variance[1:3]);

call SSM_build(Z,d,H,T,c,R,Q,0);

astar = 100|4|4|2|3;

{y,a} = RND_SSM(astar,Nobs);  /* Simulation of the state space model */

_tsm_optmum = 0;
__output = 1;
_print = 1;

output file = bsm1.out reset;

{theta,stderr,Mcov,Logl} = BSM(y,4,sqrt(variance));

sigma2 = theta^2;

print; print;
print "     True values        Estimates ";
call printfmt(variance~sigma2,1~1);

output off;
