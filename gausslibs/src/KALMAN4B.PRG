/*
**  State space model of an ARMA(1,1) process with noise
**
**  Suppose that
**     z(t) = y(t) + e(t)
**  and
**     y(t) = phi1*y(t-1) + u(t) -theta1*u(t-1)
**
**  e(t) is the noise (the measure error for example)
**  y(t) is the ARMA(1,1) process
**  z(t) is the ARMA(1,1) process with noise
**
**  The state space form is:
**                       _      _
**              _      _ | y(t) |
**     z(t)  =  | 1  0 | |      |  + e(t)
**              -      - | u(t) |
**                       -      -
**
**    _      _    _               _  _        _     _   _
**    | y(t) |    | phi1  -theta1 |  | y(t-1) |     | 1 |
**    |      |  = |               |  |        |  +  |   | u(t)
**    | u(t) |    |   0       0   |  | u(t-1) |     | 1 |
**    -      -    -               -  -        -     -   -
**
**  Conditional Maximum Likelihood with the Kalman Filter
*/

new;
library tsm,optmum,pgraph;

/*
**  Simulation of the process with
**
**  var[e(t)] = 0.25, phi1 = 0.95, theta1 = 0.5 and var[u(t)] = 1
**  y(0) = 40;
*/

rndseed 123456;

t_ = seqa(1,1,500);

u = rndn(500,1)*sqrt(1);     /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,40,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*sqrt(0.25);
zt = y + e;                 /* Simulate the z(t) process */


Z = {1 0}; d = 0;
c = {0,0}; R = {1,1};

a0 = {40,0};       /* y(0) = 40 and u(0) = 0 */
P0 = zeros(2,2);

proc ml(theta);
  local T,H,Q,LogL;
  T = (theta[1]~-theta[2])|(0~0);
  Q = theta[3]^2;
  H = theta[4]^2;
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(zt,a0,P0);
  LogL = KF_ml;
  retp(LogL);
endp;

_print = 1;
_tsm_optmum = 1;
_tsm_Mcov = 1;
_tsm_parnm = "phi1"|"theta1"|"sig_e"|"sig_u";

output file = kalman4b.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,0.95|0.5|sqrt(1|0.25));

output off;



