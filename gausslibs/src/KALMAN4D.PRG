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
**  Restricted Maximum Likelihood in the time domain
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

u = rndn(500,1)*sqrt(1);    /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,40,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*sqrt(0.25);
zt = y + e;                 /* Simulate the z(t) process */

/*
**  One restriction: THETA1 = 0
*/

RR_ = design(1|0|2|3); r_ = zeros(4,1);

proc sgf_arma_plus_noise(beta,lambda);
  local phi1,theta1,sig_u,sig_e;
  local w,w1,w2,g;
  phi1 = beta[1];
  theta1 = beta[2];
  sig_u = beta[3];
  sig_e = beta[4];
  w = cos(lambda);
  w1 = 1 - 2*theta1*w + theta1^2;
  w2 = 1 - 2*phi1*w + phi1^2;
  g = (w1./w2)*sig_u^2 + sig_e^2;
  retp(g);
endp;

sv = 0.95|sqrt(1|0.25);

output file = kalman4d.out reset;

_tsm_parnm = "phi1"|"theta1"|"sig_e"|"sig_u";
_tsm_Mcov = 1;
_tsm_optmum = 1;
_print = 1;
{theta1,stderr,Mcov,Logl} =  FD_cml(zt,&sgf_arma_plus_noise,sv,RR_,r_);

call pause(2);

Z = {1 0}; d = 0;
c = {0,0}; R = {1,1};
sv = theta1;

T = (sv[1]~-sv[2])|(0~0);
Q = sv[3]^2;
H = sv[4]^2;
call SSM_build(Z,d,H,T,c,R,Q,0);

/* Use the FDML estimates as initial conditions of the filter */

{a0,P0} = SSM_ic;

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

sv = sv[1 3 4];

{theta2,stderr,Mcov,Logl} = TD_cml(&ml,sv,RR_,r_);

output off;



