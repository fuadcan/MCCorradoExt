/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 172-173
**
**  Basic structural model (s = 4)
**
**  Maximum likelihood in the time domain
*/

new;
library tsm,optmum,pgraph;

rndseed 123;

variance =  1|0.25|0.25|3;
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

seasonality = a[.,3];

output file = BSM3.out reset;

{beta,stderr,Mcov,Logl} = BSM(y,4,sqrt(variance));

/* There are many ways to initialize the filter... */

P0 = zeros(5,5); a0 = y[1]|zeros(4,1);

/*
   Try for example
   sigma2 = beta^3;
   P0 = diagrv(eye(5),sigma2[1 2 3]|sigma2[3]|sigma2[3]);
   a0 = zeros(5,1);
*/

proc ml(theta);
  local H,Q,Logl;
  theta = theta^2;
  H = theta[4];
  Q = diagrv(zeros(3,3),theta[1:3]);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  call KFiltering(y,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 1;

/* Use the FDML estimates as starting values */

{theta,stderr,Mcov,Logl} = TD_ml(&ml,beta);

output off;

sigma2 = theta^2;
H = sigma2[4];
Q = diagrv(zeros(3,3),sigma2[1:3]);
call SSM_build(Z,d,H,T,c,R,Q,0);
call KFiltering(y,a0,P0);
a = KF_matrix(3);

graphset;
  _pdate=""; _pnum = 2;
  fonts("simplex simgrma");
  title("BASIC STRUCTURAL MODEL"\
        "\Ly]t[=\202m\201]t[+\202g\201]t[+\202e\201]t[  "\
        "\202m\201]t[=\202m\201]t-1[+\202b\201]t-1[+\202c\201]t[  "\
        "\202b\201]t[=\202b\201]t-1[+\202z\201]t[   "\
        "\202S\201]0[[3]\202g\201]t-i[=\202w\201]t[");
  _plegstr = "\201True seasonal component\000Estimated seasonal component";
  _plegctl = {2 4.5 1 4.8};
  ytics(-12,14,2,0);
  xtics(0,130,10,0); xlabel("time");
  xy(t_,seasonality~a[.,3]);

