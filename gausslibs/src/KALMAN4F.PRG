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
**  Smoothing using the Kalman filter
**
**  If var[e(t)] = 0, z(t) = y(t)
**  then the smoothed estimate of y(t), that is y(t|T), is equal to z(t)
**
**  We can verify this property with an example
*/

new;
library tsm,optmum,pgraph;

/*
**  Simulation of the process with
**
**  var[e(t)] = 0, phi1 = 0.95, theta1 = 0.5 and var[u(t)] = 1
**  y(0) = 40;
*/

rndseed 123;

t_ = seqa(1,1,500);

u = rndn(500,1)*1;          /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,40,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*0;
zt = y + e;                 /* Simulate the z(t) process */

Z = {1 0}; d = 0;
c = {0,0}; R = {1,1};
T = (0.95~-0.5)|(0~0);
Q = 1;
H = 0;
call SSM_build(Z,d,H,T,c,R,Q,0);

/* Initiate the filter with a prior state */

a0 = 40|0; P0 = zeros(2,2);

/* Or use the procedure SSM_ic */

{a0,P0} = SSM_ic;

call KFiltering(zt,a0,P0);

{as,Ps} = KSmoothing;

graphset;
  _pdate = ""; _pnum = 2;
  begwind;
  makewind(9,6.855,0,0,0);
  makewind(9/2,6.855/2,3,2,0);
  setwind(1);
    title("ARMA(1,1) process"\
          "\L \LSMOOTHING");
    xlabel("Observation");
    xy(t_,as[.,1]);
  setwind(2);
    _pnumht = 0.20; _paxht = 0.20; _ptitlht = 0.25;
    title("original data");
    xy(t_,zt);
  endwind;


