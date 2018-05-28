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
**  Forecast using the Kalman filter
*/

new;
library tsm,optmum,pgraph;

/*
**  Simulation of the process with
**
**  var[e(t)] = 0.25, phi1 = 0.95, theta1 = 0.5 and var[u(t)] = 1
**  y(0) = 40;
*/

rndseed 123;

t_ = seqa(1,1,500);

u = rndn(500,1)*sqrt(1);    /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,40,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*sqrt(0.25);
zt = y + e;                 /* Simulate the z(t) process */


/* You can use the TDML or FDML estimates to forecast the process */

Z = {1 0}; d = 0;
c = {0,0}; R = {1,1};
T = (0.95~-0.5)|(0~0);
Q = 1;
H = 0.25;
call SSM_build(Z,d,H,T,c,R,Q,0);

/* Initiate the filter with a prior state */

a0 = 40|0; P0 = zeros(2,2);

/* Or use the procedure SSM_ic */

{a0,P0} = SSM_ic;

call KFiltering(zt,a0,P0);

{z_forecast,MSE,a_forecast,P_forecast} = KForecasting(10); /* 10 periods */

XX = seqa(451,1,60);
YY = miss(zeros(60,4),0);
YY[1:50,1] = zt[451:500];
YY[51:60,2] = z_forecast;
YY[51:60,3:4] = z_forecast+(-1.96~1.96).*sqrt(MSE);

graphset;
  _pdate = ""; _pnum = 2;
  title("ARMA(1,1) process with noise"\
        "\L \LFORECASTING");
  xlabel("Observation");
  _pltype = 6|1|3|3;
  xy(XX,YY);



