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
*/

new;
library tsm,optmum,pgraph;

/*
**  Simulation of the process with
**
**  var[e(t)] = 0.25, phi1 = 0.95, theta1 = 0.5 and var[u(t)] = 2
**  y(0) = 20;
*/

rndseed 123;

t_ = seqa(1,1,500);

u = rndn(500,1)*sqrt(2);          /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,20,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*sqrt(0.25);
zt = y + e;                 /* Simulate the z(t) process */


Z = {1 0}; d = 0; H = 0.25;
T = {0.95 -0.5,0 0}; c = {0,0}; R = {1,1}; Q = 2;

a0 = {20,0};       /* y(0) = 20 and u(0) = 0 */
P0 = zeros(2,2);

/*  Build the state space form and run the Kalman filter */

call SSM_build(Z,d,H,T,c,R,Q,0);
call KFiltering(zt,a0,P0);

yc = KF_matrix(1);     /* y[t|t-1]      */
v = KF_matrix(2);      /* v[t]          */
a = KF_matrix(3);      /* a[t]          */
P = KF_matrix(4);      /* P[t]          */
ac = KF_matrix(5);     /* a[t|t-1]      */
Pc = KF_matrix(6);     /* P[t|t-1]      */
F = KF_matrix(7);      /* F[t|t-1]      */
invF = KF_matrix(8);   /* F[t|t-1]^(-1) */


/*  yc corresponds to the estimate of z(t|t-1) */

graphset;
  _pdate = ""; _pnum = 2; _pltype = 6|6; _pcolor = 10|13;
  title("ARMA(1,1) process with noise"\
        "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
        "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
  xlabel("Observation");
  _plegstr = "z]t[ process\000Estimated z]t|t-1[ by the Kalman filter";
  _plegctl = {2 5 4 4};
  xy(t_,zt~yc);

/*  v corresponds to the innovations */

title("ARMA(1,1) process with noise"\
      "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
      "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
_plegstr = "Innovations";
_plegctl = 1;
xy(t_,v);

/*  a[.,1] corresponds to the estimate of y(t) */

title("ARMA(1,1) process with noise"\
      "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
      "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
_plegstr = "y]t[ process\000Estimated y]t[ by the Kalman filter";
_plegctl = {2 5 4 4};
xy(t_,y~a[.,1]);

/*  a[.,2] corresponds to the estimate of u(t) */

title("ARMA(1,1) process with noise"\
      "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
      "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
_plegstr = "u]t[ process\000Estimated u]t[ by the Kalman filter";
_plegctl = {2 5 4 4};
xy(t_,u~a[.,2]);

/* ac[.,1] corresponds to the estimate of y(t|t-1) */

title("ARMA(1,1) process with noise"\
      "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
      "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
_plegstr = "y]t[ process\000Estimated y]t|t-1[ by the Kalman filter";
_plegctl = {2 5 4 4};
xy(t_,y~ac[.,1]);

/* ac[.,2] corresponds to the estimate of u(t|t-1) */

title("ARMA(1,1) process with noise"\
      "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
      "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
_plegstr = "u]t[ process\000Estimated u]t|t-1[ by the Kalman filter";
_plegctl = {2 5 4 4};
xy(t_,u~ac[.,2]);


