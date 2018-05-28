/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 89-90
**
**  Kalman Filtering and the Local Level Model
*/

new;
library tsm,optmum,pgraph;

load purse[] = purse.asc;
Nobs = rows(purse);

purse[5] = miss(0,0);      /* Including a missing value */

{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

Z = 1; d = 0; H = beta[2]^2;
T = 1; c = 0; R = 1; Q = beta[1]^2;

call SSM_build(Z,d,H,T,c,R,Q,0);

a_0 = 0;
P_0 = zeros(1,1);

call KFiltering(purse,a_0,P_0);

yc = KF_matrix(1);     /* y[t|t-1]      */
v = KF_matrix(2);      /* v[t]          */
a = KF_matrix(3);      /* a[t]          */
P = KF_matrix(4);      /* P[t]          */
ac = KF_matrix(5);     /* a[t|t-1]      */
Pc = KF_matrix(6);     /* P[t|t-1]      */
F = KF_matrix(7);      /* F[t|t-1]      */
invF = KF_matrix(8);   /* F[t|t-1]^(-1) */

t_ = seqa(1,1,Nobs);

graphset;
  fonts("simplex simgrma");
  _pdate = ""; _pnum = 2;
  begwind;
  makewind(9,6.855,0,0,0);
  makewind(4,3,5,2.5,1);
  setwind(1);
    _pframe = 0;
    title("LOCAL LEVEL MODEL"\
          "\Ly]t[=\202m\201]t[+\202e\201]t[  "\
          "\202m\201]t[=\202m\201]t-1[+\202c\201]t[");
    _plegctl = {2 6 1 4}; _pltype = 6;
    _plegstr = "y]t[\0y]t|t-1[";
    xy(t_,purse~yc);
  nextwind;
    title(""); _pnum = 0; _pframe = 0; _paxes = 0; _pbox = 15;
    _plegctl = {2 9 6 4};
    _plegstr = "\202m\201]t[";
    _pltype = 6; _pcolor = 14;
    xy(t_,a);
  endwind;
