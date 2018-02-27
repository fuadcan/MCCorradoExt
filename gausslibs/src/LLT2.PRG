/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 90-93
**
**  Kalman filtering and the Local Linear Trend Model
*/

new;
library tsm,optmum,pgraph;
load gnp[] = gnp.asc;
Nobs = rows(gnp);

_tsm_optmum = 0;   /* Method of scoring */
{beta,stderr,Mcov,Logl} =  sm_LLT(gnp,1|1|1);

Z = 1~0; d = 0;
T = {1 1,0 1}; c = 0|0; R = {1 0,0 1};
theta = beta^2;  /* sigma^2 */
H = theta[3];
Q = (theta[1]~0)|(0~theta[2]);

call SSM_build(Z,d,H,T,c,R,Q,0);

a0 = gnp[1]|0; P0 = zeros(2,2);
call KFiltering(gnp,a0,P0);

a = KF_matrix(3);      /* a[t]          */

t_ = seqa(1909,1,61);

graphset;
  fonts("simplex simgrma");
  _pdate = ""; _pnum = 2;
  begwind;
  makewind(9,6.855,0,0,0);
  makewind(4,3,1,2.5,1);
  setwind(1);
    _pframe = 0;
    title("LOCAL LINEAR TREND MODEL"\
          "\Ly]t[=\202m\201]t[+\202e\201]t[  "\
          "\202m\201]t[=\202m\201]t-1[+\202b\201]t-1[+\202c\201]t[  "\
          "\202b\201]t[=\202b\201]t-1[+\202z\201]t[");
    _plegctl = {2 6 6.5 1};
    _plegstr = "gnp\000\202m\201]t[";
    xy(t_,gnp~a[.,1]);
  nextwind;
      _pnum = 2; _pframe = 1|1; _paxes = 1;
      _ptitlht = 0.25; _pnumht = 0.20;
      title("\202b\201]t[");
      _plegstr = ""; _plegctl = 0;
      xy(t_,a[.,2]);
  endwind;
