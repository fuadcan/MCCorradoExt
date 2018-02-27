/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 172-173
**
**  Basic structural model (s = 4)
**
**  MONTE CARLO EXPERIMENTS
*/

new;
library tsm,optmum,pgraph;

rndseed 123456;

sigma2 = {};
Nsimulation = 100;

i = 1;
do until i>Nsimulation;

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

   _tsm_optmum = 1;
   __output = 1;
   _print = 0;
   __title = ftos(i,"Simulation nø %lf",1,0);
   {theta,stderr,Mcov,Logl} = BSM(y,4,sqrt(variance));

   sigma2 = sigma2|(theta'^2);

   i = i+1;
endo;

output file = bsm2.out reset;

cls;
locate 1,1;
print; print;
print "    True values       MC Estimates ";
call printfmt(variance~meanc(sigma2),1~1);
print; print;
output off;


x = {}; d = {};
_Kernel = {0,.,.,6,0,.};
i = 1;
do until i>4;
  {xi,di,Fi,retcode} = Kernel(sigma2[.,i]);
  x = x~xi; d = d~di;
  i = i + 1;
endo;

let string nom = "c" "z" "w" "e"; /* eta zeta omega epsilon */


graphset;
  fonts("simplex simgrma"); _pdate = ""; _pnum = 2;
  _pnumht = 0.20; _ptitlht = 0.25;
  begwind;
  window(2,2,0);
  i = 1;
  do until i > 4;
    setwind(i);
    str = "\201Density of the estimates \202s\201[2]";
    str = str $+"\202]"$+nom[i]$+"[";
    title(str);
    _pline = 1~6~variance[i]~0~variance[i]~1000~1~15~0;
    call xy(x[.,i],d[.,i]);
    i = i+1;
  endo;
  endwind;

