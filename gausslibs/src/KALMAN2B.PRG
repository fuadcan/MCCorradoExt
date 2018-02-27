new;
library tsm,optmum,pgraph;
declare external PARAMETERS;


load ENDOG = kal2a;
load EXOG = kal2b;

/* State space model */

proc Z(t);
  local Zt;
  Zt = EXOG[t,.];
  retp(Zt);
endp;

proc d(t);
  local dt;
  dt = 0;
  retp(dt);
endp;

proc H(t);
  local Ht;
  Ht = PARAMETERS[1]^2;
  retp(Ht);
endp;

proc T(t);
  local Tt;
  Tt = eye(3);
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = 0|0|0;
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = eye(3);
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = diagrv(eye(3),PARAMETERS[2:4]^2);
  retp(Qt);
endp;

X = EXOG[1:3,.];
Y = ENDOG[1:3];
beta = invpd(X'X)*X'Y;

a0 = beta;
P0 = zeros(3,3);

Y = miss(zeros(3,1),0)|ENDOG[4:200];

proc ml(theta);
  local Logl;

  PARAMETERS = theta;

  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  call KFiltering(Y,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;


_tsm_Mcov = 1;
_tsm_optmum = 1;

sv = 1|0.1*ones(3,1);

output file = kalman2b.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;


PARAMETERS = theta;

call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
call KFiltering(y,a0,P0);

at = KF_matrix(3);                 /* Read a(t)  */
Pt = KF_matrix(4);                 /* Read P(t); */

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  _pnumht = 0.25; _ptitlht = 0.30; _paxht = 0.25;
  _pltype = 6|1|1; _pcolor = 10|11|11;
  begwind;
  window(3,1,0);
  setwind(1);
    title("\201Evolution of the parameter \202b\201]0,t[");
    xlabel("Observation");

    ESTIMATE = at[.,1];
    STDERR = Pt[.,1];
    Bound = ESTIMATE + (-1.96~1.96).*sqrt(STDERR);
    Bound = ESTIMATE~Bound;
    xy(seqa(1,1,200),Bound);

  setwind(2);
    title("\201Evolution of the parameter \202b\201]1,t[");

    ESTIMATE = at[.,2];
    STDERR = Pt[.,3];
    Bound = ESTIMATE + (-1.96~1.96).*sqrt(STDERR);
    Bound = ESTIMATE~Bound;
    xy(seqa(1,1,200),Bound);

  setwind(3);
    title("\201Evolution of the parameter \202b\201]2,t[");

    ESTIMATE = at[.,3];
    STDERR = Pt[.,6];
    Bound = ESTIMATE + (-1.96~1.96).*sqrt(STDERR);
    Bound = ESTIMATE~Bound;
    xy(seqa(1,1,200),Bound);

  endwind;

