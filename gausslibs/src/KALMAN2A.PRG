new;
library tsm,optmum,pgraph;
declare external PARAMETERS;

cls;

rndseed 156;

x1 = 10+rndn(200,1)*3;
trend = seqa(1,1,200);
constant = 2*ones(100,1)|12*ones(100,1);
noise = rndn(200,1);

y = constant + 3*x1 + trend + noise;

save kal2a = y;
EXOG = ones(200,1)~x1~trend;
save kal2b = EXOG;

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


a0 = zeros(3,1);
P0 = zeros(3,3);

proc ml(theta);
  local Logl;

  PARAMETERS = theta;

  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  call KFiltering(y,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;


_tsm_Mcov = 1;
_tsm_optmum = 1;

sv = 1|0.1*ones(3,1);

output file = kalman2a.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;

PARAMETERS = theta;

call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
call KFiltering(y,a0,P0);

at = KF_matrix(3);                 /* Read a(t) */

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  _pnumht = 0.25; _ptitlht = 0.30; _paxht = 0.25;
  begwind;
  window(3,1,0);
  setwind(1);
    title("\201Evolution of the parameter \202b\201]0,t[");
    xlabel("Observation");
    xy(seqa(1,1,200),at[.,1]);
  setwind(2);
    title("\201Evolution of the parameter \202b\201]1,t[");
    xy(seqa(1,1,200),at[.,2]);
  setwind(3);
   title("\201Evolution of the parameter \202b\201]2,t[");
   xy(seqa(1,1,200),at[.,3]);
  endwind;

