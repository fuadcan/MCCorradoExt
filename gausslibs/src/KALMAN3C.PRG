new;
library tsm,optmum;
TSMset;

declare external PARAMETERS;

load xt = kal3a; load yt = kal3b;

/* State space model */

proc Z(t);
  local Zt;
  Zt = xt[t]~0~1;
  retp(Zt);
endp;

proc d(t);
  local dt;
  dt = 0;
  retp(dt);
endp;

proc H(t);
  local Ht;
  Ht = 0;
  retp(Ht);
endp;

proc T(t);
  local Tt;
  Tt = (PARAMETERS[2]~-PARAMETERS[3]~0)|zeros(2,3);
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = 0|0|0;
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = (1|1|0)~(0|0|1);
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = diagrv(zeros(2,2),PARAMETERS[4 1]^2);
  retp(Qt);
endp;



a0 = 10|0|0;
P0 = zeros(3,3);

proc ml(theta);
  local Logl;

  PARAMETERS = theta;

  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  call KFiltering(yt,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;


_tsm_Mcov = 1;
_tsm_optmum = 1;

sv = ones(4,1)/2;

output file = kalman3c.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;



