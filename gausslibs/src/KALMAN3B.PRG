new;
library tsm,optmum;
TSMset;

declare external PARAMETERS;

load xt = kal3b; load yt = kal3a;

/* State space model */

proc Z(t);
  local Zt;
  Zt = xt[t]~0;
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
  Tt = (PARAMETERS[2]~-PARAMETERS[3])|(0~0);
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = 0|0;
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = 1|1;
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = PARAMETERS[4]^2;
  retp(Qt);
endp;



a0 = 10|0;
P0 = zeros(2,2);

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

output file = kalman3b.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;


