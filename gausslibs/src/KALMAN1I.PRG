
new;
library tsm,optmum,pgraph;
TSMset;

/*
** The two important lines are signalled by  <<<<<<<<<<<<<<<<
*/

DECLARE EXTERNAL PARAMETERS;      /*  <<<<<<<<<<<<<<<<<  */

load y = kal1a;
a0 = zeros(2,1);
P0 = zeros(2,2);

proc Z(t);
  local Zt;
  Zt = (PARAMETERS[1]~0)|(0~PARAMETERS[2]);
  retp(Zt);
endp;

proc d(t);
  local dt;
  dt = PARAMETERS[3]|0;
  retp(dt);
endp;

proc H(t);
  local Ht;
  Ht = (PARAMETERS[4]~0)|(0~PARAMETERS[5]);
  retp(Ht);
endp;

proc T(t);
  local Tt;
  Tt = (PARAMETERS[6]~PARAMETERS[7])|(0~PARAMETERS[8]);
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = PARAMETERS[9 10];
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = 1|1;
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = PARAMETERS[11];
  retp(Qt);
endp;


proc ml(theta);
  local Logl;

  PARAMETERS = theta;             /*  <<<<<<<<<<<<<<<<<  */

  call SSM_build(&Z,&d,&H,&T,&c,&R,&Q,1);
  call KFiltering(y,a0,P0);
  Logl = KF_ml;
  retp(Logl);
endp;


let sv = 1 1 10 2 1 0.5 0.3 0.2 1 4 1;

_tsm_Mcov = 3;
_tsm_optmum = 1;

output file = kalman1i.out reset;

{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;





