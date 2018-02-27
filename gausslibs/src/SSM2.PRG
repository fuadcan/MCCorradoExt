new;
library tsm,optmum;
TSMset;

/* Time-variant SSM */

proc Z(t);
  local Zt;
  Zt = eye(2);
  retp(Zt);
endp;

proc d(t);
  local dt;
  dt = 10|0;
  retp(dt);
endp;

proc H(t);
  local Ht;
  Ht = {2 0, 0 1};
  retp(Ht);
endp;

proc T(t);
  local Tt;
  Tt = {0.5 0.3, 0 0.2};
  retp(Tt);
endp;

proc c(t);
  local ct;
  ct = 1|4;
  retp(ct);
endp;

proc R(t);
  local Rt;
  Rt = 1|1;
  retp(Rt);
endp;

proc Q(t);
  local Qt;
  Qt = 1;
  retp(Qt);
endp;

call SSM(&Z,&d,&H,&T,&c,&R,&Q,1);   /* cn = 1 */
