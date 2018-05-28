new;
library tsm,optmum;
TSMset;

/* Time-invariant SSM */

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM(Z,d,H,T,c,R,Q,0);   /* cn = 0 */


