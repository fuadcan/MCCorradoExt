
new;
library tsm,optmum;
TSMset;

/* Time-invariant SSM */

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */


{y,a} = RND_SSM(miss(0,0),10);  /* RND_SSM computes the initial values */

mask = ones(1,5);
let fmt[5,3]=
    "*.*lf" 10 0
    "*.*lf" 14 3
    "*.*lf" 14 3
    "*.*lf" 14 3
    "*.*lf" 14 3;

output file = kalman1c.out reset;

j = seqa(1,1,10);
print "      Simulation of the SSM ";
print chrs(45*ones(75,1)); print;
print "  Observation       y1            y2           alpha         beta";
print;
call printfm(j~y~a,mask,fmt);

output off;

rndseed 123456;

{y,a} = RND_SSM(miss(0,0),100);
save kal1a = y;
save kal1b = a;
