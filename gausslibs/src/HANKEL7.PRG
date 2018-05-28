new;
library tsm,optmum;
TSMset;


cls;

output file = hankel7.out reset;

_svdtol = 1e-5;

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
H = SSM_Hankel(10);
M = rank(H);
print ftos(M,"McMillan order: %lf",1,0);
print; print;


Z = {1 0}; d = 0; H = 0.25;
T = {0.95 -0.5,0 0}; c = {0,0}; R = {1,1}; Q = 2;

call SSM(Z,d,H,T,c,R,Q,0);
call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
H = SSM_Hankel(10);
M = rank(H);
print ftos(M,"McMillan order: %lf",1,0);
print; print;


Z = {1 1 0}; d = 0; H = 0.25;
T = {0.95 -0.5 0,0 0.5 0.25,0.17 0.13 0.16};
c = {0,0,5}; R = {1,1,1}; Q = 2;

call SSM(Z,d,H,T,c,R,Q,0);
call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
H = SSM_Hankel(10);
M = rank(H);
print ftos(M,"McMillan order: %lf",1,0);
print; print;



Z = {1 1 0,0 4 2}; d = 0|10; H = {0.25 0,0 10};
T = {0.95 -0.5 0,0 0.5 0.25,0.17 0.13 0.16};
c = {0,0,5}; R = {1,1,1}; Q = 2;

call SSM(Z,d,H,T,c,R,Q,0);
call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
H = SSM_Hankel(10);
M = rank(H);
print ftos(M,"McMillan order: %lf",1,0);
print; print;



Z = {1 1,4 2,2 -3}; d = 0|0|0; H = {5 1 0,1 4 0,0 0 8};
T = {0.5 0.45,-0.5 0.8};
c = {0,0}; R = eye(2); Q = {2 0.5,0.5 1};

call SSM(Z,d,H,T,c,R,Q,0);
call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */
H = SSM_Hankel(10);
M = rank(H);
print ftos(M,"McMillan order: %lf",1,0);



output off;