new;
library tsm,optmum;


Z = {1 0}; d = 0; H = 0.25;
T = {0.95 -0.5,0 0}; c = {0,0}; R = {1,1}; Q = 2;

call SSM(Z,d,H,T,c,R,Q,0);

