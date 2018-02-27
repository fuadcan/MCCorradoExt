new;
library optmum,tsm,pgraph;
TSMset;

Z = {1 1,4 2,2 -3}; d = 0|0|0; H = {5 1 0,1 4 0,0 0 8};
T = {0.5 0.45,-0.5 0.8};
c = {0,0}; R = eye(2); Q = {2 0.5,0.5 1};

call SSM_build(Z,d,H,T,c,R,Q,0);

{a0,P0} = SSM_ic;

{ys,as} = RND_SSM(a0,100);

call KFiltering(ys,a0,P0);

K = KF_gain;

at_1 = KF_matrix(5);
vt = KF_matrix(2);

i = 1;
do until i > 99;

  Kt = reshape(K[i,.],3,2)';
  a = T*at_1[i,.]' + c + Kt*vt[i,.]';

  print ftos(i,"iteration nø %lf",1,0); print; print a~at_1[i+1,.]'; print; call pause(2);
  
  i = i + 1;

endo;


Pbar = ARE(0.00001,100);
Fbar = Z*Pbar*Z'+H;
Kbar = T*Pbar*Z'*invpd(Fbar);   /* Harvey's formula 3.3.12c */

print "The steady-state gain matrix (computed with ARE procedure) is:";
print Kbar;
print;
print "The steady-state gain matrix (computed with KF_gain procedure) is:";
print reshape(K[100,.],3,2)';


