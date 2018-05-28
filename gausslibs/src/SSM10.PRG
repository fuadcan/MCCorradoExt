new;
library tsm,optmum,pgraph;
TSMset;

Z = {1 1,4 2,2 -3}; d = 0|0|0; H = {5 1 0,1 4 0,0 0 8};
T = {0.5 0.45,-0.5 0.8};
c = {0,0}; R = eye(2); Q = {2 0.5,0.5 1};

call SSM_build(Z,d,H,T,c,R,Q,0);

Nr = 15;
s = seqa(1,1,Nr);

Omega = SSM_fevd(Nr);

output file = ssm10.out reset;

mask = 1~1~1;
let fmt[3,3]=  "*.*lf" 13 0 "*.*lf" 18 4 "*.*lf" 18 4;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in the first variable (percent)";
print;
print chrs(45*ones(79,1));
print "                         First            Second      ";
print "         periods      innovations       innovations   ";
print chrs(45*ones(79,1));
i = 1;
do until i > Nr;
  rsp = reshape(Omega[i,.],2,3)';
  rsp = rsp[1,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in the second variable (percent)";
print;
print chrs(45*ones(79,1));
print "                         First            Second      ";
print "         periods      innovations       innovations   ";
print chrs(45*ones(79,1));
i = 1;
do until i > Nr;
  rsp = reshape(Omega[i,.],2,3)';
  rsp = rsp[2,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

print; print;

print;
print "FORECAST ERROR VARIANCE DECOMPOSITION in the third variable (percent)";
print;
print chrs(45*ones(79,1));
print "                         First            Second      ";
print "         periods      innovations       innovations   ";
print chrs(45*ones(79,1));
i = 1;
do until i > Nr;
  rsp = reshape(Omega[i,.],2,3)';
  rsp = rsp[1,.];
  call printfm(i~(rsp*100),mask,fmt); print;
  i = i +1;
endo;

output off;
