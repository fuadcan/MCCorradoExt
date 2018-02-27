
new;
library tsm,optmum,pgraph;
TSMset;

load y = kal1a;
load ESTIMATES = kal1e;

a0 = zeros(2,1);
P0 = zeros(2,2);


Z = (ESTIMATES[1]~0)|(0~ESTIMATES[2]);
d = ESTIMATES[3]|0;
H = (ESTIMATES[4]~0)|(0~ESTIMATES[5]);
T = (ESTIMATES[6]~ESTIMATES[7])|(0~ESTIMATES[8]);
c = ESTIMATES[9 10];
R = 1|1;
Q = ESTIMATES[11];

call SSM_build(Z,d,H,T,c,R,Q,0);
call KFiltering(y,a0,P0);

{yF,mse,aF,PF} = KForecasting(20);

s = seqa(1,1,120);
Y1 = y[.,1]|miss(zeros(20,1),0);
YF1 = miss(zeros(100,1),0)|YF[.,1];
PF1 = miss(zeros(100,1),0)|PF[.,1];
Bound = YF1 + (-1.96~1.96).*sqrt(PF1);

graphset;
_pdate = ""; _pnum = 2;
title("Forecasting y[1]]t[");
_pltype = 6|1|3|3;
xlabel("Observation");

xy(s,Y1~YF1~Bound);


