

#include eqsolvemt.sdf


proc fct(struct PV par, struct DS d);
    local f1,f2,f3,x;
    x = pvUnpack(par,"x");
    f1 = 3*x[1]^3 + 2*x[2]^2 + 5*x[3] - 10;
    f2 = -x[1]^3 - 3*x[2]^2 + x[3] + 5;
    f3 = 3*x[1]^3 + 2*x[2]^2 -4*x[3];
    retp(f1|f2|f3);
endp;

proc fjc(struct PV par, struct DS d);
    local fjc1,fjc2, fjc3, x;
    x = pvUnpack(par,"x");
    fjc1 = 9*x[1]^2 ~ 4*x[2] ~ 5;
    fjc2 = -3*x[1]^2 ~ -6*x[2] ~ 1;
    fjc3 = 9*x[1]^2 ~ 4*x[2] ~ -4;
    retp(fjc1|fjc2|fjc3);
endp;

start = { -1, 12, -1 };

struct PV par;
par = pvPack(pvCreate,-1|12|-1,"x");

struct eqSolvemtControl c;
c.jacobianProc = &fjc;
c.printIters = 1;

struct eqSolvemtOut out;
out = eqSolvemt(&fct,par,dsCreate,c);



