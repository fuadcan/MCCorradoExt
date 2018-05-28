
print "qnewtonmt1.e example";
print;
print "This example is taken from G.P.Y. Clarke";
print "Approximate Confidence Limits for a Parameter Function";
print "In a Nonlinear Regression, JASA, 82:221-230, 1987.";
print;
print "The data is the  weight of cut grass from 10 randomly ";
print "sited quadrants taken each week for 13 weeks from a ";
print "grazing pasture.";
print;

#include qnewtonmt.sdf

struct DS d0;
d0 = reshape(dsCreate,2,1);

d0[1].dataMatrix = { 3.183,
           3.059,
           2.871,
           2.622,
           2.541,
           2.184,
           2.110,
           2.075,
           2.018,
           1.903,
           1.770,
           1.762,
           1.550 };

  d0[2].dataMatrix = seqa(1,1,13);

/*  fitting the Mitcherlitz model  */

proc logl(struct PV par, struct DS d);
    local dev,s,b;
    b = pvUnpack(par,"b");
    dev = d[1].dataMatrix - b[3] - b[2] * exp(-b[1]*d[2].dataMatrix);
    s = dev'dev / rows(dev);
    retp(-lnpdfmvn(dev,s));
endp;

struct QNewtonmtControl c;
c.PrintIters = 0;
c.output = 1;
c.covType = 2;

struct PV par;
par = pvPack(pvCreate,1|1|0,"b");

struct QNewtonmtOut out1;
out1 = QNewtonmt(&logl,par,d0,c);

