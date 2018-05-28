cls;

print "qnewton1.e example";
print;
print "This example is taken from G.P.Y. Clarke";
print "Approximate Confidence Limits for a Parameter Function";
print "In a Nonlinear Regression, JASA, 82:221-230, 1987.";
print;
print "The data is the  weight of cut grass from 10 randomly ";
print "sited quadrants taken each week for 13 weeks from a ";
print "grazing pasture.";
print;

print "Press any key to continue"; call keyw;
cls;

   wgt = { 3.183,
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

   nobs = 10;

   week = seqa(1,1,13);

/*  fitting the Mitcherlitz model  */

proc logl(b);
    local dev,s;
    dev = wgt - b[3] - b[2] * exp(-b[1]*week);
    s = dev'dev / rows(dev);
    retp(-sumc(lnpdfmvn(dev,s)));
endp;

start = { 1, 1, 0};

__output = 1;
_qn_PrintIters = 1;

{ bh,fmin,g,retcode } = QNewton(&logl,start);

