new;
cls;
library pgraph;
graphset;

rndseed 3456456;

print "Scatter2.e example demonstrates a basic 3d scatter plot";
print "using random number 3d matrix";

x = rndn(100,10);
y = rndn(100,10);
z = rndn(100,10);


_plctrl = -1;
_pstype = 8;
_psymsiz = 1;
xyz(x,y,z);
