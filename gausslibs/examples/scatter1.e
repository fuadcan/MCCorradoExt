new;
cls;

rndseed 343677;

print "Scatter1 example demonstrates a basic 2 demisional graph";



library pgraph;
x = seqa(5,4,87);
y = rndn(87,1);


_plctrl = -1;
_pstype = 8;
_psymsiz = 1;
xy(x,y);
