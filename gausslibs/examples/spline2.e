library pgraph;
#include pgraph.ext
graphset;


rndseed 456356475;

x = seqa(1,1,5);
y = 1 + .2 * x + .2 * rndn(rows(x),1);

#ifUNIX
vv = { 0,0,640,480,40,80,1,6,15,0,0,2,2 };
call WinSetActive(WinOpenPQG(vv,"",""));
#endif

Sigma = 1;
GridFactor = 5;
s = 0.01;
d = 1;


{ xx,yy } = spline1D(x,y,d,s,Sigma,GridFactor);



dd = miss(zeros(rows(xx)-rows(x),1),0);

x = x | dd;
y = y | dd;

xy(xx~x,yy~y);

#ifUNIX
call WinSetActive(1);
#endif


