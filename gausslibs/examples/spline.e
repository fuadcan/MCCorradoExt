library pgraph;
#include pgraph.ext
graphset;

x = seqa(1,1,5);
y = seqa(1,1,6);

z = { 3 3 2 4 5 4,
      3 3 2 3 4 5,
      4 3 2 2 4 4,
      4 4 3 3 5 6,
      5 6 4 4 6 7 };



#ifUNIX
vv = { 0,0,640,480,40,80,1,6,15,0,0,2,2 };
call WinSetActive(WinOpenPQG(vv,"z","z"));
#endif

call surface(x',y,z');



Sigma = 1;
GridFactor = 3;

{ xx,yy,zz } = spline2D(x,y,z,Sigma,GridFactor);

#ifUNIX
vv[1] = 100;
vv[2] = 100;
call WinSetActive(WinOpenPQG(vv,"smoothed z","smoothed z"));
#endif

call surface(xx',yy,zz');

#ifUNIX
call WinSetActive(1);
#endif


