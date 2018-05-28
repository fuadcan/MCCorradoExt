library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 40 80 1 6 15 0 0 2 2;
pcw = WinOpenPQG( v, "Contour Map", "Contour" );
call WinSetActive( pcw );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/2.5,6.855/2.5,.5,.5,0);

/* create data */
xl = 51;
yl = 51;
pix = pi/(xl-1);
piy = pi/(yl-1);
x = seqa(-88,2.5,xl)'*pix;
y = seqa(-88,2.5,yl)*piy;
z = sin(x).*cos(x-y);
_plegctl = { 1 3.5 };
_plegstr="";

title("Contour w/Surface Inset\lsin(x)cos(x-y)");
_pxpmax = 3;
_pypmax = 3;
_pzpmax = 3;
_pnum = 0;
ztics(-1.2,1.2,0.07,0);
xtics(-5.6,2.4,1,0);
ytics(-5.6,2.4,1,0);
_pzclr = { 12, 4, 5, 9, 3, 11 };
contour(x,y,z);

nextwind;
_pbox = 15;
_paxes= 0;
_pframe=0;
_plegctl = 0;
_psurf = { 0, 0 };
title("");
xl = 41;
yl = 41;
pix = pi/(xl-1);
piy = pi/(yl-1);
x = seqa(-51,2.5,xl)'*pix;
y = seqa(-51,2.5,yl)*piy;
z = sin(x).*cos(x-y);
margin(0,0,.1,.1);
volume(1,.7,1);
view(3,2.5,2);
surface(x,y,z);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
