library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps2 = WinOpenPQG( v, "Surface Chart", "Surf2" );
call WinSetActive( wps2 );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/3,6.855/3,0,6.855/3*2,0);
makewind(9/3,6.855/3,6,6.855/3*2,0);

_pzpmax = 4;
_pdate  = "";
_paxes  = 1;
_paxht  = .1;
_pnumht = .1;
_pframe = 0;
_psurf  = { 0, 0 };
_pzclr  = { 1, 9, 5, 4, 12, 15 };
_plegctl= { 1 4 };
_plegstr="";
title("Surface and Contour   ");

#ifDOS

xlabel("\212Radial Displacement - X");
ylabel("\212Radial Displacement - Y");

#else

xlabel("Radial Displacement - X");
ylabel("Radial Displacement - Y");

#endif

x = seqa(-9.3,.3,61)';
y = seqa(-12.6,.35,71);
z = sin(sqrt(x^2+y^2)) ./ sqrt(x^2+y^2);
z = z .* sin(x/3);
ztics(0,.3,.1,1);
volume(2.2,2.2,.5);
view(-1,-3,3);
surface(x,y,z);

nextwind;
_pbox = 7;
_paxes = 0;
xlabel("");
ylabel("");
title("");
_plegctl = 0;
ztics(0,1,.5,1);
volume(2,2,.2);
view(-1,-3,6);
surface(x,y,z);

nextwind;
_pnum = 0;
_pbox = 7;
ztics(-.3,.3,.03,0);
xtics(-9,9,3,0);
ytics(-12,12,3,0);
contour(x,y,z);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
