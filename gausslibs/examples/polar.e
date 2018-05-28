library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 0 15 0 2 2;
wpol = WinOpenPQG( v, "Polar Coordinates", "Polar" );
call WinSetActive( wpol );

#ENDIF

/* create data */
n = seqa(0,1,101);
theta = n.*pi/50;
radius = sin(0.5 .* theta);
theta1 = seqa(0,pi/100,101);
radius1 = sin(theta);
radius2 = cos(theta);

/* create windows, 1 large, 4 small transparent windows */
begwind;
sx = 9/3;
sy = 6.885/3;
makewind(9*.9,6.855*.9,.45,.34,0);
makewind(sx,sy,0,sy*2,1);
makewind(sx,sy,sx*2,sy*2,1);
makewind(sx,sy,0,0,1);
makewind(sx,sy,sx*2,0,1);

_plwidth = 6;
_pltype = 6;
_ptitlht = .2;

#IFUNIX

_pmcolor = 0;

#ENDIF

_pcolor = { 10, 11, 12, 13 };
_pgrid = 2;
fonts("simplex complex");
title("\202GAUSS\l\201Multiple Polar Windows");
polar(radius~radius~radius1~radius1,theta~theta1~theta1~radius2);

nextwind;
fonts("simplex");
_pgrid = 0;
_plwidth = 0;
_pcolor = 10;
title("");
polar(radius,theta);

nextwind;
_pcolor = 11;
polar(radius,theta1);

nextwind;
_pcolor = 12;
polar(radius1,theta1);

nextwind;
_pcolor = 13;
polar(radius1,radius2);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF


