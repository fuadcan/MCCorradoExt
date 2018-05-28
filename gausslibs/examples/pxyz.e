library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wxyz = WinOpenPQG( v, "XYZ Plot", "XYZ" );
call WinSetActive( wxyz );

#ENDIF

begwind;
makewind(9,6.855,0,0,1);
makewind(9/2.9,6.855/2.9,0,0,0);
makewind(9/2.9,6.855/2.9,0,3.8,0);
_psurf = 0;
title("\202XYZ Curve - \201Toroidal Spiral");
fonts("simplex complex");
xlabel("X");
ylabel("Y");
zlabel("Z");

setwind(1);
t = seqa(0,.0157,401);
a = .2; b=.8; c=20;
x = 3*((a*sin(c*t)+b) .* cos(t));
y = 3*((a*sin(c*t)+b) .* sin(t));
z = a*cos(c*t);
margin(.5,0,0,0);
ztics(-.3,.3,.3,0);
_pcolor = 10;
view(-3,-2,4);
volume(1,1,.7);
_plwidth = 5;
xyz(x,y,z);

nextwind;
margin(0,0,0,0);
title("");
x = x .* (sin(z)/10);
_paxes = 0;
_pframe = 0;
_pbox = 13;
_pcolor = 11;
_plwidth = 0;
view(15,2,10);
xyz(x,y,z);

nextwind;
_pcolor = 9;
a = .4; b=.4; c=15;
x = 3*((a*sin(c*t)+b) .* cos(t));
y = 3*((a*sin(c*t)+b) .* sin(t));
z = a*cos(c*t);
volume(1,1,.4);
xyz(x,y,z);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
