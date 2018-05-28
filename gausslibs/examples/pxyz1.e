library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wxyz1 = WinOpenPQG( v, "XYZ Plot", "XYZ1" );
call WinSetActive( wxyz1 );

#ENDIF

title("XYZ - 3d Curves");
xlabel("X Axis");
ylabel("Y Axis");
zlabel("Z Axis");
_plwidth = { 2 2 0 };
_pframe = 0;
_pcolor = { 13 12 9 };
ztics(0,1,.5,0);
view(5,3,5);
volume(1,1,1);
_plegctl = 1;
_plegstr = "Spherical Sin Wave\0Conical Sin Wave";

t = seqa(0,.0157,400);
a=10; b=1; c=.3;
x = 1.5*(b*(1-c^2*cos(a*t)^2)^.5 .* cos(t));
y = 1.5*(b*(1-c^2*cos(a*t)^2)^.5 .* sin(t));
z = c*cos(a*t) + .5;

a=10; c=.5;
x1 = c*(1+cos(a*t)) .* cos(t) * 1.2;
y1 = c*(1+cos(a*t)) .* sin(t) * 1.2;
z1 = c*(1+cos(a*t)) * .5;
xyz(x~x1,y~y1,z~z1);

#IFUNIX

call WinSetActive( 1 );

#ENDIF

