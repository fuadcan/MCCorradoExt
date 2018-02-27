library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps5 = WinOpenPQG( v, "Surface Chart", "Surf5" );
call WinSetActive( wps5 );

#ENDIF

_pframe = 0;
_psurf = { 0,0 };
_pdate = "";
_pzclr = 9;
_pnumht = .1;
_paxht = .12;
_ptitlht = .14;
title("Matrix Surface Example");
xlabel("Vector X");
ylabel("Vector Y");
fonts("complex");
margin(0,1.4,0,0);

/* main surface data */
x = seqa(-4,.2,81)';
y = seqa(-8,.2,81);

z = sin(x/5) .* cos(y/5);
z = z + sin(y/2)*.3;
z = z + sin(x/2)*.3;
z = z + sin(x/6)*.5;

/* create hole in surface */
flat = minc(minc(z));
x1 = seqa(-1.0,.04,81)';
y1 = seqa(-3.7,.06,81);
z2 = cos(sqrt((x1/2)^2+(y1/2)^2)) ./ sqrt(x1^2+y1^4);
z2 = -z2 * .25;
z2 = (z2 .> .2) * .2 + (z2 .<= .2 ) .* z2;
z = z + z2';

flat = -.8;
z = (z .< flat) * flat + (z .>= flat ) .* z;
clear x1,y1;

/* draw grid lines on xy plane */
_pline3d={
0 -8 -.8 0 8 -.8 1 6 0 1,
4 -8 -.8 4 8 -.8 1 6 0 1,
8 -8 -.8 8 8 -.8 1 6 0 1,
12 -8 -.8 12 8 -.8 1 6 3 1,
-4 -6 -.8 12 -6 -.8 1 6 0 1,
-4 -4 -.8 12 -4 -.8 1 6 0 1,
-4 -2 -.8 12 -2 -.8 1 6 0 1,
-4  0 -.8 12  0 -.8 1 6 0 1,
-4  2 -.8 12  2 -.8 1 6 0 1,
-4  4 -.8 12  4 -.8 1 6 0 1,
-4  6 -.8 12  6 -.8 1 6 0 1,
-4  8 -.8 12  8 -.8 1 6 3 1
 };

/* draw curve on yz plane */
o = ones(81,1);
y2 = shiftr(y',-1,y[81])';
n = 68;
z1 = z[.,n];
z2 = shiftr(z[.,n]',-1,z[81,n])';
_pline3d = _pline3d | (13*o~y~z1~13*o~y2~z2~o*14~o*6~o*3~o);

_paxes = { 1,1,1};
_pnum = _paxes;
view(-2,-4,2);
ztics(-.8,1.6,.4,0);
volume(1,1,.8);
volume(2.5,2.5,2);
surface(x,y,z);

#IFUNIX

call WinSetActive( 1 );

#ENDIF
