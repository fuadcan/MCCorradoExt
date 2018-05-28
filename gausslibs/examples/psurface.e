library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps = WinOpenPQG( v, "Surface with Contour Inset","Surface" );
call WinSetActive( wps );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/2.5,6.855/2.5,5.4,0,0);
makewind(2.6,.7,.5,1.5,0);

/* setup globals */
_pdate="";
_pframe=0;
_paxes=0;
_psurf= { 0, 0 };
_pzclr={ 13, 4, 7, 9, 11, 15 };

/* create data */
x = seqa(-10.6,.3,71)';
y = seqa(-12.4,.35,71);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt(x^2+y^4);
z = z .* sin(x/3);

_pmsgctl = { .9 3.8 .12 15 2 14 0 };
_pmsgstr = "3D Reference Lines, Symbols, Arrows";

/* add 3d reference lines */
_pline3d = {
  -2.7   5.7  0  -2.7  -5.7   0 7 6 2,
  -2.7  12.4 .1  -2.7 -12.4  .1 7 6 2,
  -2.7  12.4  0  -2.7  12.4  .1 7 6 2,
  -2.7 -12.4  0  -2.7 -12.4  .1 7 6 2,
 -10.6     0 .1  -2.7     0  .1 7 6 2,
 -10.6     0 .1 -10.6     0   0 7 6 2,
  -2.7     0 .1  -2.7     0   0 7 6 2,
 -10.6     0  0  -9.8     0   0 7 6 2,
 -5.45     0  0  -2.7     0   0 7 6 2};

/* add 3d arrows */
_parrow3 = {
 10.6 -6 .1 10.6  6 .1 5 .2 11 12 6 4,
 10.6 -6 .1    4 -6 .1 5 .2 11 12 6 4 };

/* add 3D symbols */
_psym3d = {
  -2.7  12.4 .1 8 3 11 0,
  -2.7 -12.4 .1 8 3 11 0,
  -2.7     0 .1 8 3 11 0,
 -10.6     0 .1 8 3 11 0};

volume(1,1,1.2);
view(-3,-6,5);
title("Surface with Contour Inset");
surface(x,y,z);

nextwind;
title("");
_pmsgctl=0;
_pbox = 7;
_pline = {
    1  6  -2.7 -12.6 -2.7 12.6 1 7,
    1  6 -10.6   0.0 -2.7  0.0 1 7
    };
_psym = {
    -2.7   12.0 8 5 11 1 0,
    -2.7  -12.0 8 5 11 1 0,
    -2.7      0 8 5 11 1 0,
    -10.6     0 8 5 11 1 0
    };

xtics(-12,12,6,2);
ytics(-12,12,6,0);
ztics(-.3,.3,.008,0);
contour(x,y,z);

setwind(3);
graphset;
_pmsgctl = {
    0.6 4.0 1.75 0 2 15 2,
    1.1 1.5 1.1  0 2 11 0
    };
_pmsgstr = "GAUSS\0Publication Graphics";
_paxes = 0;
_pbox = 15;
fonts("complex");
draw;

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
