library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps3 = WinOpenPQG( v, "Surface Chart", "Surf3" );
call WinSetActive( wps3 );

#ENDIF

/* setup globals */
_pdate="";
_pframe=0;
_paxes=0;
_psurf= { 0, 0 };

/* create windows */
begwind;
makewind(9/2,6.855/2,0,3.43,1);
makewind(9/2,6.855/2,4.5,3.43,1);
makewind(9/2,6.855/2,0,0,1);
makewind(9/2,6.855/2,4.5,0,1);
makewind(9,6,0,0,1);

/* graph #1 */
_paxes = 1;
_pzclr={ -.3 15, -.28 10, -.2 2, 0. 3, .2 11, .28 15 };
x = seqa(-8.25,.3,72)';
y = seqa(-9.25,.35,72);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt(x^2+y^4);
z = z .* sin(x/3);
ztics(0,.3,.1,1);
volume(1,1,.7);
view(10,16,7);
surface(x,y,z);

/* graph #2 */
nextwind;
graphset;
_paxes = 0;
_pzclr={ -.3 13, -.02 5, 0. 4, .02 12 };
z = z .* sin(x)*.5 .* sin(y)*.5;
volume(1,1,.8);
view(6,16,7);
surface(x,y,z);

/* graph #3 */
nextwind;
_pframe = 0;
_psurf = { 0, 3 };
_pzclr={ -.3 9, -.02 1, 0. 3, .02 11 };
x = seqa(-10.25,.3,72)';
y = seqa(-10.25,.35,72);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt(x^2+y^4);
z = z .* sin(x)*.5 .* (1/cos(y)*.5);
volume(4,4,4);
view(9,12,7);
surface(x,y,z);

/* graph #4 */
nextwind;
_psurf = { 0, 0 };
_pzclr={ -.3 10, -.02 2, 0. 4, .02 12 };
x = seqa(-21.5,.6,71)';
y = seqa(-23.25,.65,71);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt(x^2+y^4);
z = z .* (1/cos(y)) .* (1/cos(x));
z = z .* (1/sin(x/2));
volume(1,1,1);
view(9,16,6);
surface(x,y,z);

/* window #5 uses draw to add text */
setwind(5);
_pmsgctl = {
 2.73 3.43 .4 0 2 9 8,
 2.70 3.45 .4 0 2 15 4};
_pmsgstr = "GAUSS\0GAUSS";
fonts("complex");
draw;

/* end windows, display graph */
endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
