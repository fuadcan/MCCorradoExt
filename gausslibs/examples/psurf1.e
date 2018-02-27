library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps1 = WinOpenPQG( v, "Surface Chart", "Surf1" );
call WinSetActive( wps1 );

#ENDIF

/* setup globals for all windows */
_pdate="";
_pframe=0;
_paxes=0;
_psurf= { 0, 0 };

/* create windows */
begwind;
makewind(9,6.855,0,0,0);
makewind(9/2.2,6.855/2.2,0,3.42,1);
makewind(9/2.2,6.855/2.2,5.2,3.46,1);

/* graph #1 */
_pzclr={ -.3 2, 0. 10, .06 14 };
x = seqa(-42.25,1.2,72)';
y = seqa(-42.25,1.2,72);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt((x/2)^8+y^2);
volume(1,1,1);
view(10,-10,10);
surface(x,y,z);

/* graph #2 */
nextwind;
_psurf = { 0, 3 };
_pzclr={ -.3 4, .9 12, 1.7 13 };
x = seqa(-3.05,.06,102)';
y = seqa(-.95,.02,42);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt((x/2)^8+y^2);
volume(1,1,1);
view(7,10,9);
surface(x,y,z);

/* graph #3 */
nextwind;
x = seqa(-3.25,.1,72)';
y = seqa(-3.25,.1,72);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt((x/2)^8+y^2);
_paxes = 1;
_psurf = { 0, 0 };
_pzclr={ -.3 1, 0.3 9, .5 11 };

/* Build additional axes and grid */
_pline3d = {
-3.25 -3.25 7     4 -3.25 7 15 6 2 1,
-3.25 -3.25 0 -3.25 -3.25 7 15 6 2 1,
-3.25 -3.25 1     4 -3.25 1  1 6 0 1,
-3.25 -3.25 2     4 -3.25 2  1 6 0 1,
-3.25 -3.25 3     4 -3.25 3  1 6 0 1,
-3.25 -3.25 4     4 -3.25 4  1 6 0 1,
-3.25 -3.25 5     4 -3.25 5  1 6 0 1,
-3.25 -3.25 6     4 -3.25 6  1 6 0 1 };

/* Build extra X,Z curve */
o = ones(rows(z),1);
z1 = z[33,.]';
z2 = shiftr(z[33,.],-1,0)';
x2 = shiftr(x,-1,0)';
_pline3d=
_pline3d | (x'~o*-3.25~z1~x2+.01~o*-3.25~z2~o*12~o*6~o*0~o);

/* Build extra Y,Z curve */
z1 = z[.,25];
z2 = shiftr(z[.,25]',-1,0)';
y2 = shiftr(y',-1,4)';
_pline3d=
_pline3d | (o*-3.25~y~z1~o*-3.25~y2~z2~o*12~o*6~o*0~o);
clear z2,y2,z1,z2;

xlabel("X");
ylabel("Y");
zlabel("Z");
volume(1,1,1);
view(6,14,7);
surface(x,y,z);

/* end windows, display graph */
endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
