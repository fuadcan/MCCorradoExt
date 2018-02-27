library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps6 = WinOpenPQG( v, "Surface Chart", "Surf6" );
call WinSetActive( wps6 );

#ENDIF

/* setup globals */
@ _pscreen = 0;    tell rerun not to display graph @
_psilent = 1;
_pdate="";
_pframe=0;
_psurf= { 0, 0 };

/* graph #4 (from psurf3.e) */
_pzclr={ -.3 10, -.02 2, 0. 4, .02 12 };
x = seqa(-21.5,.6,71)';
y = seqa(-23.25,.65,71);
z = sin(sqrt((x/2)^2+(y/2)^2)) ./ sqrt(x^2+y^4);
z = z .* (1/cos(y)) .* (1/cos(x));
z = z .* (1/sin(x/2));
ztics(0,100,20,0);
volume(6,6,5);
view(9,16,6);
surface(x,y,z);

#IFUNIX

call WinSetActive( 1 );

#ELSE

#ifDLLCALL
#else
/* now display graph with zoom of 2x centered at x=50% y=40% */
 dos pqgrun graphic -z=2,50,30 /nol;
#endif

#ENDIF
