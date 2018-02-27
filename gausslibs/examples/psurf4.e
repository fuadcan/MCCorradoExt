library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wps4 = WinOpenPQG( v, "Surface Chart", "Surf4" );
call WinSetActive( wps4 );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/2.8,6.855/2.8,0,4.3,1);

title("Artificial Mountain Range");
_pdate = "";
_pframe = 0;
_paxes = 0;
_psurf = { 0,0 };
_pzclr = { 1,2,10,11,15 };
_pmcolor = { "","","","","","","","",8 };

x = seqa(-4,.125,161)';
y = seqa(-8,.125,161);
z = sin(x) .* cos(y) * .5;
z = z .* sin(x/3) .* cos(y/3);
z = z .* sin(x/5) + sin(y/2.5)/3 + sin(x/2.5)/3;
_pnum = { 2, 2, 2 };
view(2,-2,3);
surface(x,y,z);

nextwind;
title("");
x = seqa(-4,.125*2,161/2)';
y = seqa(-8,.125*2,161/2);
z = sin(x) .* cos(y) * .5;
z = z .* sin(x/3) .* cos(y/3);
z = z .* sin(x/5) + sin(y/2.5)/3 + sin(x/2.5)/3;
view(2,-2,3);
volume(1,1,.01);
surface(x,y,z);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
