library pgraph;
graphset;

rndseed 263456;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wlog = WinOpenPQG( v, "Log - Log Graph", "LogLog" );
call WinSetActive( wlog );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/2.3,6.855/3,1.3,1.45,0);

n = 500;
x = seqa(.01,.02,n);
r = rndu(n,1)*.5 + .75;
y = ( abs(sin(x)^2 .* cos(x)^3 ./ (sin(x)+1)) ) * 1.0e3;
y = y .* r;

/* detail location... */
_pline ={
1 6 3 100 7 100 1 11 0,
1 6 7 100 7 900 1 11 0,
1 6 7 900 3 900 1 11 0,
1 6 3 900 3 100 1 11 0 };

_pgrid = { 1, 0 };
fonts("simplex complex");
title("\202LOGLOG with Detail Inset\201");
xlabel("Logarithmic X");
ylabel("Logarithmic Y");
loglog(x,y);

nextwind;
graphset;
_pnum  = 0;
_paxes = 0;
_pbox  = 11;
scalx = { 3, 7 };
scaly = { 100, 1000 };
scale(scalx,scaly);
loglog(x,y);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
