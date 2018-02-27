library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wlx1 = WinOpenPQG( v, "Log-X Graph", "LogX1" );
call WinSetActive( wlx1 );

#ENDIF

x = seqa(.06,.08,500);
y = sin(x)^2 .* cos(x)^3 ./ (sin(x)+1);

title("Logarithmic Graph with XY Inset");
xlabel("Time - sec.");
ylabel("Amplitude ");
_pnumht = .1;
_pnum = 2;
_pgrid = { 2, 2 };
xtics(.04,40,1,0);
ytics(-.8,.8,.4,2);

begwind;
makewind(0,0,0,0,0);
makewind(9/2.8,6.855/3,1.2,1,0);

logx(x,y);

nextwind;
title("");
xlabel("");
ylabel("");
_pnumht = .25;
_pgrid = { 1, 2 };
xtics(0,40,10,2);
ytics(-.8,.8,.4,2);

xy(x,y);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
