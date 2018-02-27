library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wlx = WinOpenPQG( v, "Log-X Graph", "LogX" );
call WinSetActive( wlx );

#ENDIF

/* generate data */
x = seqa(.03,.04,400);
y = sin(x)^2 .* cos(x)^3 ./ (sin(x)+1);
y1 = sin(x/2);
y2 = cos(x/2);

/* setup globals */
_pnumht = .1;
_pnum = 2;
ytics(-1,1,.5,0);
_pltype = { 6, 1, 3 };
_pcolor = { 10, 12, 13 };
_plwidth = 5;
_pgrid = { 2, 1 };
fonts("simplex complex");

/* create windows */
begwind;
makewind(9,6.855,0,0,0);
makewind(9/3.4,6.855/3,1.5,3.8,0);
makewind(9/3.4,6.855/3,1.5,1,0);

title("\202LOGX Graph with Insets\201");
xlabel("Logarithmic X");
_plegctl = 1;
_plegstr = "main\0sin(x/2)\0cos(x/2)";
logx(x,y~y1~y2);

nextwind;
_plegctl = 0;
_pnum = 0;
_pgrid = 0;
_ptitlht = .35;
_plwidth = 0;
_pcolor = 12;
xlabel("");
title("main * sin(x/2)");
margin(.2,.2,0,.5);
fonts("simplex");
logx(x,y.*y1);

nextwind;
_pcolor = 13;
title("main * cos(x/2)");
logx(x,y.*y2);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
