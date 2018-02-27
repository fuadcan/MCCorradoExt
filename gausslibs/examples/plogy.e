library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wly = WinOpenPQG( v, "Log-Y Graph", "LogY" );
call WinSetActive( wly );

#ENDIF

begwind;
makewind(9,6.855,0,0,0);
makewind(9/3,6.855/3,3.0,1.7,1);

/* create data */
x = seqa(.01,.04,250);
y = abs(sin(x)^2 .* cos(x)^2 ./ (sin(x)+1));

title("LOGY Graph w/Polar Inset\labs(sin(x)[2] * cos(x)[2] / (sin(x)+.8))");
xlabel("Sequence of X");
ylabel("Logarithmic Y Scale");

_plwidth = { 0 1 2 4 };
_pltype = 6;
_pcolor = { 11, 3, 9, 1 };
_plegctl = { 1 4 };
_plegstr = "f(x)\0f(x)*.4\0f(x)*.1\0f(x)*.02";
logy(x,y~y*.4~y*.1~y*.02);

nextwind;
_plwidth = 2;
xlabel("");
ylabel("");
_plegctl = 0;
title("Polar Representation");
_ptitlht = .25;
polar(x,y);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
