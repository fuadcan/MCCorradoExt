library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wxy2 = WinOpenPQG( v, "XY Plot", "XY2" );
call WinSetActive( wxy2 );

#ENDIF

begwind;
window(2,1,0);
makewind(9/3.3,6.855/3.2,5.62,1.2,0);
makewind(1.5,.38,1.7,5.1,0);
makewind(1.5,.38,1.7,2.1,0);

setwind(1);
n = 500;
x = seqa(.01,.02,n);
y = ( abs(sin(x)^2 .* cos(x)^3 ./ (sin(x)+1)) ) * 1.0e3;
r = rndu(n,1)*.6 + .75;
y1 = y .* r;
tmax = y * 1.2;
tmin = y * .8;

title("\202GAUSS\l\201Input Signal Analysis");
_ptitlht = .28;
ylabel("Amplitude");
xtics(0,10,2,2);
ytics(0,1000,500,5);
_pmsgstr = "Acceptable Tolerance";
_pmsgctl = { 6.5 600 .2 0 1 12 0 };
_parrow = { 6.4 620 5.5 740 3 .15 01 12 1 6 0 };
_pcolor = { 12, 12, 10 };
_pltype = 6;
_pgrid = { 1, 0 };
_pframe = 0;
_pnumht = .21;
_paxht = .21;
_pnum = { 1 2 };
_paxes = { 0 1 };
axmargin(1,.5,1.9,.1);
xy(x,tmin~tmax~y);

nextwind;
axmargin(1,.5,0,.8);
xlabel("Time - milliseconds");
title(" ");
_pmsgctl = 0;
_parrow = 0;
_pcolor = { 12, 12, 11 };
_pltype = 6;
_pnum = 2;
_paxes = 1;
xy(x,tmin~tmax~y1);

nextwind;
graphset;
margin(.5,.5,0,0);
xtics(0,10,2,2);
title("Deviation");
_ptitlht = .37;
_pbartyp = { 6 13 };
_pbox = 15;
_pframe = 0;
_paxes = 0;
_pnum = 0;
bar(x,y1-y);

nextwind;
graphset;
title("Pure Signal");
_ptitlht = 2.2;
_pmcolor = 10;
_pbox = 15;
_paxes = 0;
draw;

nextwind;
title("Distortion");
_pmcolor = 11;
draw;

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
