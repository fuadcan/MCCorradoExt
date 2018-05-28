library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 40 80 1 6 15 0 0 2 2;
pbw = WinOpenPQG( v, "Bar Chart", "Bar" );
call WinSetActive( pbw );

#ENDIF

begwind;
window(2,1,0);

nbars = 15;
fonts("complex");
t = seqa(0,1,nbars);
x = sin(t/2.5) + .3;
x1 = x * .8;
x2 = x * .6;
x3 = x * .4;
x4 = x * .2;
ylabel("Height");
xlabel("Bar Index");
ytics(-.5,1,.5,1);
_ptitlht = .27;
_paxht = .2;
_pbartyp = { 0 9, 1 10, 3 11, 4 12, 6 13 };
_pnum = 2;

title("Overlapping Bars");
_pbarwid = 0.6;       /* overlapping */
bar(0,x~x1~x2~x3~x4);

nextwind;
title("Side-by-side bars");
_plegctl = { 1 7 13 .7 };
_plegstr = "Accnt 1\0Accnt 2\0Accnt 3\0Accnt 4\0Accnt 5";
_pbarwid = -0.6;      /* side-by-side */
bar(0,x~x1~x2~x3~x4);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF

