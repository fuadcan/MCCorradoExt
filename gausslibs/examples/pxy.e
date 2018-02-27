library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wxy = WinOpenPQG(v,"XY Plot","XY");
call WinSetActive(wxy);

#ENDIF

/* setup windows */
begwind;
makewind(9,6.855,0,0,0);    @ Main graph window @
makewind(4,.7,.8,2.5,0);     @ Text Window       @

setwind(1);
/* generate data */
x = seqa(0,.06,250);
y = sin(x)^2 .* cos(x)^3 ./ (sin(x)+1);
y = y ~ y*.8 ~ y*.6 ~ y*.4;

/* scale axes */
xtics(0,15,5,5);
ytics(-.8,.8,.4,4);

/* set globals */
_plwidth = { 0, 2, 4, 6 };
_pnum = 2;
_pltype = 6;
_pcolor = { 10, 2, 12, 4 };
_pgrid = { 1, 2 };
_plctrl = 193;
_psymsiz = 4;

/* setup legend */
_plegctl = { 1 4 };
_plegstr = "f(x)\0f(x)*.8\0f(x)*.6\0f(x)*.4";

/* draw circles in plot coordinates */
csiz = seqa(.1,.3,5);
cclr = { 9, 2, 3, 4, 13 };
linsiz = { 0, 1, 2, 4, 6 };
lintyp = { 6, 6, 6, 6, 3 };
l = ones(5,1);
_pline = l*4 ~ lintyp ~ l*2.5 ~ l*.4 ~ csiz ~ l*0 ~ l*6.28 ~ cclr ~ linsiz;
clear csiz,cclr,l,linsiz;

/* draw arrows in plot coordinates */
a = ones(4,1);
ay = y[193,.]';
ltyp = { 6, 1, 3, 2 };
_parrow = a*7.5~a*.4~a*x[193]~ay~a*5~a*.2~a*21~a*14~a*1~ltyp~a*0;

/* draw an error bar */
_perrbar = { 13.5 13.3 13.7 .65 .55 .75 6 3 0 };

/* extra symbols in plot coordinates */
stype = seqa(1,1,14);
sclr = seqa(9,1,7) | seqa(9,1,7);
sx = seqa(5.5,.67,7) | seqa(5.5,.67,7);
sy = ones(7,1) * -.46 | ones(7,1) * -.54;
s = ones(14,1);
_psym = sx ~ sy ~ stype ~ s*4 ~ sclr ~ s*1;
clear stype,sclr,sx,sy,s;

/* add messages located in plot coordinates */
_pmsgctl = {
1.5 .67 .12 0 1 11 0,
6.0 -.65 .12 0 1 11 0,
8.6  .28 .12 0 1 11 0,
12 -.33 .12 0 1 11 0,
1.1 .11 .12 0 1 11 0,
12.7 .45 .12 0 1 11 0 };
_pmsgstr = "Circles, Arcs\0Extra Symbols\0Arrows\0Auto Legend"\000
"\0Line Width Control\000Error Bars";

fonts("simplex complex");

#ifDOS

/* Use font thickness in title */
title("\213\202Publication Quality Graphics\L\210Example XY graph\201");

#else

title("\202Publication Quality Graphics\lExample XY graph\201");

#endif

xlabel("X Axis Title");
ylabel("Y Axis Title");
xy(x,y);

setwind(2);
graphset;
_pnotify = 0;
_paxes = 0;
_pbox = 15;
fonts("complex");
_pmsgstr = "GAUSS\0Publication Quality Graphics";
_pmsgctl = {
1.5 3.5 2.0 0 2 11 0,
.35 1.0 1.5 0 2 10 0 };
draw;

endwind;

#IFUNIX

call WinSetActive(1);

#ENDIF
