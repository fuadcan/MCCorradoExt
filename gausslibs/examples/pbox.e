library pgraph;
graphset;

rndseed 436566;

#IFUNIX

let v = 100 100 640 480 40 80 1 6 15 0 0 2 2;
wbx = WinOpenPQG( v, "Box Plots", "Boxes" );
call WinSetActive( wbx );

#ENDIF

/* Create random data */
nboxes = 10;
npoints = 20;
xtics(0,nboxes+1,1,1);
x=seqa(1,1,nboxes)';
y=rndn(npoints,nboxes);

/* setup globals common to all plots */
_plctrl = 1;            @ Symbols outside whisker limits will be visible  @
_pstype = 3;            @ All symbols will be triangles                   @
_pcolor = 3;            @ All symbol colors will be cyan                  @
_pnum = 0;              @ No numbers on axes                              @
_ptitlht = .2;          @ Title character height                          @
margin(.1,.1,.2,.1);    @ Margin around each window                       @

@ init windows @
begwind;
window(2,2,0);
makewind(9,6.855,0,0,1);

@ plot 1: use max-min for whisker limits, norm width, lt green boxes... @
title("Max-min whiskers\lWhisker limits are actual limits of data");
_pboxctl = { .4 10 1 };
box(x,y);

@ plot 2: statist. standard whisker limits, small width, yellow boxes... @
nextwind;
title("Statistical standard\lOnly symbols outside whisker range plotted");
_plctrl = 1;
_pboxctl = { .25 14 2 };
box(x,y);


@ plot 3: Plot symbol also, norm width, lt red boxes... @
nextwind;
title("With symbols\lWhisker limits are at 10th, 90th percentiles");
_plctrl = 2;
_pboxctl = { .5 12 3 10 90 };
box(x,y);

@ plot 4: max-min for whisker limits, minimized ink width, varied color boxes @
nextwind;
title("Minimized ink method\lZero width boxes");
_pboxctl = { 0 0 1 };
_pcolor = 0;
_plctrl = 0;
box(x,y);
setwind(5);

@ add a Main title to page @
fonts("complex");
_ptitle = "";
_paxes = 0;
_pmsgctl = { 3.24 6.6 .17 0 2 10 0 };
_pmsgstr = "\201Box/Whisker Plots";
draw;

@ window manip done, show tekfile... @
endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF
