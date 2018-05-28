library pgraph;
graphset;

#IFUNIX

let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
wxy1 = WinOpenPQG( v, "Trochiods", "Trochs" );
call WinSetActive( wxy1 );

#ENDIF

/* create 2 windows side by side */
begwind;
window(1,2,0);

/* setup globals for both graphs */
_pltype = 6;
_paxes = 0;
_pdate = "";
_ptitlht = .25;
_plegctl = { 1 7 };
_plwidth = { 4, 8, 6 };
t = seqa(0,.0314,201);
setwind(1);
_ptitle = "HYPOTROCHOIDS"\0
"\lx = d*((a-b)*cos(t) + c*cos((a-b)t/b))"\0
"\ly = d*((a-b)*sin(t) - c*sin((a-b))t/b))";
_plegstr = "a=4, b=1, c=2, d=.15\0a=4, b=1, c=3, d=.15\0a=4, b=1, c=1, d=.25";
_pcolor = { 11, 9, 13 };
{ x, y} = hypo(4,1,2,.15);
{x1,y1} = hypo(4,1,3,.15);
{x2,y2} = hypo(4,1,1,.25);
x = x ~ x1 ~ x2;
y = y ~ y1 ~ y2;
xy(x,y);
clear x,x1,x2,y,y1,y2;

nextwind;
_ptitle = "EPITROCHOIDS"\0
"\lx = d*((a+b)*cos(t) - c*cos((a+b)t/b))"\0
"\ly = d*((a+b)*sin(t) - c*sin((a+b)t/b))";
_plegstr = "a=6, b=1, c=1, d=.2\0a=6, b=1, c=4, d=.25\0a=6, b=2, c=1, d=.25";
_pcolor = { 14, 2, 6 };
{ x, y} = epit(6,1,1,.2);
{x1,y1} = epit(6,1,4,.25);
{x2,y2} = epit(6,2,1,.25);
x = x ~ x1 ~ x2;
y = y ~ y1 ~ y2;
xy(x,y);

endwind;

#IFUNIX

call WinSetActive( 1 );

#ENDIF

proc (2) = hypo(a,b,c,d);
    local x0,y0;
    x0= d*( (a-b)*cos(t) + c*cos( (a-b)*t/b ) );
    y0= d*( (a-b)*sin(t) - c*sin( (a-b)*t/b ) );
    retp(x0,y0);
endp;

proc (2) = epit(a,b,c,d);
    local x0,y0;
    x0 = d*( (a+b)*cos(t) - c*cos( (a+b)*t/b ) );
    y0 = d*( (a+b)*sin(t) - c*sin( (a+b)*t/b ) );
    retp(x0,y0);
endp;

