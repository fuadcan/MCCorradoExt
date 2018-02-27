// This program creates 9 separate VWR windows in specific
// locations on the screen.

// To close the graphs go to the main GAUSS menu "Window"
// and select "Close all Graphics".

library pgraph;
graphset;

pqgwin many;

x = seqa(0,.1,600);
y = sin(x);

// set length of windows in pixels
xlen = 245;
ylen = 220;

// cycle through rows of vwr windows
for i (0,2,1);

    ypos = i * (ylen+10) + 20;

    // cycle through columns of vwr windows
    for j (0,2,1);

        // set window params, x,y,xlen,ylen
        xpos = j * (xlen+10) + 20;

        // build cmd arg string
        cmdarg = "-wp=" $+ vwrwndpos( xpos, ypos, xlen, ylen );
        graphprt(cmdarg);
        y = y ~ y*1.4;
        xy(x,y);
    endfor;
endfor;

// this proc places the vwr window coordinates into a cmd string
proc (1) = vwrwndpos( xpos, ypos, xlen, ylen );
    local str1, str2, str3, str4;

    str1 = ftos(xpos,"%*.*lf",1,0) $+ ",";
    str2 = ftos(ypos,"%*.*lf",1,0) $+ ",";
    str3 = ftos(xlen,"%*.*lf",1,0) $+ ",";
    str4 = ftos(ylen,"%*.*lf",1,0);
	retp( str1 $+ str2 $+ str3 $+ str4 );
endp;

