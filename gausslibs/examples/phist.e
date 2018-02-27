library pgraph;
graphset;

rndseed 454356;

#IFUNIX

let v = 100 100 640 480 40 80 1 6 15 0 0 2 2;
phw = WinOpenPQG( v, "Histogram", "Hist" );
call WinSetActive( phw );

#ENDIF

_pcolor = 12;
_pmcolor = 3;
_pcsel = 1;
x = round(rndn(4000,1)*100);
call hist(x,50);

#IFUNIX

call WinSetActive( 1 );

#ENDIF

