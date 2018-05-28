cls;
print;
print "This example calculates the surface area of a spherical cap of";
print "height h (measured from top of sphere) of a sphere of radius r.";
print "The example uses the function INTGRAT2.";
print;
print "Press any key to continue"; call keyw;
print;


h = 0.02;
r = 5;

proc f(x,y);
   local px2,py2;
   px2 = x^2./(r^2-x^2-y^2);
   py2 = y^2./(r^2-x^2-y^2);
   retp(sqrt(1 + px2 + py2));
endp;

proc g1(x);
    retp(sqrt((2*r*h-h^2) - x^2));
endp;

proc g2(x);
    retp(-g1(x));
endp;

_intord = 64;
_intrec = 0;
xl =  sqrt(2*r*h - h^2)|-sqrt(2*r*h - h^2);
y = intgrat2(&f,xl,&g1|&g2);
format /ld 3,3;
print "SURFACE AREA OF SPHERICAL SEGMENT";
print "       Height of Segment: " h;
print "       Radius of Sphere:  " r;
print;
format /ld 8,4;
print "Answer calculated using INTGRAT2:  " y;
print "True answer using 2(pi)Rh:         " 2*pi*r*h;

