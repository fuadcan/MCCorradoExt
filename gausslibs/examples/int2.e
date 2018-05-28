cls;
print " This example int2.e Calculates the Volume of a sphere ";
print;

radius = 3;

proc f(x,y,z);
   retp(1);
endp;

proc h1(x,y);
    retp(sqrt((radius^2) - (x^2) - (y^2)));
endp;

proc h2(x,y);
    retp(-h1(x,y));
endp;

proc g1(x);
    retp(sqrt((radius^2) - x^2));
endp;

proc g2(x);
    retp(-g1(x));
endp;

xl = radius|-radius;
_intord = 40;
_intrec = 0;
vol = intgrat3(&f,xl,&g1|&g2,&h1|&h2);
format /ld 10,4;
print "VOLUME OF A SPHERE";
print "       Radius:   " radius;
print;
print "Answer calculated by INTGRAT3: " vol;
print "True answer using 4/3*PI*R^3:  " 4/3*pi*radius^3;

