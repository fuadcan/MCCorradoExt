new;
cls;

rndseed 3665346;


print "Scatter4.e example demonstrates a basic 3d scatter plot";
print "using random numbers any values less than zero show as triangles";


library pgraph;
graphset;
z = rndn(200,1);
z1 = selif(z,(z .>= 0));
z1 = z1[1:45];
z2 = selif(z, z .< 0);
z2 = z2[1:45];

x1 = rndn(rows(z1),1);
x2 = rndn(rows(z2),1);

y1 = rndn(rows(z1),1);
y2 = rndn(rows(z2),1);

_pcolor = 2|14;
_pmcolor = 3;
_plctrl = -1|-1;
_pstype = 6|3;
_psymsiz = 5|5;
xyz(x1~x2,y1~y2,z1~z2);
