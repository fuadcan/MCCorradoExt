new;
library tsm,optmum;
format /m1 /rd 10,2;

cls;

x = rndn(3,3);
y = x'x;         /* Symmetric matrix */

print "                y"; print y; print; call pause(3);
v1 = vech(y);
v2 = vech_(y);
print "   vech(y)    vech_(y) "; print v1~v2; print; call pause(3);

y1 = xpnd(v1);
y2 = xpnd_(v2);
print "   xpnd(vech(y))    "; print y1; print;
print "   xpnd_(vech_(y))    "; print y2; print;

