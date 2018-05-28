cls;
print "Gauss example sqpsolve.e";
print;

sqpSolveSet;

proc fct(x);
   retp( (x[1] + 3*x[2] + x[3])^2 + 4*(x[1] - x[2])^2 );
endp;

proc ineqp(x);
    retp(6*x[2] + 4*x[3] - x[1]^3 - 3);
endp;

proc eqp(x);
    retp(1-sumc(x));
endp;

_sqp_Bounds = { 0 1e256 };

start = { .1, .7, .2 };

_sqp_IneqProc = &ineqp;
_sqp_EqProc = &eqp;

{ x,f,lagr,ret } = sqpSolve( &fct,start );

print;
print "published solution";
print " 0      0       1";
output off;

