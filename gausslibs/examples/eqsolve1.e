eqsolveset;

proc fct(x);
    local f1,f2,f3;
    f1 = 3*x[1]^3 + 2*x[2]^2 + 5*x[3] - 10;
    f2 = -x[1]^3 - 3*x[2]^2 + x[3] + 5;
    f3 = 3*x[1]^3 + 2*x[2]^2 -4*x[3];
    retp(f1|f2|f3);
endp;

proc fjc(x);
    local fjc1,fjc2, fjc3;
    fjc1 = 9*x[1]^2~4*x[2]~5;
    fjc2 = -3*x[1]^2~-6*x[2]~1;
    fjc3 = 9*x[1]^2~4*x[2]~-4;
    retp(fjc1|fjc2|fjc3);
endp;

start = { -1, 12, -1 };

_eqs_JacobianProc = &fjc;

output file = eqsolve1_out.txt reset;

{ x,tcode }  = eqSolve(&fct,start);

output off;

print "The results of this run will be found in your working";
print "directory in the file eqsolve1_out.txt";

