/*
**  Using vread to obtain the gradient and the Hessian
**  of the likelihood function
*/

new;
library optmum,tsm;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

_print = 0;
_tsm_optmum = 1;
_tsm_gtol = 0.001;

w = 1|0|0|2|3|0|0|4;
RR = design(w); r = zeros(8,1);

load sv = arma2; sv = sv[1 4 5 8];

{beta,stderr,Mcov,LogL} = arma_CML(y,1,1,sv,RR,r);
sigma = _arma_sigma;


outwidth 256;

output file = arma1g.out reset;

print "SIGMA: "; print sigma; print;

Jacobian = vread(_ml_derivatives,"J_matrix");
print "Jacobian: "; print Jacobian; print;

Gradient = vread(_ml_derivatives,"G_matrix");
print "Gradient: "; print Gradient; print;

Hessian = vread(_ml_derivatives,"H_matrix");
print "Hessian: "; print Hessian; print;

Information = vread(_ml_derivatives,"I_matrix");
print "Information: "; print Information; print;

output off;




