new;
library tsm,optmum;
TSMset;

rndseed 123456;

sv = 0.5|0.25|0.7|1;
z = Process;

_tsm_Mcov = 0;
_tsm_parnm = "phi1"|"sigma_u"|"theta1"|"sigma_v";
_fourier = 0;
_sgf_Jacobian_Proc = &sgf_Jacobian;


/* Unrestricted model */

{coeff1,stderr1,Mcov1,Logl1} = FD_ml(z,&sgf,sv);

/* Restricted model */

RR = design(1|2|0|3); r = 0|0|0.7|0;
sv = inv(RR'RR)*RR'sv;

{coeff2,stderr2,Mcov2,Logl2} = FD_cml(z,&sgf,sv,RR,r);

output file = fdml2b.out reset;

LR = 2*(logl1-logl2);      /* Likelihood ratio                              */
pvalue = cdfchic(LR,1);    /* Approximating the noncentral chi-squared CDF */

print ftos(LR,"Likelihood ratio statistic: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

{J_matrix,G_matrix,H_matrix,I_matrix} = FDml_derivatives(z,&sgf,coeff2);

LM = G_matrix'*inv(I_matrix)*G_matrix;
pvalue = cdfchic(LM,1);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;

proc sgf(coeff,lambda);
  local phi1,theta1,sigma_u,sigma_v;
  local w,g;
  phi1 = coeff[1];
  sigma_u = coeff[2];
  theta1 = coeff[3];
  sigma_v = coeff[4];
  w = cos(lambda);
  g = (sigma_u^2)./(1-2*phi1*w+phi1^2) +
      (sigma_v^2).*(1-2*theta1*w+theta1^2);
  retp(g);
endp;

proc sgf_Jacobian(coeff,lambda);
  local phi1,theta1,sigma_u,sigma_v;
  local w,g;
  local w1,w2,J1,J2,J3,J4,J;

  phi1 = coeff[1];
  sigma_u = coeff[2];
  theta1 = coeff[3];
  sigma_v = coeff[4];
  w = cos(lambda);

  w1 = 1-2*phi1*w+phi1^2;
  w2 = 1-2*theta1*w+theta1^2;

  J1 = 2*(w-phi1).*(sigma_u^2)./(w1^2);
  J2 = 2*sigma_u./w1;
  J3 = 2*(theta1-w).*(sigma_v^2);
  J4 = 2*sigma_v.*w2;
  J = J1~J2~J3~J4;

  retp(J);
endp;

proc process;
  local Nobs,u,x,v,y,z;
  Nobs = 1000;
  u = rndn(Nobs,1)*0.25;
  x = recserar(u,0,0.50);
  v = rndn(Nobs,1)*1;
  y = v - 0.7*(0|trimr(v,0,1));
  z = x + y;
  retp(z);
endp;

