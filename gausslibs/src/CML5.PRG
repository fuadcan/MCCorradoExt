new;
library tsm,optmum,CML;

TSMset; CMLset;

/*  Generate an ARMA(1,1) process */

rndseed 123;

et = rndn(101,1)*1.5;

yt = recserar(et[2:101]+0.6*et[1:100],0|0,0.7|0.2);


/*  Build a procedure for the likelihood function */

proc ml_(coeff,y);
  local beta,sigma,Pchol,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  beta = coeff[1:3];
  SIGMA = coeff[4]^2;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(y,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

proc ml(coeff);
  retp(ml_(coeff,yt));
endp;


_cml_C = {  1  1 0 0,
           -1 -1 0 0};
_cml_D = { 0.8, -1};

{beta,f,g,cov,ret} = CMLPrt(CML(yt,0,&ml_,0.7|0.2|-0.6|1.5));

{J_matrix,G_matrix,H_matrix,I_matrix} = TDml_derivatives(&ml,beta);

/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/

I1 = -H_matrix;          /* First approximation of the Information matrix  */
I2 = J_matrix'J_matrix;  /* Second approximation of the Information matrix */

LM1 = G_matrix'*inv(I1)*G_matrix;
pvalue = cdfchic(LM1,15);

print;
print ftos(LM1,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


LM2 = G_matrix'*inv(I2)*G_matrix;
pvalue = cdfchic(LM2,15);

print;
print ftos(LM2,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);


