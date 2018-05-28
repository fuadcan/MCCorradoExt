new;
library tsm,optmum;
declare external _probit_data;

rndseed 123456;

x = rndu(900,10)*5;
beta = rndn(10,1);
y = x*beta+rndn(900,1);
y = y.>-13;

y[4] = miss(0,0);
_tsm_optmum = 0;  /* BHHH algorithm */

output file = tdml4b.out reset;

{theta,stderr,Mcov,Logl} = PROBIT_ml(y,x);

output off;

proc (4) = PROBIT_ml(y,x);
  local k,data,retcode,y_,x_,sv,theta,stderr,Mcov,Logl;

  k = cols(x);

  _probit_data = y~x;
  {data,retcode} = Missing(_probit_data,0);
  y_ = data[.,1];
  x_ = data[.,2:k+1];

  sv = y_/x_;

  __title = "Probit Model";

  _tsm_parnm = 0$+"X"$+ftocv(seqa(1,1,k),2,0);

  /* Analytical Jacobian */

  _ml_Jacobian_proc = &_probit_Jacobian;

  {theta,stderr,Mcov,Logl} = TD_ml(&_probit_ml,sv);

  if _print == 1;
    print;
    call PROBIT_normality(y,x,theta);
  endif;

  retp(theta,stderr,Mcov,Logl);
endp;


proc _probit_ml(theta);
  local k,x,y,w,p,q,logl;

  k = cols(_probit_data)-1;
  y = _probit_data[.,1];
  x = _probit_data[.,2:k+1];

  w = x*theta;
  p = cdfn(w);
  q = 1 - p;
  logl = y.*ln(p)+(1-y).*ln(q);

  retp(logl);
endp;


proc _probit_Jacobian(theta);
  local k,x,y,w,p,q,dp,dq,J;

  k = cols(_probit_data)-1;
  y = _probit_data[.,1];
  x = _probit_data[.,2:k+1];

  w = x*theta;
  p = cdfn(w);
  q = 1 - p;
  dp = pdfn(w).*x;
  dq = -dp;

  J = y.*(dp./p)+(1-y).*(dq./q);

  retp(J);
endp;


/*
                 Bera, A.K., C.M. Jarque and L.F. Lee [1984]
       "Testing the normality assumption in limited dependent variable models"
                 International Economic Review, 25, pagers 563-578
*/


proc (2) = PROBIT_normality(y,x,theta);

  local k,data,T,xb,PHI1,PHI2,w,Xi,M1,M2,M3,Jmatrix,M4;
  local M5,i,M6,M7,M8,M9,test,pvalue;

  k = cols(x);
  data = packr(y~x); y = data[.,1]; x = data[.,2:k+1];
  T = rows(y);
  xb = x*theta;
  PHI1 = pdfn(xb); PHI2 = cdfn(xb);
  w = phi1./(phi2.*(1-phi2));

  Xi = zeros(k+2,T);

                    Xi[1:k,.] = x';
  M1 = xb^2;        Xi[k+1,.] = (M1-1)';
  M2 = xb.*(3+M1);  Xi[k+2,.] = M2';

  M3 = eye(k+2);
  Jmatrix = M3[.,(k+1):(k+2)];

  M4 = zeros(T,2); M4[.,1] = xb; M4[.,2] = 3+M1;
  M5 = PHI1.*w;

  M6 = 0;
  i = 1;
  do until i>T;
    M7 = Xi[.,i]; M6 = M6+M5[i]*(M7*M7');
    i = i+1;
  endo;

  M6 = inv(M6);
  M8 = w.*(y-PHI2).*xb;
  M8 = M4'M8;
  M9 = Jmatrix*M8;


  test = M9'*M6*M9;
  pvalue = cdfchic(test,2);

  if _print == 1;
    print ftos(test,"Normality test : %lf",10,5);
    print ftos(pvalue,"P-value : %lf",10,5);
  endif;

  retp(test,pvalue);
endp;



