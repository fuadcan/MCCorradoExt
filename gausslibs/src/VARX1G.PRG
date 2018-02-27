/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  VAR order selection
**  See LUTKEPOHL, chapter 4
**
*/

new;
library tsm,optmum,pgraph;

cls;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
data = dinv~dinc~dcons;

_print = 1;

output file = varx1g.out reset;

{CR,p} = criteria(data,1,5);

output off;

/*
**  A procedure to select the order of a VARX model
*/

proc (2) = criteria(Y,X,p_max);
  local old,K,T,m,CR,L,D,SIGMA;
  local omat,fmt,mask,p;

  old = _print;
  K = cols(Y);
  CR = zeros(p_max+1,7);

  m = 0;
  do until m>p_max;

    _print = 0;
    call varx_ML(Y,X,m);
    T = rows(packr(Y))-m;

    SIGMA = _varx_SIGMA;
    D = det(SIGMA);
    L = ln(D);

    CR[1+m,1] = L+(K^2)*m*ln(T)/T;                  /* BIC */
    CR[1+m,2] = L+3*m*(K^2)/T;                      /* AIC alpha = 3 */
    CR[1+m,3] = L+(1+m*(K^2)/T)/(1-(m*(K^2)-2)/T);  /* AICc */
    CR[1+m,4] = D*(1+2*(m*(K^2)+1)/T)^K;            /* SIC  */
    CR[1+m,5] = D*(((T+K*m+1)/(T-K*m-1))^K);        /* MFPE */
    CR[1+m,6] = L+2*(K^2)*m/T;                      /* AIC  */
    CR[1+m,7] = L+2*(K^2)*m*ln(ln(T))/T;            /* HQ   */

    m = m+1;
  endo;

  P = minindc(CR)-1;
  _print = old;

  if _print == 1;

    m = seqa(0,1,p_max+1);

    print chrs(45*ones(79,1));
    print "Lags      BIC       AICa      AICc        SIC         FPE       AIC"\
           "       HQ";
    print chrs(45*ones(79,1));

    omat = "Opt."~p';
    mask = 0~1~1~1~1~1~1~1;
    let fmt[8,3] =
       "-*.*s"  4 4
       "*.*lf" 8  0
       "*.*lf" 10 0
       "*.*lf" 10 0
       "*.*lf" 12 0
       "*.*lf" 12 0
       "*.*lf" 10 0
       "*.*lf" 10 0;
    call printfm(omat,mask,fmt); print; print chrs(45*ones(79,1));

    omat = m~CR;
    mask = 1;
    let fmt[8,3] =
       "*.*lf"  4 0
       "*.*lf" 10 3
       "*.*lf" 10 3
       "*.*lf" 10 3
       "*.*lf" 12 3
       "*.*lf" 12 3
       "*.*lf" 10 3
       "*.*lf" 10 3;
    call printfm(omat,mask,fmt);
  endif;
  retp(CR,p);
endp;





