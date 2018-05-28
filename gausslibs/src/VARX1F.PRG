/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Impulse Response Analysis
**  See LUTKEPOHL, section 3.7.3
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

_print = 0;
{theta,stderr,Mcov,LOGL} = varx_LS(data,1,2);  /* Constant and 2 lags */

beta = theta[1:18];  /* Estimates of the AR part */
SIGMA = _varx_SIGMA;

Nr = 8;
t = seqa(0,1,Nr+1);

X1 = zeros(1+Nr,2);  /* Responses of Consumption to an impulse in Income */
X2 = zeros(1+Nr,2);  /* Responses of Investment to an impulse in Consumption */
X3 = zeros(1+Nr,2);  /* Responses of Consumption
                        to an orthogonal impulse in Income */

/* Responses to Forecast Errors */

call arma_impulse(beta,2,0,Nr);

i = 0;
do until i>Nr;
  x = varget("IMPULSE"$+ftos(i,"%lf",1,0));
  y = varget("_IMPULSE"$+ftos(i,"%lf",1,0));
  X1[1+i,.] = x[3,2]~y[3,2];
  X2[1+i,.] = x[1,3]~y[1,3];
  i = i+1;
endo;

let m = {.,.};

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.20; _ptitlht = 0.25;
  begwind;
  window(2,1,0);
  setwind(1);
    title("Estimated responses of consumption to a forecast error impulse in "\
          "income"\
          "\LSee LUTKEPOHL, page 107");
    xy(t,X1[.,1]);
  nextwind;
    title("Accumulated responses of consumption to a forecast error impulse in"\
          " income"\
          "\LSee LUTKEPOHL, page 109");
    xy(t,X1[.,2]);
  endwind;


begwind;
window(2,1,0);
setwind(1);
  title("Estimated responses of investment to a forecast error impulse in "\
        "consumption"\
        "\LSee LUTKEPOHL, page 107");
  xy(t,X2[.,1]);
nextwind;
  title("Accumulated responses of investment to a forecast error impulse in "\
        "consumption"\
        "\LSee LUTKEPOHL, page 109");
  xy(t,X2[.,2]);
endwind;


/* Responses to Orthogonal Impulses */

call arma_orthogonal(beta,2,0,SIGMA,Nr);

i = 0;
do until i>Nr;
  x = varget("IMPULSE"$+ftos(i,"%lf",1,0));
  y = varget("_IMPULSE"$+ftos(i,"%lf",1,0));
  X3[1+i,.] = x[3,2]~y[3,2];
  i = i+1;
endo;

begwind;
window(2,1,0);
setwind(1);
  title("Estimated responses of consumption to an orthogonalized impulse in "\
        "income"\
        "\LSee LUTKEPOHL, page 112");
  xy(t,X3[.,1]);
nextwind;
  title("Accumulated responses of consumption to a orthogonalized impulse in "\
        "income"\
        "\LSee LUTKEPOHL, page 112");
  xy(t,X3[.,2]);
endwind;



