/*
** Press, W.H., S.A. Teukolsky, W.T. Vetterling and B.P. Flannery [1992],
** Numerical Recipes in Fortran, second edition, Cambridge University
** Press, Cambridge
**
** See Press and al., section 14.8
**
*/

new;
library tsm,pgraph;

rndseed 123;

t = seqa(1,1,300);
xx = cos(1+t/300).*(0.5+rndu(300,1));
xx = recserar(xx,0,0.8);

xx[8 25 56 98] = miss(zeros(4,1),0);

x = xx + rndn(300,1)*0.1;

xs = Savitzky_Golay(x,2,5,5);


graphset;
  begwind;
  makewind(9,6.855,0,0,0);
  makewind(9/2-1,6.855/2-1,1,1,0);
  setwind(1);
    _pdate = ""; _pnum = 2;
    title("Filtered series");
    xy(t,xs);
  setwind(2);
    _pnumht = 0.30; _ptitlht = 0.3;
    title("Function plus noise");
    xy(t,x);
  endwind;


