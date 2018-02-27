/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, Cambridge
**
**  Draw the power spectra for stochastic cycles
**  Harvey, pages 61-62
*/


new;
library tsm,optmum,pgraph;

_cycle_prmt = 0;              /* No parametrization */

lambda = seqa(0,pi/1000,1001);

/* rho = 0.7 and lambda_c = pi/4 */

g = _cycle_sgf(0.7|pi/4|1,lambda);

graphset;
  _pdate = "";_pnum = 2; fonts("simplex simgrma");
  title("\201Power spectra for stochastic cycles"\
        "\L\202l\201]c[ = \202p\201/4    \202r\201 = 0.7");
  xlabel("frequency (radian)");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda,g/(2*pi));


/* rho = 0.9 and lambda_c = pi/4 */

g = _cycle_sgf(0.9|pi/4|1,lambda);

  title("\201Power spectra for stochastic cycles"\
        "\L\202l\201]c[ = \202p\201/4    \202r\201 = 0.9");
  xy(lambda,g/(2*pi));

/* rho = 0.99 and lambda_c = pi/4 */

g = _cycle_sgf(0.99|pi/4|1,lambda);

  title("\201Power spectra for stochastic cycles"\
        "\L\202l\201]c[ = \202p\201/4    \202r\201 = 0.99");
  xy(lambda,g/(2*pi));

/*
   First cycle:  rho = 0.8 and lambda_c = pi/10
   Second cycle: rho = 0.9 and lambda_c = pi/2
*/

g = _cycle_sgf(0.8|pi/10|1,lambda) + _cycle_sgf(0.9|pi/2|1,lambda);

  title("\201Power spectra for the sum of two stochastic cycles"\
        "\LFirst cycle: \202l\201]c[ = \202p\201/10    \202r\201 = 0.8"\
        "\LSecond cycle: \202l\201]c[ = \202p\201/2    \202r\201 = 0.9");
  xy(lambda,g/(2*pi));

