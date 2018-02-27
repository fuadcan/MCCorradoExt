/*
** BROCKWELL and DAVIS [1991], Time Series: Theory and Methods,
** second edition, Springer-Verlag, New-York
**
** Example 11.7.1, pages 451-453
*/


new;
library tsm,optmum,pgraph;
TSMset;

let lead =
     0.060       0.25      -0.57       0.58      -0.20       0.23     -0.040
     -0.19      0.030       0.42      0.040       0.24       0.34      -0.46
     -0.18     -0.080       0.29       0.56      -0.37       0.20       0.54
     -0.31      0.030       0.52      -0.70       0.35      -0.63       0.44
     -0.38     -0.010       0.22       0.10      -0.50      0.010       0.30
     -0.76       0.52       0.15      0.060      -0.10       0.21     -0.010
      0.70      -0.22      -0.76      0.060      0.020      -0.17     -0.080
     0.010       0.11      -0.39      0.010       0.50     -0.020      -0.37
     -0.13      0.050       0.54      -0.46       0.25      -0.52       0.44
     0.020      -0.47       0.11      0.060       0.25      -0.35        0.0
    -0.060       0.21     -0.090       0.36      0.090     -0.040      -0.20
      0.44      -0.23       0.40     -0.010       0.17      0.080       0.58
     -0.27       0.79      -0.21      0.020       0.30       0.28      -0.27
    -0.010      0.030       0.16      -0.28       0.15       0.26      -0.36
      0.32      -0.11       0.22      -0.65        0.0       0.47       0.16
     -0.19       0.48      -0.26       0.21        0.0      -0.20       0.35
      0.38      -0.48       0.20      -0.32       0.43      -0.50       0.12
     -0.17       0.15      -0.36       0.35     -0.030      -0.18       0.16
     0.070       0.21      -0.50       0.23      -0.13       0.14      -0.15
      0.19      -0.24       0.26      -0.22       0.17       0.37     -0.060
      0.29      -0.34      -0.12      -0.16       0.25      0.080     -0.070
      0.26      -0.37;

let sales =
     -0.60      -0.10      -0.50       0.10        1.2       -1.6        1.4
      0.30       0.90       0.40      -0.10        0.0        2.0        1.4
       2.2        3.4        0.0      -0.70       -1.0       0.70        3.7
      0.50        1.4        3.6        1.1       0.70        3.3       -1.0
       1.0       -2.1       0.60       -1.5       -1.4       0.70       0.50
      -1.7       -1.1      -0.10       -2.7       0.30       0.60       0.80
       0.0        1.0        1.0        4.2        2.0       -2.7       -1.5
     -0.70       -1.3       -1.7       -1.1      -0.10       -1.7       -1.8
       1.6       0.70       -1.0       -1.5      -0.70        1.7      -0.20
      0.40       -1.8       0.80       0.70       -2.0      -0.30      -0.60
       1.3       -1.4      -0.30      -0.90        0.0        0.0        1.8
       1.3       0.90      -0.30        2.3       0.50        2.2        1.3
       1.9        1.5        4.5        1.7        4.8        2.5        1.4
       3.5        3.2        1.5       0.70       0.30        1.4      -0.10
      0.20        1.6      -0.40       0.90       0.60        1.0       -2.5
      -1.4        1.2        1.6       0.30        2.3       0.70        1.3
       1.2      -0.20        1.4        3.0      -0.40        1.3      -0.90
       1.2      -0.80       -1.0      -0.80      -0.10       -1.5       0.30
      0.20      -0.50      -0.10       0.30        1.3       -1.1      -0.10
     -0.50       0.30      -0.70       0.70      -0.50       0.60      -0.30
      0.20        2.1        1.5        1.8       0.40      -0.50       -1.0
      0.40       0.50;


_fourier = 1;
{lambda,I1} = PDGM(lead);
{lambda,I2} = PDGM(sales);
{lambda,I12} = CPDGM(lead,sales);

/* User lag window generator */

/* The weight function is W(k) = 1/13 */
/* The bandwidth is M = 6 */

_smoothing = {6,0,0,0.25};
_window = ones(13,1)/13;

I1 = smoothing(I1);
I2 = smoothing(I2);
I12 = smoothing(I12);


c = lambda/(2*pi);

{w,alpha,phi} = Cspectrum(I1,I2,I12);

phi = phi[1:75];
alpha = alpha[1:75];
c = c[1:75];
lambda = lambda[1:75];
I1 = I1[1:75];
I2 = I2[1:75];
I12 = I12[1:75];

graphset;
  _pdate = "";
  xlabel("c");
  xtics(0,0.5,0.1,10);

  title("spectral density estimate"\
        "\LFigure 11.5 of Brockwell and Davis @[1991@]");
  ytics(0,0.05,0.005,0);
  xy(c,I1/(2*pi));

  title("spectral density estimate"\
        "\LFigure 11.6 of Brockwell and Davis @[1991@]");
  ytics(0,1.6,0.1,0);
  xy(c,I2/(2*pi));

  title("cospectrum");
  ytics(-0.5,0.4,0.1,0);
  xy(c,real(I12));

  title("quadrature spectrum");
  ytics(-0.4,0.4,0.1,0);
  xy(c,-imag(I12));

  title("absolute coherency or cross-amplitude spectra"\
        "\LFigure 11.7 of Brockwell and Davis @[1991@]");
  ytics(0,1,0.1,0);
  xy(c,alpha);

/* The slope at the low frequencies is 4.1 */
L1 = 4.1*lambda;
/* The slpoe at the other frequencies is 2.7 */
L2 = 2.7*lambda - 5.25;

  title("phase spectra"\
        "\LFigure 11.8 of Brockwell and Davis @[1991@]");
  ytics(-4,4,1,0);
  _pltype = 6|1|1;
  xy(c,phi~L1~L2);

