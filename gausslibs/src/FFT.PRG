new;
library tsm;
TSMset;

load x[] = lynx.asc;

x = log(x);
Nobs = rows(x);

d1 = dfft(x);
d2 = fft(x);

/* The first Fourier frequency corresponds to the mean of the data */

print ftos(real(d1[1]),"dfft(0) = %lf",4,3);
print ftos(real(d2[1]),"fft(0) =  %lf",4,3);
print;

_fourier = 1;
d1 = fourier(x)/Nobs;
_fourier = 0;
d2 = fourier(x)/Nobs;

print ftos(real(d1[1]),"dfft(0) = %lf",4,3);
print ftos(real(d2[1]),"fft(0) =  %lf",4,3);

