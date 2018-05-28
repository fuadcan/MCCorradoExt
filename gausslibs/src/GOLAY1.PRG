/*
** Press, W.H., S.A. Teukolsky, W.T. Vetterling and B.P. Flannery [1992],
** Numerical Recipes in Fortran, second edition, Cambridge University
** Press, Cambridge
**
** Compute the smoothed coefficients for differents values of M, nL and nR.
**
** See Press and al., page 646
**
*/

new;
library tsm,pgraph;

cls;

format /rdn 7,3;

m = miss(zeros(3,1),0);

call Savitzky_Golay(0,2,2,2);
w = m|_window;
print w';

call Savitzky_Golay(0,2,3,1);
w = m[1:2]|_window;
print w';

call Savitzky_Golay(0,2,4,0);
w = m[1]|_window;
print w';

call Savitzky_Golay(0,2,5,5);
w = _window;
print w';

call Savitzky_Golay(0,4,4,4);
w = m[1]|_window;
print w';

call Savitzky_Golay(0,4,5,5);
w = _window;
print w';



