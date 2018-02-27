new;
library tsm,optmum,pgraph;
TSMset;

rndseed 123456;

_fourier = 1;

Z = eye(2); d = 0|0;
let H[2,2] = 0.2 0 0 0.1;
let T[2,2] = 0.5 0.3 0 -0.5;
c = 0|0; R = 1|1; Q = 0.25;

call SSM_build(Z,d,H,T,c,R,Q,0);
{y,a} = RND_SSM(0|0,100);

{lambda,I} = PDGM2(y);
g = sgf_SSM(lambda);

/* first component */

CV1a = real(inverse_fourier(I[.,1]));
CV1b = real(inverse_fourier(g[.,1]));

/* second component */

CV2a = real(inverse_fourier(I[.,4]));
CV2b = real(inverse_fourier(g[.,4]));

/* first component/second component */

CV3a = real(inverse_fourier(I[.,3]));
CV3b = real(inverse_fourier(g[.,3]));

/* second component/first component */

CV4a = real(inverse_fourier(I[.,2]));
CV4b = real(inverse_fourier(g[.,2]));

graphset;
    _pdate = "";
    title("Covariogram of the first component of the multivariate process");
    _pltype = 6|1;
    _plegstr = "estimated\0theoretical";
    _plegctl = {2 6 5 3};
    xlabel("Lag");
    xy(seqa(0,1,11),CV1a[1:11]~CV1b[1:11]);

graphset;
    _pdate = "";
    title("Covariogram of the second component of the multivariate process");
    _pltype = 6|1;
    _plegstr = "estimated\0theoretical";
    _plegctl = {2 6 5 1.25};
    xlabel("Lag");
    xy(seqa(0,1,11),CV2a[1:11]~CV2b[1:11]);

graphset;
    _pdate = "";
    title("Cross-covariogram of the multivariate process"\
          "\Lcov(y]1[(t),y]2[(t-lag))");
    _pltype = 6|1;
    _plegstr = "estimated\0theoretical";
    _plegctl = {2 6 5 3.25};
    xlabel("Lag");
    xy(seqa(0,1,11),CV3a[1:11]~CV3b[1:11]);


graphset;
    _pdate = "";
    title("Cross-covariogram of the multivariate process"\
          "\Lcov(y]2[(t),y]1[(t-lag))");
    _pltype = 6|1;
    _plegstr = "estimated\0theoretical";
    _plegctl = {2 6 5 4};
    xlabel("Lag");
    xy(seqa(0,1,11),CV4a[1:11]~CV4b[1:11]);


