/*
** Monte Carlo investigations of the Hankel singular values
*/


new;
library tsm,optmum,pgraph;

svd_ = {};

i = 1;
do until i>100;

  y = rndn(100,1);
  H = Hankel(y,3);
  svd_1 = svd(H);

  y = recserar(y,1,.95);
  H = Hankel(y-meanc(y),3);
  svd_2 = svd(H);
  svd_ = svd_|(svd_1[1]~svd_2[1]~svd_2[2]);

  i = i+1;
endo;

_kernel = {0,.,.,6,0,.};
d_ = {}; x_ = {};

j = 1;
do until j>3;

  {x,d,F,retcode} = Kernel(svd_[.,j]);

  d_ = d_~d;
  x_ = x_~x;

  j = j + 1;

endo;

graphset;
  title("Monte Carlo investigations"\
        "\LSingular value decomposition of the Hankel matrix"\
        "\Lwhite noise process");
  _pdate = "";  _pnum = 2;
  xy(x_[.,1],d_[.,1]);

  title("Monte Carlo investigations"\
        "\LSingular value decomposition of the Hankel matrix (first sv)"\
        "\LAR(1) process");
  xy(x_[.,2],d_[.,2]);

  title("Monte Carlo investigations"\
        "\LSingular value decomposition of the Hankel matrix (second sv)"\
        "\LAR(1) process");
  xy(x_[.,3],d_[.,3]);


