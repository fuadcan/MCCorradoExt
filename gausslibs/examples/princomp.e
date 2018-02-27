 cls;
print "This example demonstrates principle component analysis";
print;
print;

    rndseed 5645658;

    cr = { 1  .8  .4  .3,     /* correlation matrix of */
          .8   1  .6  .2,     /* independent variables */
          .4  .6   1  .1,
          .3  .2  .3   1 };

    x = rndn(100,4)*chol(cr);

    { p,v,a } = princomp(x,4);

    print "fraction of variance explained by each component:";
    print;
    print v';

    print;
    print "loadings";
    print;
    print a;
