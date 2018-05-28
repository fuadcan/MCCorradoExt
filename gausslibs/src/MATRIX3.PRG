new;
library tsm,optmum;

M = 1;

do while M <= 10;

  p = M*(M+1)/2;

  L = elimination_(M);
  D = duplication_(M);
  K = commutation_(M,M);
  K1 = commutation_(M,1);
  K2 = commutation_(1,M);

  e1 = ( L*D == eye(p) ) ;
  e2 = ( K*D == D ) ;
  e3 = ( ( K1 == K2 ) == (K1 == eye(M) ) ) ;
  e4 = ( sumc(diag(K)) == M ) ;
  e5 = ( sumc(diag(D'D)) == M^2 ) ;
  e6 = ( L*L' == eye(p) ) ;
  e7 = ( sumc(diag(invpd(D'D))) == M*(M+3)/4 );

  e = e1|e2|e3|e4|e5|e6|e7;

  if e == ones(7,1);
    print "The propositions are verified for M = "$+ftos(M,"%lf",2,0);
  else;
    print "The propositions are not verified for M = "$+ftos(M,"%lf",2,0);
  endif;

  M = M+1;

endo;



