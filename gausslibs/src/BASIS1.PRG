new;
library tsm,pgraph;

M = 10;
Nobs = 2^M;

call varput(1|1,"Base1");
call varput(1|2|2,"Base2");
call varput(2|2|1,"Base3");
call varput(3|1|3|3|3,"Base4");
call varput(2|2|3|3|3|3,"Base5");
call varput(2|3|4|4|2|3|6|8|8|7|5|4,"Base6");
call varput(2|4|3|4|2|3|6|8|8|7,"Base7");
call varput(M|seqa(M,-1,M),"Base8");


i = 1;
do until i > 8;
  print "The vector ";
  str = "Base" $+ ftos(i,"%lf",1,0);
  B = varget(str);
  print B; print;
  if isbasis(Nobs,B) == 1;
    print "is a basis.";
  else;
    print "is not a basis.";
  endif;

  print; print; print; call pause(2);

  i = i + 1;
endo;


