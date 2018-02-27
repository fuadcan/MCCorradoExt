new;
library tsm,pgraph;

M = 6;
Nobs = 2^M;

call varput(0,"Base0");
call varput(1|1,"Base1");
call varput(1|2|2,"Base2");
call varput(2|2|1,"Base3");
call varput(3|1|3|3|3,"Base4");
call varput(2|2|3|3|3|3,"Base5");
call varput(2|3|4|4|2|3|6|6|5|4,"Base6");
call varput(2|4|3|4|2|3|6|8|8|7,"Base7");
call varput(M|seqa(M,-1,M),"Base8");
call varput(M*ones(Nobs,1),"Base9");


i = 0;
do until i > 9;
  str = "Base" $+ ftos(i,"%lf",1,0);
  B = varget(str);
  if isbasis(Nobs,B) == 0 and B/=0;
    i = i + 1;
    continue;
  endif;
    
  graphset;
    _pdate = "";
    title("The TIME/FREQUENCY representation of the basis "$+str$+" is ");
    call BasisPlot(B,M);
    draw;

  i = i + 1;
endo;



