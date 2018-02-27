new;
library tsm;

output file = qmf1.out reset;

cls;

{H,G,Htilde,Gtilde} = Haar;

print "Haar:   "; print;
call print_filter(H,G);

i = 4;
do until i > 20;

  {H,G,Htilde,Gtilde} = Daubechies(i);

  str = ftos(i,"Daubechies #%lf:   ",1,0);
  print str; print;
  call print_filter(H,G);

  i = i + 2;

endo;


i = 2;
do until i > 6;

  {H,G,Htilde,Gtilde} = Coiflet(i);

  str = ftos(i,"Coiflet #%lf:   ",1,0);
  print str; print;
  call print_filter(H,G);

  i = i + 2;

endo;


output off;


proc (0) = print_filter(H,G);
  local nH,s,h_,g_,str,omat,mask,fmt;

  nH = rows(H);
  s = seqa(0,1,nH);

  h_ = 0 $+ "h" $+ ftocv(s,1,0);
  g_ = 0 $+ "g" $+ ftocv(s,1,0);

  str = zeros(nH,1) $+ "";

  omat = h_~H~str~g_~G;
  mask = 0~1~0~0~1;

  print "  Low-pass filter                        High-pass filter";
  print "------------------------------------------------------------";
  let fmt[5,3] =  " -*.*s"  7 6
                  "*.*lf"  12 5
                  "-*.*s"  20 6
                  " -*.*s"  7 6
                  "*.*lf"  12 5;
  call printfm(omat,mask,fmt); print;
  print ftos(sumc(H),"sum of hi:    %lf",10,4);
  print ftos(H'H,"sum of hi^2:  %lf",10,4);
  print ftos(sumc(G),"sum of gi:    %lf",10,4);
  print ftos(G'G,"sum of gi^2:  %lf",10,4);
  print ftos(H'G,"sum of hi*gi: %lf",10,4);
  print;
  retp;
endp;



