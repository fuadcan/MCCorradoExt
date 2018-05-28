// This example simulates animation by updating a single
// viewer window with sequential graphs.";

library pgraph;
graphset;

// setup color table
ctable = { 15 1, 14 2, 13 3, 12 4, 10 5, 9 6, 8 7, 7 8 };

// position viewer window
graphprt("-wp=200,200,440,380");

/* create data */
i=0;
k=50;
z=20;
q=0.5;

for i (0,300,1);
    clrindex = floor(i/20)%8 + 1;
    _pcolor = ctable[clrindex,.];
  
    n = seqa(0,1,101);
    theta = n.*pi/z;
    radius = sin(q .* theta);
    theta1 = seqa(i/3,pi/k,101);
    radius1 = sin(theta);
    radius2 = cos(theta);
    _pdate="";
    _pgrid = 2;
    _plwidth = 9;
    _pltype = 0;
    _ptitlht = .2;
    _paxes = 0;
    _pgrid = 0;
    
    polar(radius~radius1,theta1~theta1);
    
    k=k+.5;
    q=q+.05;
    
endfor;

