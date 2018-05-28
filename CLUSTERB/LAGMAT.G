proc(1) = lagmat(x);

local y;

y = shiftr( x' , 1 , miss(0,0) );
y = y';

retp(y);
endp;


