/*
** Procedure cointmat:
**
** Purpose:
** Constructs matrix with cointegrating vectors for cluster.
**
** Format:
** alfa = cointmat( k );
**
** Input:
** k        scalar, number of series in cluster.
**
** Output:
** alfa     (k-1xk)-matrix, with corresponding cointegrating vectors.
*/
proc(1) = cointmat( k );

local alfa ;

alfa = zeros( k-1 , k );
alfa[.,1:k-1] = alfa[.,1:k-1]+eye(k-1);
alfa[.,2:k] = alfa[.,2:k]-eye(k-1);

retp(alfa);
endp;


