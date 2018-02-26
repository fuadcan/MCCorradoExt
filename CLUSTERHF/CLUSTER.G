/*
** Procedure cluster:
**
** Purpose:
** Performs cluster algorithm for August 1998 version of convergence paper.
**
** Format:
** { C0 , C1 } = cluster( Y , l , critpval );
**
** Input:
** Y        (TxK)-matrix, containing series in columns.
** l        scalar, bandwith parameter for Bartlett window.
** critpval scalar, critical p-value used as stopcriterion for algorithm.
**
** Output:
** C0       (Kx1)-vector, containing cluster numbers of corresponding series
**              for asymptotically perfect convergence.
** C1       (kx1)-vector, containing cluster numbers of corresponding series
**              for asymptotically relative convergence.
*/
Proc(2) = cluster( Y , l , critpval );

local asymdis0 , asymdis1 , C , K , T , M0 , M1 , p0mat , p1mat , i , j , n ,
        w0 , w1 , p0 , p1 , clussize , noclus , maxrow , maxcol ,
        cluscntr , alfa , Y1 , clusvec1 , clusvec2 , cluscomb , c0 , c1 ,
        clusmax , clusmin,p1matprint,p0matprint,maxpvalue ;

@ initializing parameters @
K = cols(Y);
T = rows(Y);
C = seqa(1,1,K);
M0 = zeros( k , k );
M1 = zeros( k , k );
p0mat = zeros( k , k );
p1mat = zeros( k , k );

@ loading matrices containing asymptotic distributions of tests @
asymdis0 = loadd("asymdst0");
asymdis1 = loadd("asymdst1");



ClS;
"Cluster algorithm v 3.0";
"=======================";
Locate 6,1; "Initializing....                       ";

i = 0;
n = 0;
do while ( i < k );
    i = i+1;

    j = 0;
    do while ( j < i );
        j = j+1;
        n = n+1;

        locate 9,1;
        format /ld 8,2;
        "combination " n " of " (0.5.*K.*(K+1));

        if ( i /= j );

            Y1 = Y[.,(i~j)];
            alfa = cointmat( cols( Y1 ) );
            Y1 = Y1*alfa';
            { w0 , w1 } = mvkpss( Y1 , l );
            { p0 , p1 } = pvalue( w0 , w1 , cols(Y1) , asymdis0 , asymdis1 );
           


            cluscomb = (i~j);
            clussize = cols(cluscomb);
            M0[cluscomb,cluscomb] = w0.*ones(clussize,clussize);
            M1[cluscomb,cluscomb] = w1.*ones(clussize,clussize);
            p0mat[cluscomb,cluscomb] = p0.*ones(clussize,clussize);
            p1mat[cluscomb,cluscomb] = p1.*ones(clussize,clussize);


p1matprint=p1mat; 
p0matprint=p0mat; 

            M0 = M0.*(eye(k).==0);
            M1 = M1.*(eye(k).==0);
            p0mat = p0mat.*(eye(k).==0);
            p1mat = p1mat.*(eye(k).==0);

        endif;

    endo; @ j @
endo; @ i @



maxpvalue=maxc(p1matprint');
print "p1mat" p1matprint[1:5,.];
print "maxpvalue" maxpvalue; 
print "maxindc" maxindc(p1matprint');

 
print "p0mat" p0matprint[1:5,.];
   
print "p1mat" p1matprint[1:5,.];


Locate 6,1; "Clustering... Asymptotically perfect convergence";
locate 9,1;
"                                                                         ";

noclus = k ;
do while (maxc(maxc(p0mat))>critpval) and ( noclus > 1 );
    noclus = noclus-1;

    locate 9,1;
    format /ld 8,2;
    noclus " clusters, out of " k " series ";

    if ( noclus > 1 );
    { maxrow , maxcol } = locmin( -p0mat );
    cluscntr = ( c .== c[maxrow] )+( c.==c[maxcol] );
    cluscntr = cluscntr.*seqa(1,1,k);
    cluscntr = delif( cluscntr , (cluscntr.==0) );
    clussize = rows(cluscntr);
    clusmax  = maxc(c[maxrow]|c[maxcol]);
    clusmin  = minc(c[maxrow]|c[maxcol]);

    c[cluscntr] = clusmin.*ones(clussize,1);
    c = c-(c.>clusmax);

    i = 0;
    do while ( i < noclus );
        i = i+1;

        if ( i /= clusmin );
            cluscntr = ( c.==clusmin )+
                       ( c.==i );
            cluscntr = cluscntr.*seqa(1,1,k);
            cluscntr = delif( cluscntr , (cluscntr.==0) );
            clussize = rows(cluscntr);

            y1 = y[.,cluscntr];
            alfa = cointmat( cols(y1) );
            y1 = y1*alfa';

            { w0 , w1 } = mvkpss( Y1 , l );
            { p0 , p1 } = pvalue( w0 , w1 , cols(Y1) , asymdis0 , asymdis1 );

            clusvec1 = ( c.==clusmin ).*seqa(1,1,k);
            clusvec1 = delif( clusvec1 , (clusvec1.==0) );
            clusvec2 = ( c.==i ).*seqa(1,1,k);
            clusvec2 = delif( clusvec2 , (clusvec2.==0) );

            M0[clusvec1,clusvec2] = w0.*ones(rows(clusvec1),rows(clusvec2));
            M0[clusvec2,clusvec1] = w0.*ones(rows(clusvec2),rows(clusvec1));
            M1[clusvec1,clusvec2] = w1.*ones(rows(clusvec1),rows(clusvec2));
            M1[clusvec2,clusvec1] = w1.*ones(rows(clusvec2),rows(clusvec1));
            p0mat[clusvec1,clusvec2]
                = p0.*ones(rows(clusvec1),rows(clusvec2));
            p0mat[clusvec2,clusvec1]
                = p0.*ones(rows(clusvec2),rows(clusvec1));
            p1mat[clusvec1,clusvec2]
                = p1.*ones(rows(clusvec1),rows(clusvec2));
            p1mat[clusvec2,clusvec1]
                = p1.*ones(rows(clusvec2),rows(clusvec1));
            M0[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            M1[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            p0mat[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            p1mat[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));

        endif;
    endo;
    else;
        c = ones(k,1);
    endif;
endo;

c0 = c;

Locate 6,1; "Clustering... Asymptotically relative convergence";
locate 9,1;
"                                                                         ";

do while (maxc(maxc(p1mat))>critpval) and ( noclus > 1 );
    noclus = noclus-1;

    locate 9,1;
    format /ld 8,2;
    noclus " clusters, out of " k " series ";

    if ( noclus > 1 );
    { maxrow , maxcol } = locmin( -p1mat );
    cluscntr = ( c .== c[maxrow] )+( c.==c[maxcol] );
    cluscntr = cluscntr.*seqa(1,1,k);
    cluscntr = delif( cluscntr , (cluscntr.==0) );
    clussize = rows(cluscntr);
    clusmax  = maxc(c[maxrow]|c[maxcol]);
    clusmin  = minc(c[maxrow]|c[maxcol]);

    c[cluscntr] = clusmin.*ones(clussize,1);
    c = c-(c.>clusmax);

    i = 0;
    do while ( i < noclus );
        i = i+1;

        if ( i /= clusmin );
            cluscntr = ( c.==clusmin )+
                       ( c.==i );
            cluscntr = cluscntr.*seqa(1,1,k);
            cluscntr = delif( cluscntr , (cluscntr.==0) );
            clussize = rows(cluscntr);

            y1 = y[.,cluscntr];
            alfa = cointmat( cols(y1) );
            y1 = y1*alfa';

            { w0 , w1 } = mvkpss( Y1 , l );
            { p0 , p1 } = pvalue( w0 , w1 , cols(Y1) , asymdis0 , asymdis1 );

            clusvec1 = ( c.==clusmin ).*seqa(1,1,k);
            clusvec1 = delif( clusvec1 , (clusvec1.==0) );
            clusvec2 = ( c.==i ).*seqa(1,1,k);
            clusvec2 = delif( clusvec2 , (clusvec2.==0) );

            M0[clusvec1,clusvec2] = w0.*ones(rows(clusvec1),rows(clusvec2));
            M0[clusvec2,clusvec1] = w0.*ones(rows(clusvec2),rows(clusvec1));
            M1[clusvec1,clusvec2] = w1.*ones(rows(clusvec1),rows(clusvec2));
            M1[clusvec2,clusvec1] = w1.*ones(rows(clusvec2),rows(clusvec1));
            p0mat[clusvec1,clusvec2]
                = p0.*ones(rows(clusvec1),rows(clusvec2));
            p0mat[clusvec2,clusvec1]
                = p0.*ones(rows(clusvec2),rows(clusvec1));
            p1mat[clusvec1,clusvec2]
                = p1.*ones(rows(clusvec1),rows(clusvec2));
            p1mat[clusvec2,clusvec1]
                = p1.*ones(rows(clusvec2),rows(clusvec1));

p1matprint=p1mat; 
p0matprint=p0mat; 

            M0[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            M1[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            p0mat[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));
            p1mat[clusvec1,clusvec1] = zeros(rows(clusvec1),rows(clusvec1));

        endif;
    endo;
    else;
        c = ones(k,1);
    endif;
endo;

c1 = c;

retp( c0 , c1 );
endp;

