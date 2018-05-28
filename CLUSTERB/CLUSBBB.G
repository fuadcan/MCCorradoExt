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
** Output:  @

** C0       (Kx1)-vector, containing cluster numbers of corresponding series
**              for asymptotically perfect convergence.
** C1       (kx1)-vector, containing cluster numbers of corresponding series
**              for asymptotically relative convergence.
*/
proc(2) = clusbbb( Y , l , critpval );

local W0PH , W1PH ,asymdis0, asymdis1, C , K , T , M0 , M1 , p0mat , p1mat , i , j , n ,
        w0 , w1 , p0 , p1 , clussize , noclus , maxrow , maxcol ,
        cluscntr , alfa , Y1 , clusvec1 , clusvec2 , cluscomb , c0 , c1 ,
        clusmax , clusmin, maxpvalue, p1matprint,p0matprint ;

@ initializing parameters @
K = cols(Y);


T = rows(Y);
C = seqa(1,1,K);        /*creates an addtive sequence from 1 to K with increments 1*/
M0 = zeros( k , k );    /*CREATES A MATRIX K BY K OF ZERO WHERE K IS THE NUMBER 
 OF REGIONS*/
M1 = zeros( k , k );    /*CREATES A MATRIX K BY K OF ZERO WHERE K IS THE NUMBER
 OF REGIONS*/
p0mat = zeros( k , k ); /*CREATES A MATRIX K BY K OF ZERO WHERE K IS THE NUMBER 
OF REGIONS*/
p1mat = zeros( k , k ); /*CREATES A MATRIX K BY K OF ZERO WHERE K IS THE NUMBER 
OF REGIONS*/

/*@ loading matrices containing asymptotic distributions of tests- we do not
 need this with bootstrapped values @
@asymdis0 = loadd("asymdst0");@
@asymdis1 = loadd("asymdst1");@
*/


ClS;


"Cluster algorithm TEST TTTTTTT v 3.0";
"=======================";
Locate 6,1; "Initializing....Test                       ";

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

            Y1 = Y[.,(i~j)];              /*THIS EXCTRACTS FROM Y THE SERIES 
                                                OF OUTPUT for REGIONS I-J 
                                                AND COMPACT THEM COLUMNWISE*/

            

alfa = cointmat( cols( Y1 ) );     /*CONSTRUCTS MATRIX WITH COINTEGRATING
                                                                    VECTORS*/



            Y1 = Y1*alfa';            /*MULTIPLIES THE PAIRWSIE SERIES for
                                                 THE COINTEGRATING VECTOR*/
              
            
            
            { w0 , w1 } = mvkpss( Y1 , l );    /*CALCULATES THE MKPSS TEST
                                                                      for BANDWIDTH L*/

            {W0PH, W1PH} = kpssbbb( Y1 , l );  /*CALCULATES THE CRITICAL VALUES for 
                                                        W0  and W1 VIA BOOTSTRAPPING*/
  
            { p0 , p1 } = pvalue( W0, W1 , cols(Y1) , W0PH , W1PH ); 
                                           /*RETURNS THE P-VALUES OF THE TEST*/
            
           cluscomb = (i~j);        /*CONSIDERS THE COMBINATION OF REGIONS I and J*/
       
           clussize = cols(cluscomb);  /*GIVES THE NUMBER OF COLUMNS OF CLUSCOMB. 
                                             AS THEY ARE PAIRWISE COMPARISON SIZE
                                                      IS 2*/
           
            M0[cluscomb,cluscomb] = W0.*ones(clussize,clussize); 
                  /*CREATES A MATRIX 2 BY 2 WITH EQUAL VALUES W0
                and INSERTS THEM IN M0 AT POSITIONS II,IJ,JI,JJ*/

            M1[cluscomb,cluscomb] = W1.*ones(clussize,clussize);     /*IDEM for M1*/
            p0mat[cluscomb,cluscomb] = p0.*ones(clussize,clussize);  /*IDEM for P0*/
            p1mat[cluscomb,cluscomb] = p1.*ones(clussize,clussize);  /*IDEM for P1*/


p1matprint=p1mat; 
p0matprint=p0mat; 


            M0 = M0.*(eye(k).==0);    
      /*eye(k).==0 CRETAES A MATRIX OF ONES  EXCEPT THAT ON THE DIAGONAL 
       AND  MULTIPLIES THAT for MO.   IT RETUNR A MATRIX OF THE TEST W0 
      AT POSITION  IJ and JI (symmetric)*/                    
        
           
            M1 = M1.*(eye(k).==0);                    /*IDEM for M1*/

            p0 = p0mat.*(eye(k).==0);       /*IDEM for P0MAT*/
         
            p1 = p1mat.*(eye(k).==0);        /*IDEM for P1MAT*/
          
 
        endif;

    endo; @ j @
endo; @ i @



/*maxpvalue=maxc(p1matprint');

print "p1mat" p1matprint[1:4,.];
print "maxpvalue" maxpvalue; 
print "maxindc" maxindc(p1matprint');

 
print "p0mat" p0matprint[1:4,.];
   
print "p1mat" p1matprint[1:4,.];*/



Locate 6,1; "Clustering... Asymptotically perfect convergence";
locate 9,1;
"                                                                         ";
noclus = k ;
do while (maxc(maxc(p0mat))>critpval) and ( noclus > 1 ); 

      /*takes the max of p0mat and if inside critical region 
       you continue the pairwise comparison.  it exists at 
       least one pair i and j  which are clustering*/

print  "rowsp0mat colsp0mat" rows(p0mat)~cols(p0mat);

    noclus = noclus-1;
    locate 9,1;
    format /ld 8,0;
    noclus " clusters, out of " k " series ";
    if ( noclus > 1 );

    /*{ maxrow , maxcol } = locmin( -p0mat );*/
    { maxrow , maxcol } = locmin( (-p0mat + eye(k).*999));  

    /* Gives row and column of minimum element in matrix 

       If minimum is not unique it returns the minimum element with the
       lowest column and if this column contains 2 elements equal to the
       minimum it returns the lowest row */

    cluscntr = ( c .== c[maxrow] )+( c.==c[maxcol] ); 

     /*identifies elements with minimum p0mat. 
       c is a row vector containing the indices of different regions

      It will create a vector with 2 elements = to 1 and 0 else.
      if maxrow=maxcol then there willbe one element = to 2 else 0
      This will then generate an error im code below
      */

    cluscntr = cluscntr.*seqa(1,1,k);
 
   /*creates a vector of zeros except for the region numbers corresponding 
   to pmin*/

    cluscntr = delif( cluscntr , (cluscntr.==0) ); 

    /*delete rows from cluscntr which are equal to zero*/

    clussize = rows(cluscntr); 

    /*gives the cluster size*/

    clusmax  = maxc(c[maxrow]|c[maxcol]); 
/*identifies element with max index*/

    clusmin  = minc(c[maxrow]|c[maxcol]);  
/*identifies element with min index*/

 print "clusmin clusmax clussize" clusmin~clusmax~clussize;

 print "rows cols c, rows cols cluscntr " rows(c)~cols(c)~rows(cluscntr)~cols(cluscntr);

    c[cluscntr] = clusmin.*ones(clussize,1); 
/*attributes to regions  the same cluster number*/
     
     
    
    c = c-(c.>clusmax);  
  /*creates a vector where regions belonging to the same cluster have 
   same indicator number*/


    
    i = 0;
    do while ( i < noclus );
        i = i+1;

        if ( i /= clusmin );  
         /*it considers single regions not initially clustered and we pose 
          i in clusmin*/

            cluscntr = ( c.==clusmin )+  ( c.==i );
                  /* creates a new cluster matching regions in 
                   clusmin with region i*/
                      
            cluscntr = cluscntr.*seqa(1,1,k); 
                   /*gives row of new bivariate comparison*/
             
            cluscntr = delif( cluscntr , (cluscntr.==0) );
                 /*deletes rows which are not zero*/
            
             clussize = rows(cluscntr);          
                /*gives cluster size*/
             
           y1 = y[.,cluscntr];   
              /*extracts series corresponding to the new pairwise comparison*/
           
          alfa = cointmat( cols(y1) );             
             /*creates cointegration vector*/
          
           y1 = y1*alfa';                        
          /*maps into cointegrating series*/

           { w0 , w1 }    = mvkpss( Y1 , l );       

            /*calculates MKPSS test*/

           { W0PH, W1PH } = kpssbbb( Y1 , l );        
            /*gives bootstrapped critical values*/
          
            { p0 , p1 } = pvalue( W0, W1, cols(Y1) , W0PH , W1PH ); 
            /*calculates p-values*/
          

            clusvec1 = ( c.==clusmin ).*seqa(1,1,k);       
             /*identifies elements from clusmin and i*/

            clusvec1 = delif( clusvec1 , (clusvec1.==0) );
            clusvec2 = ( c.==i ).*seqa(1,1,k);
            clusvec2 = delif( clusvec2 , (clusvec2.==0) );

            M0[clusvec1,clusvec2] = W0.*ones(rows(clusvec1),rows(clusvec2));
            M0[clusvec2,clusvec1] = W0.*ones(rows(clusvec2),rows(clusvec1));
            M1[clusvec1,clusvec2] = W1.*ones(rows(clusvec1),rows(clusvec2));
            M1[clusvec2,clusvec1] = W1.*ones(rows(clusvec2),rows(clusvec1));
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
    format /ld 8,0;
    noclus " clusters, out of " k " series ";

    if ( noclus > 1 );
  /*  { maxrow , maxcol } = locmin( -p1mat );*/

{ maxrow , maxcol } = locmin( (-p1mat + eye(k).*999)); 

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

            { W0PH, W1PH } = kpssbbb( Y1 , l );

            { p0 , p1 } = pvalue( W0, W1, cols(Y1) , W0PH , W1PH );

            clusvec1 = ( c.==clusmin ).*seqa(1,1,k);
            clusvec1 = delif( clusvec1 , (clusvec1.==0) );
            clusvec2 = ( c.==i ).*seqa(1,1,k);
            clusvec2 = delif( clusvec2 , (clusvec2.==0) );

            M0[clusvec1,clusvec2] = W0.*ones(rows(clusvec1),rows(clusvec2));
            M0[clusvec2,clusvec1] = W0.*ones(rows(clusvec2),rows(clusvec1));
            M1[clusvec1,clusvec2] = W1.*ones(rows(clusvec1),rows(clusvec2));
            M1[clusvec2,clusvec1] = W1.*ones(rows(clusvec2),rows(clusvec1));
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




c1 = c;

retp( c0 , c1 );


endp;


/*print "c" c;
print "c .== c[maxrow]" (c .== c[maxrow]);
print " c.==c[maxcol]"  (c.==c[maxcol]);
print "cluscntr" cluscntr;
print "seqa(1,1,k)" seqa(1,1,k);
print "cluscntr = delif( cluscntr , (cluscntr.==0) )"  delif( cluscntr , (cluscntr.==0) ); 
print "clussize" clussize;
print "clusmax" clusmax;
print "clusmin" clusmin;
print "c[cluscntr]" c[cluscntr];
print "c.>clusmax" c.>clusmax;
print "c-(c.>clusmax)" c-(c.>clusmax);
print "maxc(maxc(p0mat))" maxc(maxc(p0mat));
print "locmin( -p0mat )" locmin( -p0mat );
print "M0[cluscomb,cluscomb]" M0[cluscomb,cluscomb];    */

