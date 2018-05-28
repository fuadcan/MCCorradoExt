   /*
** Procedure KPSSB:
**
** Purpose:
** Calculates the multivariate KPSS-test via bootstrapping.
**
** Format:
** { W0 , W1 } = KPSSB( Y , l );
**
** Input:
** Y    (TxN)-matrix, containing series.
** l    scalar, truncation parameter for Bartlett window.
**
** Output:
** W0   scalar, value of multivariate KPSS test for null hypothesis of zero
**          mean stationarity.
** W1   scalar, value of multivariate KPSS test for null hypothesis of level
**          stationarity.
*/


proc(2) = KPSSBBB( Y , l );

local psi1_0,psi1_1,T,N,S_0,S_1,p,Sdemean_0,Sdemean_1,W0B,W1B,W0,W1,i,j,e_0,e_1,eb_0,eb_1,R,W0BS,W1BS,W0PH,W1PH,ind,YB_0,YB_1,vind,ve_0,ve_1,veb_0,veb_1,
X, B,eOrig, SOrig, SdemeanOrig, psi1Orig, X_1, B_1,veyb_0,d,dmin,vey,eb_b,Yb,TT,dd, B_0,ind_0,ind_1,DY,vind_0,vind_1,ll,vcb,aic,sbc,B0,B1,e0,e1,Bb_1,Bb_0,c,ar,sbcm,ma,
eb_b1,Yb_L,Yb_LP,Yb_D,w,rows_storeB0,store_B0,c1,h;


library arima, pgraph;
arimaset;


@ ******************************** @
@ setting optimisation routine for arima@
@ ******************************** @

_am_opts[1]=0;

_am_opts[2]=0;


_am_itol[2] =1e-8;
_am_itol[3] =1e-6;



format /ld 8,4;

T = rows(Y);
N = cols(Y);

/*
LET h[3,2] = 1, 2, 3, 4, 5, 6;

h =   (abs(h) .>= 2).*0.00 + 
                 (abs(h) .<= 2).*h;

print "h" h; 
*/





@ ******************************** @
@ FIT FIRST DIFFERENCE OF Y WITH ARIMA(P,1,1) TO OBTAIN e_0 @
@ ******************************** @
ar=2;

c=1;

ma=1;

/*B_0=zeros(ar+ma+c,N); */

/* in the multivariate version we may have up to N cross-sectional units
so we need to create a homogeneous matrix for AR parameters. columns will contain coefficients
for different lag orders */

B_0=zeros(ar,N); 

e_0=zeros(T-1,N);

sbcm=zeros(ar,N);



/*this parts gives ARIMA estimates for all regions and up to max lag order ar.
This info is needed to select optmal lag length */

i=1;

j=1;

do until i gt N;
 do until j gt ar;


{B0,ll,e0,vcb,aic,sbc} = arima(0,Y[.,i],j,1,ma,c);

/*we store aic for different lag orders across cross-sectional units*/

sbcm[j,i]=aic;

j=j+1;



endo;

i=i+1;

endo;



@ ********************************                                  @
@ FINDS OPTIMAL LAG LENGTH  and GETS COEFFICINENTS B0 and ERROR E0  @
@ ********************************                                  @

/* sbcm is a AR*N so we find the minimum aic value across columns -that is 
for each cointegrating relationship of the multivariate version. We take transpose to
get a column vector
;
 */



d=(minindc(sbcm))';

/*
print "sbc" sbc;
print "rows(d)~cols(d)" rows(d)~cols(d);
print "rows(Y)~cols(Y)" rows(Y)~cols(Y);
*/


/* For each cointegrating relationship we estimate arima using optimal
AR lag length. We store the first AR coefficients for each column i.*/
/*print "B0"  B0[1:d[i],i];*/
/*store_B0 = B0[1:d[i],i];*/



i=1;

do until i gt N;

{B0,ll,e0,vcb,aic,sbc} = arima(0,Y[.,i],d[1,i],1,1,c);

/*when i is 2 then d still scalar and therefore d[1,i] inex out of range*/

/*
print "i,N" i~N;
print "d" d;
print "inside loop rows(Y)~cols(Y)" rows(Y)~cols(Y);
*/

/*storing first AR components of vector B0*/


store_B0 = B0[1:d[1,i]];
rows_storeB0 = rows(store_B0);
c1=1;
do until c1 gt rows_storeB0;

store_B0[c1] =   (abs(store_B0[c1]) .>= 0.999).*0.95 + 
                 (abs(store_B0[c1]) .<= 0.999).*store_B0[c1];
c1 = c1+1;
endo;


/*
store_B0 = B0[1:d[1,i]];
store_B0 =   (abs(store_B0) .>= 0.999).*0.00 + 
                 (abs(store_B0) .<= 0.999).*store_B0;
*/

B_0[1:d[1,i],i]=store_B0;

/*stores coeff in i-th column of B_0. If lag length p
is less than AR the last AR-p elemenst in each column 
of B_0 will be zero. This defines a matrix of coeff AR 
by N used for iteration
print "B0[1:d[1,i]]" B0[1:d[1,i]];
print "B_0[1:d[1,i],i]" B_0[1:d[1,i],i];

*/


/*print "cols(B_0)~rows(B_0)" cols(B_0)~rows(B_0);*/

e_0[.,i]=e0;


i=i+1;


endo; 


@ ******************************** @
@ Centering residuals??@
@ ******************************** @


e_0 = e_0 - (meanc(e_0)') ;




@ ******************************** @
@ Resampling via Bootstrapping@
@ ******************************** @



/*R = 200;*/

R=200;

W0B = zeros( R , 1 );
W0BS = zeros( R , 1 );
W1B = zeros( R , 1 ); 
W1BS = zeros( R , 1 );

rndseed 0.0001; 

j = 1;


@ VECTORISATION @



ve_0 = vecr(e_0);



do while (j < R);


@ GENERATING RANDOM NUMBERS for RESAMPLING @
    


ind_0 = ceil((T-1)*N*rndu(T-1,N));




vind_0 = vecr(ind_0);



@ RESAMPLING @


veb_0 = ve_0[vind_0];

/*SHALL ALSO COEFFICIENTS BE RESAMPLED??*/


@ THIS PARTS RECOVERS THE ORIGINAL MATRIX STRUCTURE for ERROR TERM  @

eb_0 = reshape(veb_0,T-1,N);


@ THIS PARTS GENERATES BOOTSTRAPPED SERIES and sets starting values @


Yb_L=zeros(T,N);
Yb_D=zeros(T-1,N);


/* starting value equal zero or sample value??*/
 
@ THIS PARTS DERIVES STATIONARY DELTAY FROM ARIMA -TREND = 0@


/* This part sets starting values for  Yb_D and Yb_L to zero */
Yb_D[1:ar,.]= (Y[2:ar+1,.]-Y[1:ar,.]).*0;
Yb_L[1:ar,.]=Yb_L[1:ar,.].*0;


/* Starting from the ar+1 component in Yb_D this part generates iteratively
Yb_D[t] by multiplying ARIMA coefficients in B_0
for elements Yb_D[t-1]up to Yb_D[t-ar]. Note that some coeff in B_0 will be zero 
as optimal lag length varies and maybe less than max lag length ar.
To correctly identify the autoregressive elements in Yb_D 
we use notation  Yb_D[i-ar:i-1,.].
The first row correponds to lag t-ar and the last to lag t-1
whereas the coeff in B_0 start from lag t-1 up to lag t-ar. 
So we need to reverse coefficients
in  B_0 */

i=ar+1;
do until i gt T-1;
Yb_D[i,.]= sumc((rev(B_0) .* Yb_D[i-ar:i-1,.]))' + 
                          (eb_0[i,.] - eb_0[i-1,.]);

i=i+1;
endo;


/* this parts recovers level of Yb_L from Yb_D by applying 
same iteration as described above*/


w=ar+1;
do until w gt T-1;
Yb_L[w,.]=Yb_L[w-1,.]+sumc((rev(B_0) .* Yb_D[w-ar:w-1,.]))'
                     + (eb_0[w,.] - eb_0[w-1,.]) ;

w=w+1;
endo;

 

/*
w=ar+1;
do until w gt T-1;
Yb_L[w,.]=Yb_L[w-1,.]+Yb_D[w-1,.];

w=w+1;
endo;
*/


/*
print "Yb_D[1:T-1,.]~Yb_L[1:T-1,.]~Y[1:T-1,.]" Yb_D[1:T-1,.]~Yb_L[1:T-1,.];

print "meanYb_L[1:T-1,.]~meanYb_L[1:T-1,.]~meanY[1:T-1,.]" meanc(Yb_D[1:T-1,.])'~meanc(Yb_L[1:T-1,.])';
*/




@ THIS PARTS REGRESSES BOOTSTRAPPED SERIES AGAINST A TREND and A CONSTANT 
(SEE PAG. 65 PAPER)  @

/*Is this correct? Shall we consider the null here or is it built in the
= test?*/


format /ld 8,4;




TT = rows(Yb_L);
X = ones(TT,1)~seqa(1,1,TT);
B = Yb_L/X;
eb_0 = Yb_L-X*B;



@ THIS PARTS BUILDS THE NEWEY-WEST ESTIMATOR, G,  UNDER THE NULL (SEE PAG.
 65 PAPER)  @

psi1_0 = NeweyWst( eb_0 , l );



@ THIS INVERTS THE MATRIX G for THE TWO TESTS DESCRIBED AT PAG 66  @

psi1_0 = invpd(psi1_0);



@ THIS DEFINES MATRIX S for PERFECT (PAG 65) and RELATIVE CONVERGENCE
 (PAG 66) @


S_0 = cumsumc(Yb_L);
Sdemean_1 = cumsumc( (Yb_L-meanc(Yb_L)') );

@ THIS BUILDS THE TESTS for  PERFECT and RELATIVE CONVERGENCE (PAG 66) @


W0B[j]=0;
W1B[j]=0;


i = 0;


do while (i < T-1);
   
i = i + 1;

locate 9,1;
format /ld 8,2;
        
   
W0B[j] = W0B[j]+S_0[i,.]*psi1_0*S_0[i,.]';

W1B[j]= W1B[j]+Sdemean_1[i,.]*psi1_0*Sdemean_1[i,.]';


endo;

W0B[j]= W0B[j]./((T-1)^2);
W1B[j]= W1B[j]./((T-1)^2);

j = j + 1;



endo;

W0BS = sortc(W0B,1);
W1BS = sortc(W1B,1);


W0PH=zeros(99,1);

W1PH=zeros(99,1);



p=1;
do while (p < 100);
W0PH[p]=W0BS[(R)*(p/99)];
W1PH[p]=W1BS[(R)*(p/99)];
p = p + 1;
endo;


/*print "WOPH~W1PH" W0PH~W1PH;*/


/*
print "WO W1 " W0~W1; 
print "WO W1 " W0~W1; 
"e eb " e[1:5,.]~eb[1:5,.];
print "WOPH~W1PH" WOPH~W1PH;
print "WOB~W1B" W0B~W1B;*/



retp(W0PH,W1PH);



endp;












