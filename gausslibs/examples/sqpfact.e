
#include sqpsolvemt.sdf



lambda = { 1.0  0.0,
           0.5  0.0,
           0.0  1.0,
           0.0  0.5 };

lmask  = {  0    0,
            1    0,
            0    0,
            0    1 };

phi = { 1.0  0.3,
        0.3  1.0 };

theta = { 0.6  0.0  0.0  0.0,
          0.0  0.6  0.0  0.0,
          0.0  0.0  0.6  0.0,
          0.0  0.0  0.0  0.6 };

tmask = {  1    0    0    0,
           0    1    0    0,
           0    0    1    0,
           0    0    0    1 };

struct PV par0;                 
par0 = pvCreate;
par0 = pvPackm(par0,lambda,"lambda",lmask);    
par0 = pvPacks(par0,phi,"phi");
par0 = pvPacksm(par0,theta,"theta",tmask);     
    

struct SQPsolveMTControl c0;
c0 = sqpSolveMTcontrolCreate;      

lind = pvGetIndex(par0,"lambda"); /* get indices of lambda parameters */
                                  /* in parameter vector              */
tind = pvGetIndex(par0,"theta");  /* get indices of theta  parameters */ 
                                  /* in parameter vector              */ 

c0.bounds = ones(pvLength(par0),1).*(-1e250~1e250);

c0.bounds[lind,1] = zeros(rows(lind),1);
c0.bounds[lind,2] = 10*ones(rows(lind),1);
c0.bounds[tind,1] = .001*ones(rows(tind),1);
c0.bounds[tind,2] = 100*ones(rows(tind),1);



c0.output = 1;
c0.printIters = 1; 
c0.trustRadius = 1;
c0.ineqProc = &ineq;    
c0.covType = 1;
                                 
struct DS d0;
d0 = dsCreate; 
d0.dataMatrix = loadd("maxfact");

output file = sqpfact.out reset;


struct SQPsolveMTOut out0;
out0 = SQPsolveMT(&lpr,par0,d0,c0);

lambdahat = pvUnpack(out0.par,"lambda");
phihat = pvUnpack(out0.par,"phi");
thetahat = pvUnpack(out0.par,"theta");               


             
print "estimates";
print;
print "lambda" lambdahat;
print;
print "phi" phihat;
print;
print "theta" thetahat;
                      

struct PV stderr;
stderr = out0.par;  

if not scalmiss(out0.moment);
    stderr = pvPutParVector(stderr,sqrt(diag(out0.moment)));
                                           
   lambdase = pvUnpack(stderr,"lambda");
   phise = pvUnpack(stderr,"phi");
   thetase = pvUnpack(stderr,"theta");        

   print "standard errors";
   print;
   print "lambda" lambdase;
   print;
   print "phi" phise;
   print;
   print "theta" thetase;
endif;


output off;


proc lpr(struct PV par1, struct DS data1);  
    local lambda,phi,theta,sigma,logl;

    lambda = pvUnpack(par1,"lambda");
    phi = pvUnpack(par1,"phi");
    theta = pvUnpack(par1,"theta");
    sigma = lambda*phi*lambda' + theta;
                  
    logl = -lnpdfmvn(data1.dataMatrix,sigma);
    retp(logl);

endp;            



proc ineq(struct PV par1, struct DS data1);  
    local lambda,phi,theta,sigma,e;

    lambda = pvUnpack(par1,"lambda");
    phi = pvUnpack(par1,"phi");
    theta = pvUnpack(par1,"theta");
    sigma = lambda*phi*lambda' + theta;
    e = eigh(sigma) - .001; /* eigenvalues of sigma */
    e = e | eigh(phi) - .001;   /* eigenvalues of phi */
    retp(e);

endp;            

