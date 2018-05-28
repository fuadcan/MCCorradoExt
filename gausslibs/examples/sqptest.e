#include sqpsolvemt.sdf



     struct DS Y;                
     Y = dsCreate;         

     Y.dataMatrix = 3.183|
                 3.059|
                 2.871|
                 2.622|
                 2.541|
                 2.184|
                 2.110|
                 2.075|
                 2.018|
                 1.903|
                 1.770|
                 1.762|
                 1.550;

     struct DS X;
     X = dsCreate;

     X.dataMatrix = seqa(1,1,13);

     struct DS Z;
     Z = reshape(Z,2,1);
     Z[1] = Y;
     Z[2] = X;

                                             
      struct SQPsolveMTControl c1;
      c1 = sqpSolveMTcontrolCreate;  /* initializes default values */
      c1.bounds = 0~100;             /* constrains parameters */
                                     /* to be positive        */

      c1.CovType = 1;

      c1.output = 1;
      c1.printIters = 0;             
      c1.gradProc = &grad;
      
      struct PV par1;   
      par1 = pvCreate;

//      start = { 1, 2.5, .1 };
      start = { 2, 4, 2 };
      par1 = pvPack(par1,start,"Parameters");




      output file = sqptest1.out reset;

      struct SQPsolveMTout out1;              
      out1 = SQPsolveMT(&Micherlitz,par1,Z,c1);
      
      estimates = pvGetParVector(out1.par);                         
      print " parameter estimates ";

      print estimates;
      print;
      print " standard errors ";
      print sqrt(diag(out1.moment));
 
      output off;

      proc Micherlitz(struct PV par1,struct DS Z);
           local p0,e,s2;
	   struct DS Y;
	   struct DS X;
           p0 = pvUnpack(par1,"Parameters");
           y = Z[1];
           x = Z[2];
           e = y.dataMatrix - p0[1] - p0[2]*exp(-p0[3]*x.dataMatrix);
           s2 = moment(e,0)/(rows(e)-1);         
           retp( (2/rows(e))*(e.*e/s2 + ln(2*pi*s2)));
      endp;


      proc grad(struct PV par1, struct DS Z);
           local p0,e,e1,e2,e3,w,g,s2;  
      
           p0 = pvUnpack(par1,"Parameters");
           y = Z[1];
           x = Z[2];

           w = exp(-p0[3]*x.dataMatrix);
           e = y.dataMatrix - p0[1] - p0[2]*w;
           s2 = moment(e,0)/rows(e);
           
           e1 = -ones(rows(e),1);
           e2 = -w;
           e3 = p0[2]*x.dataMatrix.*w;             

           w = (1 - e.*e / s2) / rows(e);
           g = e.*e1 + w*(e'e1);
           g = g ~ (e.*e2 + w*(e'e2));
           g = g ~ (e.*e3 + w*(e'e3));    
                             
           retp(4*g/(rows(e)*s2));
      endp;


