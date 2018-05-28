/*
**                   dloop2.e 
**
**
**  This example requires that You turn on dataloop translation 
**  before continuing. In Windows go to *Action* there enable
**  Translate Data Loop commands*
**
**  In tgauss type in the command *config* next select *Run options* 
**  next select *Translator* and toggle it to *On*, Then Quit config 
**  NOTE This step must be completed prior to running this example";
**
*/


cls;

disable;      /* because the data set contains missing data */

base = 1;     /* base for counter variable */

dataloop freqdata dd1;

    /* create counter variable -- the data for each iteration of the  */
    /*                            read loop is in x_x                 */
    extern base;
    make count = seqa(base,1,rows(x_x));
    # base = base + rows(x_x);

    make y1 = age^2 ;
    select 13 < count and count < 45;
    drop wt sex;
endata;

datalist dd1;


