/*
**                   dloop1.e 
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

m = {.};
i = 1;
dataloop sci dd1;
    extern m, i;
    vector cnst = i;
    select enrol $/= m;
    keep sex job cnst;
endata;
          
datalist dd1;

