/*
**     dloop4.e  example censors the variable job ifjob < 2.7 
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

dataloop sci dd1;
    vector k = 1;
    recode job with
       0 for job < 2.7;
    select female /= k;
    keep job sex female k;
endata;

