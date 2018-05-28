/*
**                   dloop3.e 
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

mv = {.};

dataloop freqdata dumy;
    extern mv;
    code cage with
        1  for  age >= 6,
        0  for age < 6,
        mv for age $== mv;

    recode  sex with
       1  for sex $== "M",
       0  for sex $== "F",
       mv for sex $== mv;
    keep cage sex pay age;
endata;

