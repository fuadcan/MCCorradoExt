if(Sys.info()[1]=="Linux"){setwd("~/Documents/MCCorradoExt/")} else {setwd("~/MCCorradoExt/")}
source("dataGen.R")
source("mcAGK.R")
source("mcHF.R")
source("mcCW.R")
source("mcCWplus.R")
source("mcHFplus.R")
source("mcCWplus.R")
# source("anlysSingle.R")
# source("anlysMulti.R")
# source("mclapplyhack.R")
library("parallel")

numofrep <- 100
# commands for generating data single club data:
# no intercept
dataGen(50,10,3,0.2,numofrep,T)
dataGen(100,10,3,0.2,numofrep,T)
dataGen(50,10,3,0.6,numofrep,T)
dataGen(100,10,3,0.6,numofrep,T)
# with intercept
dataGen(50,10,3,0.2,numofrep,F)
dataGen(100,10,3,0.2,numofrep,F)
dataGen(50,10,3,0.6,numofrep,F)
dataGen(100,10,3,0.6,numofrep,F)

## m=5 ##
dataGen(50,10,5,0.2,numofrep,T)
dataGen(100,10,5,0.2,numofrep,T)
dataGen(50,10,5,0.6,numofrep,T)
dataGen(100,10,5,0.6,numofrep,T)

# with intercept
dataGen(50,10,5,0.2,numofrep,F)
dataGen(100,10,5,0.2,numofrep,F)
dataGen(50,10,5,0.6,numofrep,F)
dataGen(100,10,5,0.6,numofrep,F)

## N = 20 ##
## m=3 ##
dataGen(50,20,3,0.2,numofrep,T)
dataGen(100,20,3,0.2,numofrep,T)
dataGen(50,20,3,0.6,numofrep,T)
dataGen(100,20,3,0.6,numofrep,T)
# with intercept
dataGen(50,20,3,0.2,numofrep,F)
dataGen(100,20,3,0.2,numofrep,F)
dataGen(50,20,3,0.6,numofrep,F)
dataGen(100,20,3,0.6,numofrep,F)
## m=5 ##
dataGen(50,20,5,0.2,numofrep,T)
dataGen(100,20,5,0.2,numofrep,T)
dataGen(50,20,5,0.6,numofrep,T)
dataGen(100,20,5,0.6,numofrep,T)
# with intercept
dataGen(50,20,5,0.2,numofrep,F)
dataGen(100,20,5,0.2,numofrep,F)
dataGen(50,20,5,0.6,numofrep,F)
dataGen(100,20,5,0.6,numofrep,F)
## m=7 ##
dataGen(50,20,7,0.2,numofrep,T)
dataGen(100,20,7,0.2,numofrep,T)
dataGen(50,20,7,0.6,numofrep,T)
dataGen(100,20,7,0.6,numofrep,T)
# with intercept
dataGen(50,20,7,0.2,numofrep,F)
dataGen(100,20,7,0.2,numofrep,F)
dataGen(50,20,7,0.6,numofrep,F)
dataGen(100,20,7,0.6,numofrep,F)
## m=10 ##
dataGen(50,20,10,0.2,numofrep,T)
dataGen(100,20,10,0.2,numofrep,T)
dataGen(50,20,10,0.6,numofrep,T)
dataGen(100,20,10,0.6,numofrep,T)
# with intercept
dataGen(50,20,10,0.2,numofrep,F)
dataGen(100,20,10,0.2,numofrep,F)
dataGen(50,20,10,0.6,numofrep,F)
dataGen(100,20,10,0.6,numofrep,F)

#### MULTI CLUB CASE ####
## k=2
# no intercept
dataGenplus(50,10,2,0.2,numofrep,T)
dataGenplus(100,10,2,0.2,numofrep,T)
dataGenplus(50,10,2,0.6,numofrep,T)
dataGenplus(100,10,2,0.6,numofrep,T)
# with intercept
dataGenplus(50,10,2,0.2,numofrep,F)
dataGenplus(100,10,2,0.2,numofrep,F)
dataGenplus(50,10,2,0.6,numofrep,F)
dataGenplus(100,10,2,0.6,numofrep,F)
## k = 3
# no intercept
dataGenplus(50,10,3,0.2,numofrep,T)
dataGenplus(100,10,3,0.2,numofrep,T)
dataGenplus(50,10,3,0.6,numofrep,T)
dataGenplus(100,10,3,0.6,numofrep,T)
# with intercept
dataGenplus(50,10,3,0.2,numofrep,F)
dataGenplus(100,10,3,0.2,numofrep,F)
dataGenplus(50,10,3,0.6,numofrep,F)
dataGenplus(100,10,3,0.6,numofrep,F)
## k=4
# no intercept
dataGenplus(50,20,4,0.2,numofrep,T)
dataGenplus(100,20,4,0.2,numofrep,T)
dataGenplus(50,20,4,0.6,numofrep,T)
dataGenplus(100,20,4,0.6,numofrep,T)
# with intercept
dataGenplus(50,20,4,0.2,numofrep,F)
dataGenplus(100,20,4,0.2,numofrep,F)
dataGenplus(50,20,4,0.6,numofrep,F)
dataGenplus(100,20,4,0.6,numofrep,F)
## k = 5
# no intercept
dataGenplus(50,20,5,0.2,numofrep,T)
dataGenplus(100,20,5,0.2,numofrep,T)
dataGenplus(50,20,5,0.6,numofrep,T)
dataGenplus(100,20,5,0.6,numofrep,T)
# with intercept
dataGenplus(50,20,5,0.2,numofrep,F)
dataGenplus(100,20,5,0.2,numofrep,F)
dataGenplus(50,20,5,0.6,numofrep,F)
dataGenplus(100,20,5,0.6,numofrep,F)

#########################################################
# Commands for detection of clubs via maximal clique algorithm
# no intercept
mcAGK(50,10,3,0.2,T)
mcAGK(100,10,3,0.2,T)
mcAGK(50,10,3,0.6,T)
mcAGK(100,10,3,0.6,T)
# with intercept
mcAGK(50,10,3,0.2,F)
mcAGK(100,10,3,0.2,F)
mcAGK(50,10,3,0.6,F)
mcAGK(100,10,3,0.6,F)
#
mcAGK(50,10,5,0.2,T)
mcAGK(100,10,5,0.2,T)
mcAGK(50,10,5,0.6,T)
mcAGK(100,10,5,0.6,T)
# with intercept
mcAGK(50,10,5,0.2,F)
mcAGK(100,10,5,0.2,F)
mcAGK(50,10,5,0.6,F)
mcAGK(100,10,5,0.6,F)

# no intercept
mcAGK(50,20,3,0.2,T)
mcAGK(100,20,3,0.2,T)
mcAGK(50,20,3,0.6,T)
mcAGK(100,20,3,0.6,T)
# with intercept
mcAGK(50,20,3,0.2,F)
mcAGK(100,20,3,0.2,F)
mcAGK(50,20,3,0.6,F)
mcAGK(100,20,3,0.6,F)
#
mcAGK(50,20,5,0.2,T)
mcAGK(100,20,5,0.2,T)
mcAGK(50,20,5,0.6,T)
mcAGK(100,20,5,0.6,T)
# with intercept
mcAGK(50,20,5,0.2,F)
mcAGK(100,20,5,0.2,F)
mcAGK(50,20,5,0.6,F)
mcAGK(100,20,5,0.6,F)
# m=7
mcAGK(50,20,7,0.2,T)
mcAGK(100,20,7,0.2,T)
mcAGK(50,20,7,0.6,T)
mcAGK(100,20,7,0.6,T)
# with intercept
mcAGK(50,20,7,0.2,F)
mcAGK(100,20,7,0.2,F)
mcAGK(50,20,7,0.6,F)
mcAGK(100,20,7,0.6,F)
# m=10
mcAGK(50,20,10,0.2,T)
mcAGK(100,20,10,0.2,T)
mcAGK(50,20,10,0.6,T)
mcAGK(100,20,10,0.6,T)
# with intercept
mcAGK(50,20,10,0.2,F)
mcAGK(100,20,10,0.2,F)
mcAGK(50,20,10,0.6,F)
mcAGK(100,20,10,0.6,F)
# Commands for detection of clubs via HF algorithm 
# no intercept
mcHF(50,10,3,0.2,T)
mcHF(100,10,3,0.2,T)
mcHF(50,10,3,0.6,T)
mcHF(100,10,3,0.6,T)
#
mcHF(50,10,3,0.2,F)
mcHF(100,10,3,0.2,F)
mcHF(50,10,3,0.6,F)
mcHF(100,10,3,0.6,F)
## m=5
mcHF(50,10,5,0.2,T)
mcHF(100,10,5,0.2,T)
mcHF(50,10,5,0.6,T)
mcHF(100,10,5,0.6,T)
#
mcHF(50,10,5,0.2,F)
mcHF(100,10,5,0.2,F)
mcHF(50,10,5,0.6,F)
mcHF(100,10,5,0.6,F)
# N=20
# m=3
# no intercept
mcHF(50,20,3,0.2,T)
mcHF(100,20,3,0.2,T)
mcHF(50,20,3,0.6,T)
mcHF(100,20,3,0.6,T)
# with intercept
mcHF(50,20,3,0.2,F)
mcHF(100,20,3,0.2,F)
mcHF(50,20,3,0.6,F)
mcHF(100,20,3,0.6,F)
# m=5
mcHF(50,20,5,0.2,T)
mcHF(100,20,5,0.2,T)
mcHF(50,20,5,0.6,T)
mcHF(100,20,5,0.6,T)
# with intercept
mcHF(50,20,5,0.2,F)
mcHF(100,20,5,0.2,F)
mcHF(50,20,5,0.6,F)
mcHF(100,20,5,0.6,F)
# m=7
mcHF(50,20,7,0.2,T)
mcHF(100,20,7,0.2,T)
mcHF(50,20,7,0.6,T)
mcHF(100,20,7,0.6,T)
# with intercept
mcHF(50,20,7,0.2,F)
mcHF(100,20,7,0.2,F)
mcHF(50,20,7,0.6,F)
mcHF(100,20,7,0.6,F)
# m=10
mcHF(50,20,10,0.2,T)
mcHF(100,20,10,0.2,T)
mcHF(50,20,10,0.6,T)
mcHF(100,20,10,0.6,T)
# with intercept
mcHF(50,20,10,0.2,F)
mcHF(100,20,10,0.2,F)
mcHF(50,20,10,0.6,F)
mcHF(100,20,10,0.6,F)

####  CW  ####
mcCW(50,10,3,0.2,T)
mcCW(100,10,3,0.2,T)
mcCW(50,10,3,0.6,T)
mcCW(100,10,3,0.6,T)
#
mcCW(50,10,3,0.2,F)
mcCW(100,10,3,0.2,F)
mcCW(50,10,3,0.6,F)
mcCW(100,10,3,0.6,F)
#
mcCW(50,10,5,0.2,T)
mcCW(100,10,5,0.2,T)
mcCW(50,10,5,0.6,T)
mcCW(100,10,5,0.6,T)
#
mcCW(50,10,5,0.2,F)
mcCW(100,10,5,0.2,F)
mcCW(50,10,5,0.6,F)
mcCW(100,10,5,0.6,F)

####### MULTI CLUB CASE ########
mcCWplus(50,10,2,.2,T,F)
mcCWplus(100,10,2,.2,T,F)
mcCWplus(50,10,2,.6,T,F)
mcCWplus(100,10,2,.6,T,F)
# with constant
mcCWplus(50,10,2,.2,F,F)
mcCWplus(100,10,2,.2,F,F)
mcCWplus(50,10,2,.6,F,F)
mcCWplus(100,10,2,.6,F,F)
# m=3
mcCWplus(50,10,3,.2,T,F)
mcCWplus(50,10,3,.6,T,F)
mcCWplus(100,10,3,.2,T,F)
mcCWplus(100,10,3,.6,T,F)
# with constant
mcCWplus(50,10,3,.2,F,F)
mcCWplus(50,10,3,.6,F,F)
mcCWplus(100,10,3,.2,F,F)
mcCWplus(100,10,3,.6,F,F)
# m=4
mcCWplus(50,20,4,.2,T,F)
mcCWplus(50,20,4,.6,T,F)
mcCWplus(100,20,4,.2,T,F)
mcCWplus(100,20,4,.6,T,F)
# with constant
mcCWplus(50,20,4,.2,F,F)
mcCWplus(50,20,4,.6,F,F)
mcCWplus(100,20,4,.2,F,F)
mcCWplus(100,20,4,.6,F,F)
# m=5
mcCWplus(50,20,5,.2,T,F)
mcCWplus(50,20,5,.6,T,F)
mcCWplus(100,20,5,.2,T,F)
mcCWplus(100,20,5,.6,T,F)
# with constant
mcCWplus(50,20,5,.2,F,F)
mcCWplus(50,20,5,.6,F,F)
mcCWplus(100,20,5,.2,F,F)
mcCWplus(100,20,5,.6,F,F)

####  HF  ####
mcHFplus(50,10,2,.2,T,F)
mcHFplus(100,10,2,.2,T,F)
mcHFplus(50,10,2,.6,T,F)
mcHFplus(100,10,2,.6,T,F)
# with constant
mcHFplus(50,10,2,.2,F,F)
mcHFplus(100,10,2,.2,F,F)
mcHFplus(50,10,2,.6,F,F)
mcHFplus(100,10,2,.6,F,F)
# m=3
mcHFplus(50,10,3,.2,T,F)
mcHFplus(50,10,3,.6,T,F)
mcHFplus(100,10,3,.2,T,F)
mcHFplus(100,10,3,.6,T,F)
# with constant
mcHFplus(50,10,3,.2,F,F)
mcHFplus(50,10,3,.6,F,F)
mcHFplus(100,10,3,.2,F,F)
mcHFplus(100,10,3,.6,F,F)
# m=4
mcHFplus(50,20,4,.2,T,F)
mcHFplus(50,20,4,.6,T,F) # from here 12:21
mcHFplus(100,20,4,.2,T,F)
mcHFplus(100,20,4,.6,T,F)
# with constant
mcHFplus(50,20,4,.2,F,F)
mcHFplus(50,20,4,.6,F,F)
mcHFplus(100,20,4,.2,F,F)
mcHFplus(100,20,4,.6,F,F)
# m=5
mcHFplus(50,20,5,.2,T,F)
mcHFplus(50,20,5,.6,T,F)
mcHFplus(100,20,5,.2,T,F)
mcHFplus(100,20,5,.6,T,F)
# with constant
mcHFplus(50,20,5,.2,F,F)
mcHFplus(50,20,5,.6,F,F)
mcHFplus(100,20,5,.2,F,F)
mcHFplus(100,20,5,.6,F,F)

####### CW #######
# k=2
mcCWplus(50,10,2,.2,T,F)
mcCWplus(100,10,2,.2,T,F)
mcCWplus(50,10,2,.6,T,F)
mcCWplus(100,10,2,.6,T,F)
# with constant
mcCWplus(50,10,2,.2,F,F)
mcCWplus(100,10,2,.2,F,F)
mcCWplus(50,10,2,.6,F,F)
mcCWplus(100,10,2,.6,F,F)
# m=3
mcCWplus(50,10,3,.2,T,F)
mcCWplus(50,10,3,.6,T,F)
mcCWplus(100,10,3,.2,T,F)
mcCWplus(100,10,3,.6,T,F)
# with constant
mcCWplus(50,10,3,.2,F,F)
mcCWplus(50,10,3,.6,F,F)
mcCWplus(100,10,3,.2,F,F)
mcCWplus(100,10,3,.6,F,F)
# m=4
mcCWplus(50,20,4,.2,T,F)
mcCWplus(50,20,4,.6,T,F)
mcCWplus(100,20,4,.2,T,F)
mcCWplus(100,20,4,.6,T,F)
# with constant
mcCWplus(50,20,4,.2,F,F)
mcCWplus(50,20,4,.6,F,F)
mcCWplus(100,20,4,.2,F,F)
mcCWplus(100,20,4,.6,F,F)
# m=5
mcCWplus(50,20,5,.2,T,F)
mcCWplus(50,20,5,.6,T,F)
mcCWplus(100,20,5,.2,T,F)
mcCWplus(100,20,5,.6,T,F)
# with constant
mcCWplus(50,20,5,.2,F,F)
mcCWplus(50,20,5,.6,F,F)
mcCWplus(100,20,5,.2,F,F)
mcCWplus(100,20,5,.6,F,F)

# Evaluating and Generating Results
suppressWarnings(overallRep())
suppressWarnings(overallRepMulti())

cat("See Results folder for resulting tables")
