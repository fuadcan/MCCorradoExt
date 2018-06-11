# n <- 20; clsize=10;Tm=100;frho=.2;noCons=T
library("rugarch")

hitRat <- function(maxScc,exc,clsize,n){
  UU <- maxScc; DU <- clsize-maxScc ; UD <- exc; DD <- n-clsize-exc
  hitR   <- (UU+DD)/(UU+DD+UD+DU)
  H      <- (UU)/(UU+DU)
  F      <- (UD)/(UD+DD)
  KS     <- H-F
  return(list(hitR,H,F,KS))
}

PTestHF<- function(dat){

  gmms    <- dat[1,]; dat<-dat[-1,]
  indlstA <- seq(1,nrow(dat),2); lstA <- dat[indlstA,]
  indlstR <- seq(2,nrow(dat),2); lstR <- dat[indlstR,]
  
  
  actua  <- (gmms==1)*1-(gmms!=1)*1
  forecA <- matrix(1,nrow(lstA),ncol(lstA)); forecA[lstA!="c"]<- -1
  forecR <- matrix(1,nrow(lstR),ncol(lstR)); forecR[lstR!="c"]<- -1   
  
  actua  <- rep(actua,nrow(dat)/2)
  dacA   <- DACTest(c(t(forecA)),actua,"PT",.05)$Stat / sqrt(nrow(forecA))

  dacR   <- DACTest(c(t(forecR)),actua,"PT",.05)$Stat / sqrt(nrow(forecR))
  dac    <- c(dacA,dacR)
    return(dac)
}

PTestAGK<- function(dat){
  
  gmms    <- dat[1,]
  dat     <- dat[-1,]
  

  ind.01  <- grepl("adf-0.01",rownames(dat))
  ind.05  <- grepl("adf-0.05",rownames(dat))
  ind.1   <- grepl("adf-0.10",rownames(dat))
  
  lists <- list()
  lists <- list(dat[ind.01,],dat[ind.05,],dat[ind.1,])
  
  actua <- (gmms==1)*1-(gmms!=1)*1
  forec <- lapply(lists, function(l) matrix(1,nrow(l),ncol(l)))
  for(i in 1:length(lists)){forec[[i]][which(lists[[i]]=="")]<- -1}

  actuaVecs <- lapply(forec, function(f) rep(actua,nrow(f)))
  forecVecs <- lapply(1:3, function(i) c(t(forec[[i]])))
  
  
  ptRess <- sapply(1:3, function(x) DACTest(forecVecs[[x]],actuaVecs[[x]],"PT",.05)$Stat/sqrt(nrow(forec[[x]])))
  
  return(ptRess)
  
}


anlysofoutptAGK<- function(Tm,n,clsize,noCons){
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  dirName  <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  calcRep<-list()
  for(frho in c(2,6)){calcRep[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,"-AGK",nocStr,"/reportAGK-",Tm,"-",frho/10,".rda")))}
  
  gmml<-list()
  for(frho in c(2,6)){gmml[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,"-AGK",nocStr,"/gmmlAGK-",Tm,"-",frho/10,".rda")))}
  
  gmmlCW <- list()
  for(frho in c(2,6)){gmmlCW[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_CW/gmmlCW-",Tm,"-",frho/10,".rda")))}
  
  
  nT     <- sum(!duplicated(rownames(calcRep[[2]])))
  trials <- nrow(calcRep[[2]])/nT
  TTm    <- trials*nT
  
  anlysOutp <- function(calcRep,clsize,frho,gmml){
    
    ind <- frho*10
    
    errorind <- apply(gmmlCW[[ind]]=="",1,sum)==n
    errorind <- matrix(errorind[-1],,2,byrow = T)[,1]
    
    indadf.01  <- grepl("adf-0.01",rownames(calcRep[[ind]])) & !errorind 
    indadf.05  <- grepl("adf-0.05",rownames(calcRep[[ind]])) & !errorind 
    indadf.1   <- grepl("adf-0.1",rownames(calcRep[[ind]])) & !errorind 
    
    adf.01  <- c(apply(calcRep[[ind]][indadf.01,],2,sum),nrow(calcRep[[ind]][indadf.01,]))
    adf.05  <- c(apply(calcRep[[ind]][indadf.05,],2,sum),nrow(calcRep[[ind]][indadf.05,]))
    adf.1   <- c(apply(calcRep[[ind]][indadf.1,],2,sum),nrow(calcRep[[ind]][indadf.1,]))
    
    adf.01prf  <- sum(calcRep[[ind]][indadf.01,1]==clsize & calcRep[[ind]][indadf.01,2]==0)
    adf.05prf  <- sum(calcRep[[ind]][indadf.05,1]==clsize & calcRep[[ind]][indadf.05,2]==0)
    adf.1prf   <- sum(calcRep[[ind]][indadf.1,1]==clsize & calcRep[[ind]][indadf.1,2]==0)
    
    prfs  <- c(adf.01prf,adf.05prf,adf.1prf)
    repss <- cbind(adf.01,adf.05,adf.1)
    prfs  <- prfs / repss[3,]
    
    hitReps <- sapply(1:3, function(x) unlist(hitRat(repss[1,x],repss[2,x],clsize*repss[3,x],n*repss[3,x])))
    rownames(hitReps) <- c("HitR", "H", "F", "KS")
    
    anlysPT <- PTestAGK(gmml[[ind]][c(T,!errorind),])
    rep     <- rbind(hitReps,anlysPT,prfs)
    
    return(rep)
  }
  
  tempRep <- lapply(c(2,6)/10, function(frho) anlysOutp(calcRep,clsize,frho,gmml))
  reps    <- lapply(1:6, function(x) rbind(tempRep[[1]][x,],tempRep[[2]][x,]))  
  
  
  
  return(reps)
  
}

anlysofoutptHF <- function(Tm,n,clsize,noCons){
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  dirName  <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  calcRep<-list()
  for(frho in c(2,6)){calcRep[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_HF/reportHF-",Tm,"-",frho/10,".rda")))}
  
  gmml<-list()
  for(frho in c(2,6)){gmml[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_HF/gmmlHF-",Tm,"-",frho/10,".rda")))}
  
  gmmlCW<-list()
  for(frho in c(2,6)){gmmlCW[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_CW/gmmlCW-",Tm,"-",frho/10,".rda")))}
  
  nT     <- sum(!duplicated(rownames(calcRep[[2]])))
  trials <- nrow(calcRep[[2]])/nT
  TTm    <- trials*nT
  
  
  anlysOutp <- function(calcRep,clsize,frho,gmml){
    ind <- frho*10
    
    errorind <- apply(gmmlCW[[ind]]=="",1,sum)==n
    calcRepA <- calcRep[[ind]][!matrix(errorind[-1],,2,byrow = T)[,1],c(1,2)]  
    calcRepR <- calcRep[[ind]][!matrix(errorind[-1],,2,byrow = T)[,1],c(3,4)]
    
    HFA.01   <- calcRepA[rownames(calcRepA) == "crit01",]; HFR.01   <- calcRepR[rownames(calcRepR) == "crit01",]
    HFA.05   <- calcRepA[rownames(calcRepA) == "crit05",]; HFR.05   <- calcRepR[rownames(calcRepR) == "crit05",]
    HFA.1    <- calcRepA[rownames(calcRepA) == "crit1", ]; HFR.1    <- calcRepR[rownames(calcRepR) == "crit1", ]
    
    hitRatA  <- sapply(list(HFA.01,HFA.05,HFA.1), function(d) c(apply(d,2,sum),nrow(d)))
    hitRatR  <- sapply(list(HFR.01,HFR.05,HFR.1), function(d) c(apply(d,2,sum),nrow(d)))
    
    perfA    <- sapply(list(HFA.01,HFA.05,HFA.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    perfR    <- sapply(list(HFR.01,HFR.05,HFR.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    
    perfA    <- perfA / hitRatA[3,]
    perfR    <- perfR / hitRatR[3,]
    
    repA <- apply(hitRatA, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*hr[3],n*hr[3])))
    repR <- apply(hitRatR, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*hr[3],n*hr[3])))
    
    gmmls   <- lapply(c("01","05","1"), function(crit) gmml[[ind]][!errorind & (rownames(gmml[[ind]]) %in% c("gammas",paste0(c("abs","rel"), crit))), ])
    ptestAR <- sapply(gmmls, PTestHF)
    
    rownames(repA) <- rownames(repR) <- c("HitRat","H","F","KS")
    
    repA <- rbind(repA,ptestAR[1,],perfA)
    repR <- rbind(repR,ptestAR[2,],perfR)
    
    anlys <- cbind(repA,repR)
    colnames(anlys) <- unlist(lapply(c("A","R"), function(l) paste0(l,c(".01",".05",".1"))))
    
    return(anlys)
  }
  
  tempRepp <- lapply((c(2,6))/10, function(d) anlysOutp(calcRep,clsize,d,gmml))
  
  reps     <- lapply(1:6, function(x) rbind(tempRepp[[1]][x,],tempRepp[[2]][x,]))
  
  return(reps)
}

anlysofoutptCW <- function(Tm,n,clsize,noCons){
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  dirName  <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  calcRep<-list()
  for(frho in c(2,6)){calcRep[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_CW/reportCW-",Tm,"-",frho/10,".rda")))}
  
  gmml<-list()
  for(frho in c(2,6)){gmml[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_CW/gmmlCW-",Tm,"-",frho/10,".rda")))}
  
  nT     <- sum(!duplicated(rownames(calcRep[[2]])))
  trials <- nrow(calcRep[[2]])/nT
  TTm    <- trials*nT
  
  anlysOutp <- function(calcRep,clsize,frho,gmml){
    ind <- frho*10
    
    # errorind <- which(apply(gmml[grepl("abs",rownames(gmml)),]=="",1,sum)==n)
    errorind <- apply(gmml[[ind]]=="",1,sum)==n
    calcRepA <- calcRep[[ind]][!matrix(errorind[-1],,2,byrow = T)[,1],c(1,2)]  
    calcRepR <- calcRep[[ind]][!matrix(errorind[-1],,2,byrow = T)[,1],c(3,4)]
    
    
    CWA.01   <- calcRepA[rownames(calcRepA) == "crit01",]; CWR.01   <- calcRepR[rownames(calcRepR) == "crit01",]
    CWA.05   <- calcRepA[rownames(calcRepA) == "crit05",]; CWR.05   <- calcRepR[rownames(calcRepR) == "crit05",]
    CWA.1    <- calcRepA[rownames(calcRepA) == "crit1", ]; CWR.1    <- calcRepR[rownames(calcRepR) == "crit1", ]
    
    hitRatA  <- sapply(list(CWA.01,CWA.05,CWA.1), function(d) c(apply(d,2,sum),nrow(d)))
    hitRatR  <- sapply(list(CWR.01,CWR.05,CWR.1), function(d) c(apply(d,2,sum),nrow(d)))
    
    perfA    <- sapply(list(CWA.01,CWA.05,CWA.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    perfR    <- sapply(list(CWR.01,CWR.05,CWR.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    
    perfA    <- perfA / hitRatA[3,]
    perfR    <- perfR / hitRatR[3,]
    
    repA <- apply(hitRatA, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*hr[3],n*hr[3])))
    repR <- apply(hitRatR, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*hr[3],n*hr[3])))
    
    gmmls   <- lapply(c("01","05","1"), function(crit) gmml[[ind]][!errorind & (rownames(gmml[[ind]]) %in% c("gammas",paste0(c("abs","rel"), crit))), ])
    ptestAR <- sapply(gmmls, PTestHF)
    
    rownames(repA) <- rownames(repR) <- c("HitRat","H","F","KS")
    
    repA <- rbind(repA,ptestAR[1,],perfA)
    repR <- rbind(repR,ptestAR[2,],perfR)
    
    anlys <- cbind(repA,repR)
    colnames(anlys) <- unlist(lapply(c("A","R"), function(l) paste0(l,c(".01",".05",".1"))))
    
    return(anlys)
  }
  
  tempRepp <- lapply((c(2,6))/10, function(d) anlysOutp(calcRep,clsize,d,gmml))
  
  reps     <- lapply(1:6, function(x) rbind(tempRepp[[1]][x,],tempRepp[[2]][x,]))
  
  return(reps)
}

# TmVec <- c(50,75,100,200); nVec <- c(10,20,30,40); clsize <- c(3,5,7,10)
# TmVec <- c(50,75,100); nVec <- c(10,20,30,40); clsize <- c(3,5,7,10)
# TmVec <- c(50,100); nVec <- c(10,20)

overallAnlys <- function(TmVec,n,clsize,noCons) {
  nlyAGK   <- lapply(TmVec, function(Tm) anlysofoutptAGK(Tm,n,clsize,noCons))
  nlyAGK   <- lapply(1:length(nlyAGK[[1]]), function(x) do.call(rbind,lapply(nlyAGK, function(n) n[[x]])))
  nlyHF    <- lapply(TmVec, function(Tm) anlysofoutptHF(Tm,n,clsize,noCons))
  nlyHF    <- lapply(1:length(nlyHF[[1]]), function(x)  do.call(rbind,lapply(nlyHF, function(n) n[[x]])))
  nlyCW    <- lapply(TmVec, function(Tm) anlysofoutptCW(Tm,n,clsize,noCons))
  nlyCW    <- lapply(1:length(nlyCW[[1]]), function(x)  do.call(rbind,lapply(nlyCW, function(n) n[[x]])))
  
  nocstr <- if(noCons){"A"} else {"R"}
  ovrNly   <- lapply(1:length(nlyAGK), function(x) cbind(nlyAGK[[x]],nlyHF[[x]][,grepl(nocstr,colnames(nlyHF[[x]]))],nlyCW[[x]][,grepl("A",colnames(nlyCW[[x]]))]))
  ovrNly   <- lapply(1:length(nlyAGK), function(x) ovrNly[[x]][order(rep(c(.2,.6),length(TmVec))),])
  for(i in 1:length(ovrNly)){colnames(ovrNly[[i]]) <- sapply(c("MCL","HF","CW"), function(m) paste0(m, c(".01",".05",".1")))}
  
  return(ovrNly)
  
}

overallAnlys(c(50,100),10,3,T)
overallAnlys(c(50,100),10,3,F)
overallAnlys(c(50,100),10,5,T)
overallAnlys(c(50,100),10,5,F)
overallAnlys(c(50,100),20,3,T)
overallAnlys(c(50,100),20,3,F)
overallAnlys(c(50,100),20,5,T)
overallAnlys(c(50,100),20,5,F)
overallAnlys(c(50,100),20,7,T)
overallAnlys(c(50,100),20,7,F)
overallAnlys(c(50,100),20,10,T)
overallAnlys(c(50,100),20,10,F)
overallAnlys(c(50,100),30,3,T)
overallAnlys(c(50,100),30,3,F)
overallAnlys(c(50,100),30,5,T)
overallAnlys(c(50,100),30,5,F)
overallAnlys(c(50,100),30,7,T)
overallAnlys(c(50,100),30,7,F)
overallAnlys(c(50,100),30,10,T)
overallAnlys(c(50,100),30,10,F)


overall <- function(TmVec,nVec,clsize,noCons) {
  
  tempList <- lapply(nVec,  function(n) overallAnlys(TmVec,n,clsize,noCons))
  tempList <- lapply(1:length(tempList[[1]]), function(y) lapply(1:length(nVec), function(x) tempList[[x]][[y]]))
  tempList <- lapply(tempList, function(tlist) do.call(rbind,tlist))
  
  return(tempList)
}


overallRep<-function(TmVec,nVec){
  rep3T  <- overall(TmVec,nVec, 3,T)
  rep5T  <- overall(TmVec,nVec, 5,T)
  rep7T  <- overall(TmVec,nVec[!nVec==10], 7,T)
  rep10T <- overall(TmVec,nVec[!nVec==10], 10,T)
  rep3F  <- overall(TmVec,nVec, 3,F)
  rep5F  <- overall(TmVec,nVec, 5,F)
  rep7F  <- overall(TmVec,nVec[!nVec==10], 7,F)
  rep10F <- overall(TmVec,nVec[!nVec==10], 10,F)
  
  noConst   <- lapply(1:length(rep3T), function(x) round(rbind(rep3T[[x]],rep5T[[x]],rep7T[[x]],rep10T[[x]]),2))
  withConst <- lapply(1:length(rep3F), function(x) round(rbind(rep3F[[x]],rep5F[[x]],rep7F[[x]],rep10F[[x]]),2))
  
  dtType3  <- sapply(nVec, function(n) 
    sapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm) paste(3,n,frho,Tm,sep = "-"))))
  dtType5  <- sapply(nVec, function(n) 
    sapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm) paste(5,n,frho,Tm,sep = "-"))))
  dtType7  <- sapply(nVec[!nVec==10], function(n) 
    sapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm) paste(7,n,frho,Tm,sep = "-"))))
  dtType10 <- sapply(nVec[!nVec==10], function(n) 
    sapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm) paste(10,n,frho,Tm,sep = "-"))))
  dtType   <-c(c(dtType3),c(dtType5),c(dtType7),c(dtType10))
  
  for(i in 1:6){rownames(noConst[[i]]) <- dtType}; noConst   <- noConst[2:6]
  for(i in 1:6){rownames(withConst[[i]])<-dtType} ;withConst <- withConst[2:6]
  
  for(i in 1:5){colnames(noConst[[i]]) <- paste0(colnames(noConst[[i]]),"_",c("H","F","KS","PT","Perf")[i])}
  for(i in 1:5){colnames(withConst[[i]]) <- paste0(colnames(withConst[[i]]),"_",c("H","F","KS","PT","Perf")[i])}
  
  
  
  KSrepNoConst <- do.call(cbind,noConst[1:3]); KSrepWithConst <- do.call(cbind,withConst[1:3])
  PTSnoConst   <- cbind(noConst[[4]],rep(NA,nrow(noConst[[4]])),noConst[[5]]); 
  PTSwithConst <- cbind(withConst[[4]],rep(NA,nrow(withConst[[4]])),withConst[[5]])
  
  dir.create("Results",recursive = TRUE)
  write.csv(KSrepNoConst, file="Results/single_kupiersScoreNoCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(PTSnoConst, file="Results/single_PTandPerfSuccNoCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(KSrepWithConst, file="Results/single_kupiersScoreWithCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(PTSwithConst, file="Results/single_PTandPerfSuccWithCons.csv",sep = ";",row.names = T, col.names = T)  
  
#   return(list(KSrepNoConst,PTSnoConst,KSrepWithConst,PTSwithConst))
}

TmVec <- c(50,100); nVec <- c(10,20,30)
overallRep(TmVec,nVec)
