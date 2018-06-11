# n<-20;Tm=100;frho=0.6;k=4; noCons<-T;nopois <- F MISSING for CW
# n<-10;Tm=100;frho=0.6;k=2; noCons<-T;nopois <- F
nlyzA <- function(Tm,n,k,frho,noCons,nopois){
  
  nocStr   <- if(noCons){"noCons"} else {"withCons"}
  nopoistr <- if(nopois){NULL} else {"pois"}
  dir <- paste0("Output/",nocStr,"/multiClub/")
  filename <- paste0("Results_",n,"-",k,"-A-",nocStr,"Mlt/reportA-",Tm,"-",frho,nopoistr,".rda")
  
  rep <- get(load(paste0(dir,filename)))
  
  nlyz <- apply(rep,1,sum)
  names(nlyz) <- c(".01",".05",".1")
  return(nlyz)
  
}


nlyzHF <- function(Tm,n,k,frho,noCons,nopois){
  
  nocStr   <- if(noCons){"noCons"} else {"withCons"}
  nopoistr <- if(nopois){NULL} else {"pois"}
  dir   <- paste0("Output/",nocStr,"/multiClub/")
  filename <- paste0(dir,"Results_",n,"-",k,"-",nocStr,"Mlt_HF/reportHF-",Tm,"-",frho,nopoistr,".rda")
  repHF <- get(load(filename))[,(!noCons)+1]
  repHF <- cbind(repHF[names(repHF) == "crit01"],repHF[names(repHF) == "crit05"],repHF[names(repHF) == "crit1"])
  nlyz  <- apply(repHF,2,sum)
  names(nlyz) <- c(".01",".05",".1")
  return(nlyz)
}

nlyzCW <- function(Tm,n,k,frho,noCons,nopois){
  
  nocStr   <- if(noCons){"noCons"} else {"withCons"}
  nopoistr <- if(nopois){NULL} else {"pois"}
  dir   <- paste0("Output/",nocStr,"/multiClub/")
  filename <- paste0(dir,"Results_",n,"-",k,"-",nocStr,"Mlt_CW/reportCW-",Tm,"-",frho,nopoistr,".rda")
  repCW <- get(load(filename))[,(!T)+1] ### CHANGED (!noCons to T)
  repCW <- cbind(repCW[names(repCW) == "crit01"],repCW[names(repCW) == "crit05"],repCW[names(repCW) == "crit1"])
  nlyz  <- apply(repCW,2,sum)
  names(nlyz) <- c(".01",".05",".1")
  return(nlyz)
}

# TmVec <- c(50,100); n <- 30; k <- 5; noCons <- 
# n<-10;Tm=100;frho=0.6;k=2; noCons<-T;nopois <- F

overallnlyz <- function(TmVec,n,k,noCons){
  
  # errorMat  <- matrix(NA,length(TmVec),3); colnames(errorMat) <- c(".01",".05",".1")
  # errorList <- list(errorMat, errorMat)
  
  anlyzA  <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzA(Tm,n,k,frho,noCons,T), error = function(e) rep(NA,3)))))
  anlyzHF <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzHF(Tm,n,k,frho,noCons,T), error = function(e) rep(NA,3)))))
  anlyzCW <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzCW(Tm,n,k,frho,noCons,T), error = function(e) rep(NA,3)))))
  
  anlyzA_P  <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzA(Tm,n,k,frho,noCons,F), error = function(e) rep(NA,3)))))
  anlyzHF_P <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzHF(Tm,n,k,frho,noCons,F), error = function(e) rep(NA,3)))))
  anlyzCW_P <- lapply(c(.2,.6), function(frho) t(sapply(TmVec, function(Tm)  tryCatch(nlyzCW(Tm,n,k,frho,noCons,F), error = function(e) rep(NA,3)))))
  
  
  nlyz   <- cbind(do.call(rbind,anlyzA),do.call(rbind,anlyzHF),do.call(rbind,anlyzCW)) / numofrep
  nlyz_P <- cbind(do.call(rbind,anlyzA_P),do.call(rbind,anlyzHF_P),do.call(rbind,anlyzCW_P)) / numofrep
  
  rownames(nlyz)   <- c(paste0("nopois.2","-",TmVec),paste0("nopois.6","-",TmVec))
  rownames(nlyz_P) <- c(paste0("pois.2","-",TmVec),paste0("pois.6","-",TmVec))
  colnames(nlyz) <- colnames(nlyz_P) <- sapply(c("MCL","HF","CW"), function(n) c(paste0(n,c(".01",".05",".1"))))
  
  nlyz <- rbind(nlyz,nlyz_P)
  
  return(nlyz)  
  
}

overallRepMulti <- function(){
a1 <- lapply(c(T,F), function(nc) lapply(c(2,3), function(k) overallnlyz(c(50,100),10,k,nc)))
a2 <- lapply(c(T,F), function(nc) lapply(c(4,5), function(k) overallnlyz(c(50,100),20,k,nc)))
a3 <- lapply(c(T,F), function(nc) lapply(c(5,6), function(k) overallnlyz(c(50,100),30,k,nc)))
a4 <- lapply(c(T,F), function(nc) lapply(c(7,8), function(k) overallnlyz(c(50,100),40,k,nc)))

noConst   <- do.call(rbind,c(a1[[1]],a2[[1]],a3[[1]],a4[[1]]))
withConst <- do.call(rbind,c(a1[[2]],a2[[2]],a3[[2]],a4[[2]]))

rep <- cbind(noConst,withConst)
colnames(rep) <- c(sapply(c("MCL","HF","CW"), function(n) c(paste0(n,c("01_noInt","05_noInt",".1_noInt")))), sapply(c("MCL","HF","CW"), function(n) c(paste0(n,c("01_wInt","05_wInt",".1_wInt")))))

dtType  <- sapply(c("2-10","3-10","4-20","5-20","5-30","6-30","7-40","8-40"), function(kn)
  sapply(c("nopois","pois"), function(pois) sapply(c(".2",".6"), function(frho) 
    sapply(c(50,100), function(Tm) paste(kn,pois,frho,Tm,sep = "-")))))

rownames(rep) <- c(dtType)

dir.create("Results",recursive = TRUE)
write.csv(rep, file="Results/multi_perfectScores.csv", sep = ";",row.names = T, col.names = T)


}

# overallRepMulti()