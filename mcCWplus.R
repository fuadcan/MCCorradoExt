# Tm=50;n=10;k=2; frho=.2; ind=1; noCons =T; nopois=F
library("stringr")

mcCWplus <- function(Tm,n,k,frho,noCons,nopois){
  
  inDirCWm  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  outDirCWm <- if(noCons){"Output/noCons/multiClub/"} else {"Output/withCons/multiClub/"}
  
  cat("initializing parameters\n")
  
  dirname  <- inDirCWm
  
  nocStr   <- if(noCons){"-noConsMlt"} else {"-withConsMlt"}
  poistr  <- if(nopois){NULL} else {"pois"}
  
  
  filename <- paste0("zZ_",Tm,"-",n,"-",k,"-",frho,nocStr,poistr,".rda")
  
  zz       <- get(load(paste0(dirname,filename)))
  lenz     <-length(zz[[1]])
  list     <- zz[[2]]
  gammas   <- zz[[3]]
  
  list <- list[sapply(list, length)!=1]
  fsts<- sapply(list, function(x) x[1])  
  list<- list[order(fsts)]
  
  
  write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  
  
  repForDat <- function(ind){
    
  z      <- zz[[1]][[ind]]
    
    
    
    write.table(z,file = "hfcodes/datt.csv",row.names = FALSE,col.names = FALSE)
    
    cat("Analyzing\n")
    
    RtoGauss <- function(crit="05"){
    ########################################################
      tempCW <- shell(paste0("C:/gauss6.0/tgauss -b ",'RUN hfcodes\\main',crit,'.gss'), intern=TRUE,wait=TRUE)
      tempCW <- tempCW[1:(which(grepl("GAUSS",tempCW))[1]-2)]
      
      aCrude<- strsplit(tempCW[1:((which(tempCW=="brkpnt"))-1)]," ")
      aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
      aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
      absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
      
      rCrude<- strsplit(tempCW[((which(tempCW=="brkpnt"))+1):length(tempCW)]," ")
      rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
      rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
      relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
      
      fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
      absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
      
      if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
      if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
      
    
      ############################## Evaluation ##############################
      # cat("Pre-Evaluation\n")
      listCode <- unlist(lapply(list, function(x) c(x,0)))
      
      sccCWabs05<- suppressWarnings(mean(unlist(lapply(absList, function(x) c(sort(x),0)))==listCode)==1)
      sccCWrel05<- suppressWarnings(mean(unlist(lapply(relList, function(x) c(sort(x),0)))==listCode)==1)
      
      repCW <- matrix(c(sccCWabs05,sccCWrel05),1)
      gmmlCW <- t(matrix(rep("",2*n),n))
      for(i in 1:length(absList)){gmmlCW[1,][absList[[i]]] <- paste0("c",i)} 
      for(i in 1:length(relList)){gmmlCW[2,][relList[[i]]] <- paste0("c",i)}
      
    
    return(list(repCW,gmmlCW))
    ############################## END REPORT ##############################
  }
    
    temp    <- lapply(c("01","05","1"), RtoGauss) 
    repCWs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlCWs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlCWs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repCWs)  <- c("crit01","crit05","crit1")
    return(list(repCWs,gmmlCWs))
  }
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating Results\n")
  
  reportCW  <- do.call(rbind,lapply(consConcs, function(c) c[[1]]))
  gmmlsCW   <- do.call(rbind,lapply(consConcs, function(c) c[[2]]))
  
  
  savedir <- outDirCWm
  outdir  <- paste0(savedir,"Results_",n,"-",k,nocStr,"_CW/")
  filedir <- paste0(Tm,"-",frho,poistr,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  save(reportCW,file = paste0(outdir,"reportCW-",filedir))
  save(gmmlsCW,file = paste0(outdir,"gmmlCW-",filedir))
  
  cat("Finished\n")
  
}
