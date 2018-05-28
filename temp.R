# Tm=50;n=10;k=2; frho=.2; ind=5; noCons =T; nopois=F
library("stringr")

mcCWplus <- function(Tm,n,k,frho,noCons,nopois){
  
  inDirCWm  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  outDirCWm <- if(noCons){"Output/noCons/multiClub/"} else {"Output/withCons/multiClub/"}
  
  cat("initializing parameters\n")
  
  dirname  <- inDirCWm
  
  nocStr   <- if(noCons){"-noConsMlt"} else {"-withConsMlt"}
  poistr  <- if(nopois){NULL} else {"pois"}
  tempdir  <- paste0(outDirCWm,"tempres_",Tm,"-",n,"-",k,"-",frho,nocStr,poistr,"_CW/")  
  
  filename <- paste0("zZ_",Tm,"-",n,"-",k,"-",frho,nocStr,poistr,".rda")
  
  zz       <- get(load(paste0(dirname,filename)))
  lenz     <-length(zz[[1]])
  list     <- zz[[2]]
  gammas   <- zz[[3]]
  
  list <- list[sapply(list, length)!=1]
  fsts <- sapply(list, function(x) x[1])  
  list <- list[order(fsts)]
  
  
  # write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  senderCW <- if(Sys.info()[1] == "Linux"){function(crit,Tm,n,datind){
    script <- paste0("cd CLUSTERB && wine ~/.wine/drive_c/gauss6.0/tgauss -o -b -e 'stopval=",crit,"; rown=",Tm,"; coln=",n,"; RUN MONTECARLOB",datind,".GSS'")
    system(script,intern=TRUE,wait=TRUE)
  }} else {function(crit,Tm,n,datind){ # check gauss version; and sys.name; and method (HF or CW)
    script <- paste0("cd C:\\Users\\", Sys.info()[7] ,"\\Documents\\MCCorradoExt\\CLUSTERB && C:\\gauss10\\tgauss -o -b -e stopval=",crit,";rown=",Tm,";coln=",n,";RUN MONTECARLOB",datind,".GSS")
    shell(script,intern=TRUE,wait=TRUE)}}  
  
  dir.create(tempdir,F)
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
    datind <- ((ind -1) %% 8) + 1
    write.table(z,file = paste0("CLUSTERB/datt",datind,".csv"),row.names = FALSE,col.names = FALSE)
    
    cat("Analyzing\n")
    
    RtoGauss <- function(crit){
      ########################################################
      
      tempCW  <- senderCW(crit,Tm,n,datind)
      if(!any(grepl("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC|ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempCW))){
        cat(c(F,Tm,n,k,frho,noCons,F,ind,crit,"\n"),file = "logs/mcerrorCW.log",append = T)
        absList <- list(rep(NA,3),rep(NA,3))
        relList <- list(rep(NA,3),rep(NA,3))
      } else {
        tempCW <- tempCW[(grep("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC",tempCW)):length(tempCW)]
        tempPC <- tempCW[(grep("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC",tempCW)):(grep("ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempCW)-1)]
        tempRC <- tempCW[(grep("ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempCW)):length(tempCW)]
        extractClubs <- function(temp){lapply(strsplit(temp,"\\s|c"), function(cl) {cl <- as.numeric(cl); cl[!is.na(cl)]})}
        clpattern <- "\\s+c[0-9]\\s+"
        tempPC  <- tempPC[grep(clpattern,tempPC)]
        absList <- extractClubs(tempPC)
        absList <- absList[sapply(absList, length)!=1]
        fsts    <- sapply(absList, function(x) x[1])  
        absList <- absList[order(fsts)]
        
        tempRC  <- tempRC[grep(clpattern,tempRC)]
        relList <- extractClubs(tempRC)
        relList <- relList[sapply(relList, length)!=1]
        fsts    <- sapply(relList, function(x) x[1])  
        relList <- relList[order(fsts)]
      }
      
      ############################## Evaluation ##############################
      # cat("Pre-Evaluation\n")
      listCode <- unlist(lapply(list, function(x) c(x,0)))
      
      sccCWabs05 <- suppressWarnings(mean(unlist(lapply(absList, function(x) c(sort(x),0)))==listCode)==1)
      sccCWrel05 <- suppressWarnings(mean(unlist(lapply(relList, function(x) c(sort(x),0)))==listCode)==1)
      
      repCW  <- matrix(c(sccCWabs05,sccCWrel05),1)
      gmmlCW <- t(matrix(rep("",2*n),n))
      for(i in 1:length(absList)){gmmlCW[1,][absList[[i]]] <- paste0("c",i)} 
      for(i in 1:length(relList)){gmmlCW[2,][relList[[i]]] <- paste0("c",i)}
      
      
      return(list(repCW,gmmlCW))
      ############################## END REPORT ##############################
    }
    
    temp    <- lapply(c(.01,.05,.1), RtoGauss) 
    repCWs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlCWs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlCWs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repCWs)  <- c("crit01","crit05","crit1")
    # tempdir <- paste0(savedir,"Results_",n,"-",k,nocStr,"_CW/")
    tempres <- list(repCWs,gmmlCWs) 
    save(tempres, file = paste0(tempdir,"res",ind,".rda"))
    return(tempres)
  }
  cat(paste0(outDirCWm,"Results_",n,"-",k,nocStr,"_CW/reportCW-",paste0(Tm,"-",frho,poistr,".rda")),"\n")
  if(!file.exists(paste0(outDirCWm,"Results_",n,"-",k,nocStr,"_CW/reportCW-",paste0(Tm,"-",frho,poistr,".rda")))){
  completed <- dir(tempdir)
  range_ind <- as.numeric(gsub("res|.rda","", completed))
  range_ind <- setdiff(1:lenz,range_ind)
  # consConcs <- mclapply(range_ind,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)},mc.cores=4)
  consConcs <- lapply(range_ind,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
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
  # return(list(reportCW,gmmlsCW))
  }
  unlink(tempdir,T)
}
mcCWplus(100,20,4,.6,F,F)
