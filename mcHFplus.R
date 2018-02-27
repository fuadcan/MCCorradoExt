# Tm=50;n=10;k=2; frho=.2; ind=1; noCons =T; nopois=F
library("stringr")

mcHFplus <- function(Tm,n,k,frho,noCons,nopois){
  
  inDirHFm  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  outDirHFm <- if(noCons){"Output/noCons/multiClub/"} else {"Output/withCons/multiClub/"}
  
  cat("initializing parameters\n")
  
  dirname  <- inDirHFm
  
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
  
  
  # write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  
  
  repForDat <- function(ind){
    
  z      <- zz[[1]][[ind]]
    
    
    
   datind <- ((ind -1) %% 8) + 1
    write.table(z,file = paste0("CLUSTERHF/datt",datind,".csv"),row.names = FALSE,col.names = FALSE)
 
    cat("Analyzing\n")
    
    RtoGauss <- function(crit){
    ########################################################
    script <- paste0("cd CLUSTERHF && wine ~/.wine/drive_c/gauss6.0/tgauss -o -b -e 'stopval=",crit,"; rown=",Tm,"; coln=",n,"; RUN MONTECARLO",datind,".GSS'")
      tempHF <- system(script, intern=TRUE,wait=TRUE)
      if(!any(grepl("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC|ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempHF))){
        cat(c(F,Tm,n,clsize,frho,noCons,F,ind,crit,"\n"),file = "logs/mcerrorHF.log",append = T)
        absList <- list(rep(NA,3),rep(NA,3))
        relList <- list(rep(NA,3),rep(NA,3))
      } else {
        tempHF <- tempHF[(grep("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC",tempHF)):length(tempHF)]
        tempPC <- tempHF[(grep("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC",tempHF)):(grep("ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempHF)-1)]
        tempRC <- tempHF[(grep("ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempHF)):length(tempHF)]
        extractClubs <- function(temp){sapply(strsplit(temp,"\\s|c"), function(cl) {cl <- as.numeric(cl); cl[!is.na(cl)]})}
        clpattern <- "\\s+c[0-9]\\s+"
        tempPC  <- tempPC[grep(clpattern,tempPC)]
        absList <- extractClubs(tempPC)
        tempRC  <- tempRC[grep(clpattern,tempRC)]
        relList <- extractClubs(tempRC)
      }
    
      ############################## Evaluation ##############################
      # cat("Pre-Evaluation\n")
      listCode <- unlist(lapply(list, function(x) c(x,0)))
      
      sccHFabs05 <- suppressWarnings(mean(unlist(lapply(absList, function(x) c(sort(x),0)))==listCode)==1)
      sccHFrel05 <- suppressWarnings(mean(unlist(lapply(relList, function(x) c(sort(x),0)))==listCode)==1)
      
      repHF  <- matrix(c(sccHFabs05,sccHFrel05),1)
      gmmlHF <- t(matrix(rep("",2*n),n))
      for(i in 1:length(absList)){gmmlHF[1,][absList[[i]]] <- paste0("c",i)} 
      for(i in 1:length(relList)){gmmlHF[2,][relList[[i]]] <- paste0("c",i)}
      
    
    return(list(repHF,gmmlHF))
    ############################## END REPORT ##############################
  }
    
    temp    <- lapply(c("01","05","1"), RtoGauss) 
    repHFs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlHFs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlHFs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repHFs)  <- c("crit01","crit05","crit1")
    return(list(repHFs,gmmlHFs))
  }
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating Results\n")
  
  reportHF  <- do.call(rbind,lapply(consConcs, function(c) c[[1]]))
  gmmlsHF   <- do.call(rbind,lapply(consConcs, function(c) c[[2]]))
  
  
  savedir <- outDirHFm
  outdir  <- paste0(savedir,"Results_",n,"-",k,nocStr,"_HF/")
  filedir <- paste0(Tm,"-",frho,poistr,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  save(reportHF,file = paste0(outdir,"reportHF-",filedir))
  save(gmmlsHF,file = paste0(outdir,"gmmlHF-",filedir))
  
  cat("Finished\n")
  
}
