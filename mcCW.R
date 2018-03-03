# Tm=100;n=10;clsize=3; frho=.2; ind=2; noCons=F
#########################

library("stringr")
mcCW <- function(Tm,n,clsize,frho,noCons){
  
  dir.create("logs",F)
  if(!file.exists("logs/mcerror.log")){
    cat(c("Single","Tm","n","clsize","frho","noCons","noPois","ind","crit\n"),file = "logs/mcerror.log",append = T)
  }

  inDirCW  <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  outDirCW <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  dirname  <- inDirCW
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  tempdir  <- paste0(outDirCW,"tempres_",Tm,"-",n,"-",clsize,"-",frho,nocStr,"_CW/")  
  
  cat("Initializing variables\n")
  
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",frho,nocStr,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz     <- length(zz[[1]])
  list     <- zz[[2]]
  gammas   <- zz[[3]][1:n]
  
  #########################
  # write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  dir.create(tempdir,F)
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
    
    datind <- ((ind -1) %% 8) + 1
    write.table(z,file = paste0("CLUSTERB/datt",datind,".csv"),row.names = FALSE,col.names = FALSE)
    ########################################################

  cat("Analyzing")
  
  RtoGauss <- function(crit=){
    script <- paste0("cd CLUSTERB && wine ~/.wine/drive_c/gauss6.0/tgauss -o -b -e 'stopval=",crit,"; rown=",Tm,"; coln=",n,"; RUN MONTECARLOB",datind,".GSS'")
    tempCW <- system(script, intern=TRUE,wait=TRUE)
    if(!any(grepl("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC|ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempCW))){
      cat(c(T,Tm,n,clsize,frho,noCons,F,ind,crit,"\n"),file = "logs/mcerror.log",append = T)
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
      tempRC  <- tempRC[grep(clpattern,tempRC)]
      relList <- extractClubs(tempRC)
    }
      
    ############################## Pre-evaluation ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxCWabs <- length(intersect(maxlAbs,list)); excCWabs<- length(setdiff(maxlAbs,list))
    maxCWrel <- length(intersect(maxlRel,list)); excCWrel<- length(setdiff(maxlRel,list))
    repCW <- matrix(c(maxCWabs,excCWabs,maxCWrel,excCWrel),1,4)
    
    gmmlCW <- t(matrix(rep("",2*n),n))
    gmmlCW[1,][maxlAbs] <- "c";  gmmlCW[2,][maxlRel] <- "c"
    return(list(repCW,gmmlCW))
  }
    temp <- lapply(c(.01,.05,.1), RtoGauss) 
    repCWs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlCWs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlCWs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repCWs)  <- c("crit01","crit05","crit1")
    tempdir <- paste0(savedir,"Results_",n,"-",clsize,nocStr,"_CW/")
    tempres <- list(repCWs,gmmlCWs) 
    save(tempres, file = paste0(tempdir,"res",ind,".rda"))
    return(tempres)
    ############################## END REPORT ##############################
  }
  
  
  
  ############################# CONSOLIDATION ######################
  completed <- dir(tempdir)
  range_ind <- as.numeric(gsub("res|.rda","", completed))
  range_ind <- setdiff(1:lenz,range_ind)
  consConcs <- lapply(range_ind, function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating\n")
  
  consConcs <- lapply(dir(tempdir), function(d) get(load(paste0(tempdir,d))))
  reportCW  <- do.call(rbind,lapply(consConcs, function(c) c[[1]]))
  gmmlsCW   <- do.call(rbind,lapply(consConcs, function(c) c[[2]]))
  gmmlsCW   <- rbind(gammas,gmmlsCW)
  savedir <- outDirCW
  outdir  <- paste0(savedir,"Results_",n,"-",clsize,nocStr,"_CW/")
  filedir <- paste0(Tm,"-",frho,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  
  save(reportCW,file = paste0(outdir,"reportCW-",filedir))
  save(gmmlsCW,file = paste0(outdir,"gmmlCW-",filedir))
  unlink(tempdir,T)
}
