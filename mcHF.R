# Tm=100;n=10;clsize=3; frho=.2; ind=17; noCons=T
#########################
dir.create("logs",showWarnings = F)
if(!file.exists("logs/mcerror.log")){
  cat(c("Single","Tm","n","clsize","frho","noCons",
        "noPois","ind","crit\n"),file = "logs/mcerrorHF.log",append = T)
}
library("stringr")

mcHF <- function(Tm,n,clsize,frho,noCons){
  
  inDirHF  <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  outDirHF <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  dirname  <- inDirHF
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  
  cat("Initializing variables\n")
  
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",frho,nocStr,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz     <- length(zz[[1]])
  list     <- zz[[2]]
  gammas   <- zz[[3]][1:n]
  
  #########################
  # write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
    
    datind <- ((ind -1) %% 8) + 1
    write.table(z,file = paste0("CLUSTERHF/datt",datind,".csv"),row.names = FALSE,col.names = FALSE)
    ########################################################
    
    cat("Analyzing")
    
    RtoGauss <- function(crit){
      script <- paste0("cd CLUSTERHF && wine ~/.wine/drive_c/gauss6.0/tgauss -o -b -e 'stopval=",crit,"; rown=",Tm,"; coln=",n,"; RUN MONTECARLO",datind,".GSS'")
      tempHF <- system(script, intern=TRUE,wait=TRUE)
      if(!any(grepl("ASYMPTOTICALLY PERFECT CONVERGENCE FOR MC|ASYMPTOTICALLY RELATIVE CONVERGENCE FOR MC",tempHF))){
        cat(c(T,Tm,n,clsize,frho,noCons,F,ind,crit,"\n"),file = "logs/mcerror.log",append = T)
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
      
      ############################## Pre-evaluation ##############################
      
      indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
      maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
      maxlAbs   <- absList[[maxIndAbs]]
      indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
      maxIndRel <- sample(which(max(indRel)==indRel),1)
      maxlRel   <- relList[[maxIndRel]]
      
      maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(maxlAbs,list))
      maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(maxlRel,list))
      repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
      
      gmmlHF <- t(matrix(rep("",2*n),n))
      gmmlHF[1,][maxlAbs] <- "c";  gmmlHF[2,][maxlRel] <- "c"
      return(list(repHF,gmmlHF))
      
    }
    
    temp <- lapply(c(.01,.05,.1), RtoGauss) 
    repHFs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlHFs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlHFs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repHFs)  <- c("crit01","crit05","crit1")
    return(list(repHFs,gmmlHFs))
    
    
    ############################## END REPORT ##############################
  }
  
  
  
  ############################# CONSOLIDATION ######################
  begin <- c(1, lenz/2-2+1, lenz-4+1)
  end   <- c(lenz/2-2, lenz-4, lenz)
  consConcs1 <- mclapply(begin[1]:end[1],function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)},mc.cores = 8)
  cat("Finished: Part 1\n")
  consConcs2 <- mclapply(begin[2]:end[2],function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)},mc.cores = 8)
  cat("Finished: Part 2\n")
  consConcs3 <- mclapply(begin[3]:end[3],function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)},mc.cores = 4)
  cat("Finished: Part 3\n")
  consConcs <- c(consConcs1,consConcs2,consConcs3)
  rm(list = c("consConcs1","consConcs2","consConcs3"))
  cat("Consolidating\n")
  
  reportHF  <- do.call(rbind,lapply(consConcs, function(c) c[[1]]))
  gmmlsHF   <- do.call(rbind,lapply(consConcs, function(c) c[[2]]))
  gmmlsHF   <- rbind(gammas,gmmlsHF)
  
  
  savedir <- outDirHF
  outdir  <- paste0(savedir,"Results_",n,"-",clsize,nocStr,"_HF/")
  filedir <- paste0(Tm,"-",frho,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  
  save(reportHF,file = paste0(outdir,"reportHF-",filedir))
  save(gmmlsHF,file = paste0(outdir,"gmmlHF-",filedir))
  
}

