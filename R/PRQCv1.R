###################################################### 
## Print-run QC
## Feb 26, 2004
##
## Two functions: PRv9mers
##              : PRvQCHyb
####################################################### 

## setwd("C:/MyDoc/Projects/SFGH/Quality/PrintRun/9mers/")
## source("C:/MyDoc/Projects/Rpackages/arrayQuality/R/PRQCv1.R")
## TOUSE: PRv9mers(fnames="9Mm.28.gpr", path="09Mm", prname="09Mm", DEBUG=TRUE)
## PRv9mers(path="08Hs", prname="08Hs", DEBUG=TRUE)

PRv9mers<-  function(fnames,
                    path=".",
                    dev = "png",  #set default to be png
                    DEBUG=FALSE,
                    prargs=NULL,
                    samepr=TRUE,
                    prname="xMm",
                    save = TRUE,
                    ...) 
{
  ## Setting defaults
  require(marray)
  require(limma)
  
  defs <- list(fill = TRUE, name.Gf="F532 Median", name.Gb="B532 Median", check.names=FALSE, as.is=TRUE,
               comment.char="", sep="\t", header=TRUE, quote = "\"", name.W="Flags")
  opt <- list(...)
  Args <- maDotsMatch(opt, defs)
  estpMatrix <- signalMatrix<- NULL
  
  if (DEBUG) print("Getting File Names")
  ## Getting File names
  if (missing(fnames))
    fnames <- file.path(path, dir(path, pattern = "*\\.gpr$"))
  

  if (is.null(path)) {
    fnames <- fnames
    resdir <- paste(prname, "PRQC", sep="")
  }else {
    resdir <- file.path(path, paste(prname,"PRQC", sep=""))
    #fnames <- file.path(path, fnames)
  }
  dir.create(resdir)
  print(resdir)
  
  for(f in fnames)
    {

      ## Set up output name
      if (DEBUG) print("Name the output file")
      tmp <- unlist(strsplit(f, "/"))
      fnopath <- tmp[length(tmp)]  ## filename without the path attached
      tmp <- unlist(strsplit(fnopath, "\\."))
      fstart <- paste(tmp[-length(tmp)], collapse=".")

      ## Set up arguments
      read.args <- maDotsMatch(Args, formals(args("read.GenePix")))
      read.args$fnames <- f
      #read.args$path <- read.args$name.Rf <- read.args$name.Rb <- NULL
       read.args$name.Rf <- read.args$name.Rb <- NULL
      read.args$path <- path
      if(DEBUG) print(read.args$path)
      read.args <- c(read.args, list(name.Rf=NULL, name.Rb=NULL))
      if(DEBUG) cat("Reading", read.args$file, "...\n")
      mraw <- do.call("read.GenePix", read.args) 
      GInfo <- maGeneTable(mraw)

      ###################
      ## EmAlgorithm split
      if(DEBUG) print("start EM")
      fg <- log(mraw@maGf, 2)
      cutoff <- quantile(log(mraw@maGb,2), prob=0.8, na.rm=TRUE)
      estp <- EMSplit(vect=fg, cutoff=cutoff)
      if(samepr){
        signalMatrix <- cbind(signalMatrix, fg)
        estpMatrix <- cbind(estpMatrix, estp)}
      else
        {
          tmp <- GInfo[,"ID"]
          write.table(cbind(tmp, estp), file=file.path(resdir, paste(fnopath, ".9mer.xls", sep="")),
                      row.names=FALSE, sep="\t")
        }
      if(DEBUG) print("Done...EM")

      
      ###################
      ## Setting up output device
      ###################
      if (DEBUG) print("Name of output device")
      plotdef <- switch(dev,
                        "bmp" = list(dev=list(width=1200, height=1000, bg="white"), suffix="bmp"),
                        "jpeg" = list(dev=list(quality=100, width=1200, height=1000, bg="white"), suffix="jpeg"),
                        "jpg" =  list(dev=list(quality=100, width=1200, height=1000, bg="white"), suffix="jpeg"),
                        "postscript" = list(dev=list(paper="special", width=12, height=10, bg="white"), suffix="ps"),
                        "png" =  list(dev=list(width=1200, height=1000, bg="white"), suffix="png"),
                        list(dev=list(width=1200, height=1000,bg="white"), suffix="png"),
                    )
      if(!is.element(dev, c("bmp", "jpeg","png","postscript","jpg")))
        print("Format error, format will be set to PNG")

      ## Names
      fname <- paste("P9mers", fstart,  plotdef$suffix, sep=".")
      plotdef <- c(plotdef, list(main=paste(fname, ": 9 mers QC")))
      
      ###################
      ## Plot
      if(DEBUG) print("start layout")
      if(save)  do.call(dev, maDotsDefaults(opt, c(list(filename=file.path(resdir, fname)), plotdef$dev)) ) 

      ## Layout 
      layout(matrix(c(7,1,2,7,3,3,7,4,4,7,5,6), 3, 4), height=c(0.5, 4, 4), width = c(12, 5, 2, 7))
        
      ## 1) Boxplot split by Plate
      if(DEBUG) print("start 1")
      par(mar=c(5, 4, 4, 2) + 0.1)
      boxplot(mraw, xvar="maPlate", yvar="maLG", ylab="Log Intensity", las=2)
      
      ## 2) Boxplot split by Print-tip
      if(DEBUG) print("start 2")
      boxplot(mraw, yvar="maLG", ylab="Log Intensity")
      

      ## 3,4) maGf
      if(DEBUG) print("start 3/4")
      qpImage(mraw, xvar="maLG", main="Spatial: Log Green")
      
      ## 5) Density Plot
      if(DEBUG) print("start 5")
      xrange <- range(log(as.numeric(c(mraw@maGf, mraw@maGb)),2), na.rm=TRUE, finite=TRUE)
      tmp <- log(mraw@maGf,2)  ## Make sure there is no inf 
      Gfden <- density(tmp[is.finite(tmp)]) ## doesn't handle inf values
      plot(Gfden, main="Foreground", xlim=xrange)
      if(!is.null(GInfo[,"ID"]))  ## make sure you there is a column name "ID"
        {
          yrange <- range(Gfden$y, na.rm=TRUE, finite=TRUE)
          tmpid <- maControls(mraw)
          tmp <- estp[tmpid == "probes"]
          text(xrange[2] - 3, yrange[2] - 0.1, paste("Probes Only", length(tmp)), cex=1.5,col=4)
          text(xrange[2] - 3, max(0, yrange[2] - 0.2), paste("Missing: ", sum(tmp <= 0.5)), cex=1.5, col=6)
        }
      
      ## 6) Density Plot Background
      if(DEBUG) print("start 6")
      tmp <- log(mraw@maGb,2)
      plot(density(tmp[is.finite(tmp)]), main="Background", xlim=xrange)
    
      ## 7) Title
      if(DEBUG) print("start 7")
      layout(1)
      par(mar=c(2,2,4,2))
      mtext(plotdef$main, line=3)
      mrawheader <- readGPRHeaders(f)
      mtext(paste("Date: ",  mrawheader$DateTime, " :: PMT", mrawheader$PMTGain), line=2, cex=0.8)
    
      
      if(DEBUG) print("Done...")
      ## Finishing
      if (save == TRUE) {
        cat(paste("save as", fname, "\n"))
        dev.off()
      }
    } ## end for

  if(samepr)
    {
      cat("Extracting missing probes \n")
      colnames(estpMatrix) <- fnames
      estpMatrixAve <- apply(estpMatrix, 1, median, na.rm=TRUE)
      signalMatrixAve <- apply(signalMatrix, 1, median, na.rm=TRUE)
      if(!is.null(GInfo[,"ID"]))
         {
           Info <- cbind(ID=GInfo[,"ID"], Name=GInfo[,"Name"],
                         estpMatrix, average=round(estpMatrixAve,3), Signal = signalMatrixAve)
           write.table(Info,file=file.path(resdir, paste(prname,"9mer.xls", sep="")), row.names=FALSE, sep="\t")
           tmpid <- maControls(mraw)
           write.table(Info[((tmpid == "probes") & (estpMatrixAve <=0.5)),],
                       file=file.path(resdir, paste(prname,"Missing.xls", sep="")), row.names=FALSE, sep="\t")
           write(as.vector(Info[,1]),
                 file=file.path(resdir, paste(prname, "QuickList.txt", sep="")), ncolumns=1)
         }
      cat("Done \n")
    }
}## end function


EMSplit <- function(vect, cutoff=9)
  {
    require(mclust)
    meV.na <- function(data, cutoff=9)
      {
        index <- (is.na(data) | !is.finite(data))
        y <- data[!index]
        testy <- me("V", y, cbind(y < cutoff, y > cutoff))
        p <- rep(NA, length(data))
        p[!index] <- testy$z[,2]
        return(p)
      }
    res <- meV.na(vect, cutoff=cutoff)
    return(res)
  }
