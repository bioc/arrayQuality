#####################################################
## Set of function written to estimate
## quality of microarray slides
##
## Author: Agnes Paquet
## Modified: 04/16/2004
##           09/21/2004
##           09/25/2004
## source("~/Projects/madman/Rpacks/arrayQuality/R/qualFunc.R")
#####################################################


###        TO DO  NEXT         ####
## 1) check that all call to HsRefenceDB are written properly
## 2) check boxplot if reference has Inf values
## 3) agQuality: pb when plotting, check isBelowLim and quality score
## 4) outputnormdata: allow choice of norm method


###################################################
## Given a gpr file
## Computes needed statistics to assess quality
###################################################

## Argument: result of readGPR
## Returns: matrix of numbers

slideQuality <- function(gprData=NULL, controlMatrix = controlCode, controlId = c("ID", "Name"), DEBUG=FALSE,...)
  {
    if (DEBUG) print("SlideQuality starting")
    
    # Check input argument
    if (is.null(gprData))
      stop("Slide information is missing:")
    
    # Read data in

    if (DEBUG) print("SlideQuality 1")

    controlId <- controlId[1]

    ############################################
    ## Channel specific data
    ############################################

    # Foreground and Background info

    Rf <- log.na(gprData[["RfMedian"]],2)
    Gf <- log.na(gprData[["GfMedian"]],2)
    Rb <- log.na(gprData[["RbMedian"]],2)
    Gb <- log.na(gprData[["GbMedian"]],2)

    rRf <- range(Rf, na.rm=TRUE)
    rangeRf <- rRf[2] - rRf[1]
    rGf <- range(Gf, na.rm=TRUE)
    rangeGf <- rGf[2] - rGf[1]
    rRb <- range(Rb, na.rm=TRUE)
    rangeRb <- rRb[2] - rRb[1]
    rGb <-range(Gb, na.rm=TRUE)
    rangeGb <- rGb[2] - rGb[1]

    RbMad <- mad(Rb, na.rm=TRUE)
    RbIqr <- IQR(Rb, na.rm=TRUE)
    GbMad <- mad(Gb, na.rm=TRUE)
    GbIqr <- IQR(Gb, na.rm=TRUE)
                                
    if (DEBUG) print("SlideQuality 3")
    
    # S2N ratios: summary of stat
    
    ifelse(length(gprData[["RfMean"]]) != 0,
           RS2N <- as.vector(log.na(gprData[["RfMean"]]/gprData[["RbMedian"]],2)),
           RS2N <- as.vector(log.na(gprData[["RfMean"]],2)))
    RS2Ninfo <- boxplot(RS2N, plot=FALSE)

    ifelse(length(gprData[["GfMean"]]) != 0,
           GS2N <- as.vector(log.na(gprData[["GfMean"]]/gprData[["GbMedian"]],2)),
           GS2N <- as.vector(log.na(gprData[["GfMean"]],2)))
    GS2Ninfo <- boxplot(GS2N, plot=FALSE)

    RS2Nmedian <- RS2Ninfo[[1]][3]
    GS2Nmedian <- GS2Ninfo[[1]][3]
    
    # Spots: summary of spot area
    spotArea <- median(gprData[["spotArea"]], na.rm=TRUE)
    spotRadius <- round(sqrt(spotArea) / pi)
    GBvar <- var(log.na(gprData[["GbMedian"]],2), na.rm=TRUE)
    RBvar <- var(log.na(gprData[["RbMedian"]],2), na.rm=TRUE)

    
    ############################################
    ## Combined Channel data
    ############################################    

    # M and A values (no background)
    if (DEBUG) print("SlideQuality 4")
    Mmean <- log.na(gprData[["RfMean"]],2) - log.na(gprData[["GfMean"]],2)
    Mmedian <- log.na(gprData[["RfMedian"]],2) - log.na(gprData[["GfMedian"]],2)

    Amean <- (log.na(gprData[["RfMean"]],2) + log.na(gprData[["GfMean"]],2))/2
    Amedian <- (log.na(gprData[["RfMedian"]],2) + log.na(gprData[["GfMedian"]],2))/2

    MmedInfo <- boxplot(Mmedian, plot=FALSE)
    AmedInfo <- boxplot(Amedian, plot=FALSE)

    
    # Comparison of mean(Mvalues) and median(Mvalues)
    # MMR amd mad

    goodM <- !(is.infinite(Mmean) | is.infinite(Mmedian) |
               is.na(Mmean) | is.na(Mmedian))
    
    MMR <- Mmean[goodM] - Mmedian[goodM]
    MMRmad <- mad(MMR, na.rm = TRUE)
    
    #number of spots with abs(mmr) > 0.5
    if (DEBUG) print("SlideQuality 5")
    numSpotOverMmrLim <- 0
    for(i in 1:length(MMR))
      {
        if (abs(MMR[i]) > 0.5)
          numSpotOverMmrLim <- numSpotOverMmrLim + 1
      }

    ifelse(length(MMR)!=0,
           percentSpotOverMmrLim <-  (numSpotOverMmrLim/length(MMR))*100,
           percentSpotOverMmrLim <- 0)
    
    #IQR(MMR)
    mmrIqr <- IQR(MMR, na.rm = TRUE)     
    
    # Normalization
    if (DEBUG) print("SlideQuality 6")

    goodAM <- !(is.infinite(Amedian) | is.infinite(Mmedian) |
                is.na(Amedian) | is.na(Mmedian))
    
    fit <- lowess(Amedian[goodAM], Mmedian[goodAM])

    MSE <- function(x, center) {
      if(missing(center))
        return(mean.na(x))
      else
        return(mean.na(sqrt(sum((x-center)^2))))
    }
    
    # MSE of lowess curve

    mseFit <- MSE(fit$y, center=0)

    # Print-tip information, mse by print-tip group
    if (DEBUG) print("SlideQuality 7")
    goodMed <-  !(is.infinite(Mmedian) | is.na(Mmedian))

    pt <- cbind(Block = as.vector(as.numeric(gprData[["Block"]])),
                M = Mmedian)
    #get the mean of Mvalues by print-tip group
    printTip <-  by(pt[,2],factor(pt[,1]), mean)
    goodPt <- !(is.infinite(printTip) | is.na(printTip))
    msePtip <- MSE(printTip[goodPt], center=0)


    if (DEBUG) print("SlideQuality 8")

    # Controls
    if (DEBUG) print("before arrayControls")
    Control <- arrayControls(gprData, SFGHControlcode = controlMatrix, id ="ID")
    numE <- numNeg <- numPos <- numProb <- 0
    emp <- neg <- pos <- NA

    for(i in 1:length(Control))
      {
        if (Control[i] == "Empty")
          {
            numE <- numE + 1
            emp <- c(emp, Amedian[i])
          }
        if (Control[i] == "Negative")
          {
            numNeg <- numNeg + 1
            neg <- c(neg, Amedian[i])
          }
        if (Control[i] == "Positive")
          {
            numPos <- numPos + 1
            pos <- c(pos, Amedian[i])
          }
      }

    EmptyMed <- median(as.numeric(emp), na.rm=TRUE)
    NegativeMed <- median(as.numeric(neg), na.rm=TRUE)
    PositiveMed <- median(as.numeric(pos), na.rm=TRUE)

    difEmptyNegative <- EmptyMed - NegativeMed
    difPositiveNegative <- PositiveMed - NegativeMed

    # Replicates
    if (DEBUG) print("SlideQuality 9")

    gId <- gprData[[controlId]]
    Replicates <- arrayReplicates(gprData, id = controlId)

    if (DEBUG) print(length(Replicates))
  
    index <- c()
    for(r in Replicates){
      for(i in 1:length(gId)){
        if (r == gId[i]) index <- c(index,i)
      }
    }

    if (length(index) == 0)
      repA <- NA
    else 
      repA <- Amedian[index]
    
    varRepA <- var(repA, na.rm=TRUE)

    # Flags
    Flags <- gprData[["Flags"]]
    numFlag <- Flags[Flags !=0]
    ifelse((length(Flags)!=0),
           percentFlag <- (length(numFlag)/length(Flags))*100,
           percentFlag <- 0)

    # Prepare for output
    # We want a list of all numbers, returned as a matrix
    # Easier to compare later

    sortedMeasures <- c("range RF", "range GF",
                        "- RB mad", "- GB mad",
                        "Median RS2N", "Median GS2N",
                        "- Median A for empty ctrl",
                        "- Median A for neg ctrl",
                        "Median A for positive ctrl",
                        "Pos ctl median A - Neg ctl median A",
                        "- Var replicated spots A values",
                        "- Mvalues MSE by print-tip",
                        "- MSE lowess",
                        "- % flagged spots",
                        "- Mvalues MMRmad",
                        "- % spots with Mvalues MMRmad>0.5",
                        )
        

    sortedRes <- c(rangeRf, rangeGf,
                   -RbMad, -GbMad,
                   RS2Nmedian, GS2Nmedian,
                   -EmptyMed, -NegativeMed, PositiveMed,
                   difPositiveNegative,
                   -varRepA,
                   -msePtip, -mseFit,
                   -percentFlag,                                      
                   -MMRmad, -percentSpotOverMmrLim,
                   )
    
    numResult <- as.matrix(sortedRes)
    rownames(numResult) <- sortedMeasures
    colnames(numResult) <- gprData[["File"]]
                
    if (DEBUG) print("SlideQuality done...")
    return(numResult)         
  }


##################################################
## scalRefTable: same as qualRefTable, but use
## reference table as is
## reference is a result from globalQuality
##################################################

scaleRefTable <- function(reference=NULL, organism=c("Mm", "Hs"))
  {
    
    organism <- organism[1]
    if(is.null(reference))
      {
        if(organism == "Mm")
          {
            if (!("MmReferenceDB" %in% ls(1)))
              data(MmReferenceDB)
            reference <- MmReferenceDB
          }
        else
          {
            if(!("HsReferenceDB" %in% ls(1)))
              data(HsReferenceDB)
            reference <- data(HsReferenceDB)
          }
      }

    ave <- apply(reference, 1, median, na.rm=TRUE)
    iqr <- apply(reference, 1, IQR, na.rm=TRUE)
    reftab <- cbind(median=ave,
                    iqr=iqr)
    return(reftab)
  }

## Takes a matrix/vector of numbers as argument
## Scales this matrix according to median and iqr of reference values
## if reference is missing, use reference values (Mm by default, or Hs)
## (globalQuality result)

arrayScal <- function(numMat, reference=NULL, organism=c("Mm", "Hs"))
{
  organism=organism[1]

  if(missing(numMat))
    stop("Input error, matrix to scale missing")
  else
    {
      if(is.null(reference))
        {
          if(organism == "Mm")
            {
              if (!("MmReferenceDB" %in% ls(1)))
                data(MmReferenceDB)
              reference <- MmReferenceDB
            }
          else
            {
              if(!("HsReferenceDB" %in% ls(1)))
                data(HsReferenceDB)
              reference <- data(HsReferenceDB)
            }
        }

      reftab <- scaleRefTable(reference=reference, organism=organism)
      
      ## Determine dimensions
      if(!is.null(dim(numMat)))
        {
          ncol <- dim(numMat)[2]
          nrow <- dim(numMat)[1]
        }
      else
        {
          ncol <- 1
          nrow <- length(numMat) ## numMat is a vector
        }
 
      #Check dimensions 
      if(nrow!=dim(reftab)[1])
        {
          stop("Error is scaling, files of different length")
        }
      else
        { 
          tmp <- numMat

          tmp2 <- sweep(tmp,1,reftab[,"median"], FUN="-")
          ##Pb: deal with iqr==0!!
          ##Test if Inf values don't crash boxplot
          tmp <- sweep(tmp2,1,reftab[,"iqr"],FUN="/")
           
          #for(j in 1:ncol)
          #  for(i in 1:nrow)
          #  {
          #    if(reftab[i,"iqr"]!=0)
          #      tmp[i,j] <- (as.numeric(tmp[i,j]) - as.numeric(reftab[i,"median"]))/
          # as.numeric(reftab[i,"iqr"])
          #    else
          #      {
          #        print("One of the ranges = 0")
          #        tmp[i,j] <- as.numeric(tmp[i,j]) - as.numeric(reftab[i,"median"])
          #      }              
          #  }
          
          return(tmp)
        }
    }
}


################################################
## HTML report
# Read in names of plots
# suppose directory contains
# same number of qualPlot and diagPlot
# based on alphabetical order

# nbBelow is a matrix: for each slide: number of measure below range of good slides and total number of measures

quality2HTML <- function(fnames=NULL, path=".", DiagPlot=NULL, QCplot=NULL, resdir=".", nbBelow=NULL)
  {
    print("starting HTML")
    HTwrap <- function(x, tag = "TD", option="align", value="center") {
      if (option == "")
        {
          paste("<", tag, "\"",  ">", x, "</", tag, ">", sep = "")
        }
      else
        paste("<", tag," ", option, "=\"",value, "\"",  ">", x, "</", tag, ">", sep = "")
    }
    
    HTimg <- function(src){
      paste("<a href= \"",src, "\">",
            "<img width=\"300\", height=\"200\", aligh=\"center\", src= \"",
            src, " \"/></a>", sep="")
    }

    HTscore <- function(filename,src){
      paste(filename, "<br>",
            "Measures below range: ", src)
    }

    
    tableTag <- function(fnames,fullDiagp,fullQCp) 
      {
        
        tab <- ""
        tr <- ""
        
        for(i in 1:min(length(fullQCp), length(fullBoxp)))
          {
            td2 <- HTwrap(HTimg(fullDiagp[i]), tag="TD")
            td1 <- HTwrap(HTimg(fullQCp[i]), tag="TD")
            td3 <- HTwrap(HTscore(fnames[i],paste(nbBelow[i,1], nbBelow[i,2], sep="/")), tag="TD",
                          option="align", value="left") 
            
            tr <- HTwrap(paste(td1, td2, td3, sep="\n"), tag="TR")
            tab <- paste(tab,tr,sep="")
          }
        return(tab)
        
      }

    if(missing(QCplot) || is.null(QCplot))
      {
        QCp <- dir(path, pattern="QCPlot*")
        fullQCp <- file.path(path, QCp)
      }
    else {
      fullQCp <- QCplot
    }

    if(missing(DiagPlot) || is.null(DiagPlot))
      {
        boxp <- dir(path, pattern="qualPlot*")
        fullBoxp <- file.path(path,boxp)
      }
    else {
      fullBoxp <- DiagPlot
    }

    outputfile <- file(file.path(resdir,"qualityReport.html"),"w") 
    datadir <- system.file("gprQCData", package="arrayQuality")
    html <- paste(readLines(file.path(datadir, "index.html")), "\n", collapse="")

    split <- unlist(strsplit(html, split="<a name=\"table\"></a>"))

    tab <- tableTag(fnames,fullQCp, fullBoxp)
    
    cat(split[1], tab,split[2], file=outputfile)     
    close(outputfile)
    print("End HTML")
  }

#slidequality is the result of slideQuality for ONE slide
qualityScore <- function(slidequality, organism=c("Mm", "Hs"), reference=NULL)
  {
    organism <- organism[1]
    slidequality <- as.vector(slidequality)

    if(is.null(reference))
      {
        if(organism == "Mm")
          {
            if (!("MmReferenceDB" %in% ls(1)))
              data(MmReferenceDB)
            reference <- MmReferenceDB
          }
        
        else
          {
            if(!("HsReferenceDB" %in% ls(1)))
              data(HsReferenceDB)
            reference <- data(HsReferenceDB)
          }
      }

    score <- matrix(0,nrow=length(slidequality), ncol=1)

    for(i in 1:length(slidequality))
      {
        vect <- reference[i,]
        score[i] <- (length(vect[vect < slidequality[i]])/length(vect))*100
      }
    return(score)

  }

##############################################
## Jean: Sep 21, 2004 : make it more general by allowing different input source
##Reads in gpr files
# returns matrix of QC measures
globalQuality <- function(fnames = NULL, path = ".",
                          organism=c("Mm", "Hs"),
                          output=FALSE,
                          resdir=".",
                          DEBUG = FALSE,
                          inputsource = "readGPR",
                          ...)
  {
    # Check input arguments
    if (DEBUG) print("Starting globalQuality")
    
    if (missing(fnames) | is.null(fnames))
      {
        if(inputsource == "readGPR") fnames <- dir(path, pattern = ".*\\.gpr$")
        if(inputsource == "readAgilent") fnames <- dir(path, pattern = ".*\\.txt$")
        if(inputsource == "readSpot") fnames <- dir(path, pattern = ".*\\.spot$")
      }
    
    organism <- organism[1]
        
    # Prepares results
    quality <- NULL
    
    # Call to slideQuality for each gpr file

    for (i in 1:length(fnames))
      {
        if (DEBUG) print("In the loop ")
        f <- fnames[i]
        gp <- do.call(inputsource, args=list(fnames=f, path=path))
        restmp <- slideQuality(gp)
        quality <- cbind(quality, restmp[,1])
        meas <- rownames(restmp)
      }
    
    colnames(quality) <- fnames
    rownames(quality) <- meas
       
    # Results
    if (output)
      write.table(quality, "quality.txt",sep="\t", col.names=NA)
    
    return(quality)
  }
  
#################################
## Plot function

# 1. Plot the boxplot
# 2. SuperImpose value for slides of interest (1 or more)
# arrayQuality and reference are results from
# slideQuality and globalQuality respectively


qualBoxplot <- function(arrayQuality=NULL,  reference=NULL, organism=c("Mm", "Hs"), DEBUG=FALSE,...)
  {
    if (is.null(arrayQuality) || missing(arrayQuality))
      stop("No data to plot")
    
    # Reference = output of globalQuality for ref slides
    # if NULL, reads in matrix store as RData
    organism<-organism[1]
    
    score <- qualityScore(arrayQuality)
    if(is.null(reference))
      {
        if(organism == "Mm")
          {
            if (!("MmReferenceDB" %in% ls(1)))
              data(MmReferenceDB)
            reference <- MmReferenceDB
          }
        
        else
          {
            if (!("HsReferenceDB" %in% ls(1)))
              data(HsReferenceDB)
            reference <- HsReferenceDB            
        }
      }

    scalref <- arrayScal(reference, reference=reference)    
    scalarray <- as.matrix(arrayScal(arrayQuality, reference=reference))
    
    if(is.null(dim(scalref))|nrow(arrayQuality)!=nrow(reference))
      stop("Input must be a matrix resulting from slideQuality.R")
    else
      {
        #Boxplot of reference arrays quality measure

        nr <- nrow(scalref)

        ##goodLim <- matrix(0,nrow=nr,ncol=1)
        ##pb, what if arrayQuality = more than 1 col!
        goodLim <- matrix(0, nrow=nr, ncol=ncol(arrayQuality))

        lim <- range(as.numeric(scalarray),as.numeric(scalref),na.rm=TRUE, is.finite=TRUE)
        plot(0:(nr+1),0:(nr+1),xlim=lim,type="n", axes=FALSE,ylab="",xlab="", xaxt="n",
             main="Array Quality Control Comparison")

        axis(2, at=1:nr, labels=rownames(arrayQuality), cex.axis=0.7, las=2)
        tmp <- c(paste(as.character(round(score)), "(", as.character(round(arrayQuality,1)), ")"), "% (value)")
        axis(4, at=1:(nr+1), labels=, tmp, cex.axis=0.8, las=2)
        axis(1, at=c(lim[1]+0.5, lim[2]-0.5), labels= c("Problematic", "Good"))

        ##No pb with Na's
        for(i in 1:nr)
          {
            if(!is.na(min(as.numeric(scalref[i,]))))  ## new line Sept, 2004
              {
                bp <- boxplot(as.numeric(scalref[i,]), at=i,
                              add=TRUE, horizontal=TRUE, axes=FALSE)
                goodLim[i,1] <- bp$stats[1]  ## pb if more than 1 col
              }
            else
              {
              goodLim[i, 1] <- NA
              }
          
            ## not ploting if NA
            tmp <- round(quantile(as.numeric(reference[i,]),
                                  probs=c(0.75, 0.25), na.rm=TRUE),1)
            text(quantile(as.numeric(scalref[i,]), probs=c(0.8, 0.2), na.rm=TRUE),rep((i+0.3),2),
                 as.character(tmp), cex=0.7, col="blue")
          }

        
        #Line for tested arrays
        #Number of measure below criteria
        nc <- ncol(scalarray)
        col <- rainbow(nc)

        isBelowLim <- matrix(FALSE, ncol=nc, nrow=nr)
        
        for(i in 1:nc)
          {
            lines(scalarray[,i], 1:nr, col=col[i], lty=2)
            points(scalarray[,i],1:nr, col=col[i], pch=15, cex=1.5)
            missingdata <- as.vector(is.na(scalarray[,i]))
            if(length(missingdata[missingdata]) > 0)
              {
                text(rep(lim[1]+ 1,length(missingdata[missingdata])),
                     c(1:nr)[missingdata],
                     rep("NA", length(missingdata[missingdata])),
                     col="green", cex=0.8)
              }
            
            for(j in 1:nr)
              {
                if(!is.na(goodLim[j,1]))
                  {
                    if(as.numeric(scalarray[j,i]) < goodLim[j,1] || is.na(scalarray[j,i]))
                      isBelowLim[j,i] <- TRUE
                  }
                else
                  {
                    isBelowLim[j,i] <- NA
                  }
              }
          }
        box()
        #legend
        leg.txt <- colnames(arrayQuality)
        legend(lim[1],nr+1,leg.txt,lty=1,col=col, cex=0.7)

        res <- character(0)
        for(i in 1:nc)
          res <- c(res, length(isBelowLim[isBelowLim[,i],i]))
        return(cbind(res, nr))
      }
  }

########################################################
## Write normalized data to .txt file
## By default, background subtraction IS performed
########################################################

outputNormData <- function(mraw=NULL, DEBUG = FALSE,...)
  {
    opt <- list(...)   
    norm.defs <- maDotsMatch(maDotsDefaults(opt, list(norm="p")),formals(args(maNorm)))

    if (DEBUG) cat("Using normalization method:  ", norm.defs$norm, "\n")
    mnorm <- do.call("maNorm", c(list(mraw), norm.defs))

    #mnorm <- maNorm(mraw, )
    write.marray(mnorm, "NormalizedData.xls")
  }


readcontrolCode <- function(file = "SpotTypes.txt", path = NULL, sep = "\t", check.names = FALSE, controlId=c("ID", "Name"), ...) 
  {
    require(limma)
    controlId <- controlId[1]
    spotTypes <- readSpotTypes(file=file, path=path, sep=sep, check.names=check.names, ...)
    controlCode <- spotTypes[, c(grep(controlId, colnames(spotTypes)), 1)]
    colnames(controlCode) <- c("Pattern", "Name")
    return(controlCode)
  }

