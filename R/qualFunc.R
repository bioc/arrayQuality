######################################
## Set of function written to estimate
## quality of microarray slides
##
## Author: Agnes Paquet
## Modified: 02/05/2004


###################################################
## READ .gpr file, returns a list of measures
## extracted from the .gpr file
###################################################

## No background substraction
## No normalization
## Suppose that fnames is a .gpr file.
## Returns a list of data needed for quality
## Combines read.marrayRaw, maQualityMain, gpTools
## One slide only!!!

readGPR <- function (fnames = NULL, path= ".", echo=TRUE, skip = 0,
                     sep ="\t", quote= "",...)
  {
    print("Starting readGPR")
    print(path)
    # Test if input data is OK
    if (is.null(path))
      path <- "."
 
    if (missing(fnames) | is.null(fnames)) 
      fnames <- dir(path, pattern = "*\\.gpr$")
    fullfnames <- file.path(path, fnames)
    f <- fullfnames[1]
    opt <- list(...)
    
    # Search in the gpr file where are the colums starting
    y <- readLines(fullfnames[1], n = 100)
    skip <- grep("F635 Median", y)[1] - 1
    print(skip)
  
    # Read data in       
    name <- id <- NULL
    block <- column <- row <- NULL
    Date <- PMTR <- PMTG <- normMethod <- normCoef <- NULL
    Gfmedian <- Gfmean <- Gfsd <- NULL
    Gbmedian <- Gbmean <- Gbsd <- NULL
    Rfmedian <- Rfmean <- Rfsd <- NULL
    Rbmedian <- Rbmean <- Rbsd <- NULL
    Gfsat <- Rfsat <- NULL
    spotDia <- spotFArea <- spotBArea <- NULL
    Flags <- NULL
    
    tmp <- readLines(f, n = 40)
    if (length(grep("DateTime", tmp)) != 0) 
      Date <- gsub("\"", "",
                   strsplit(tmp[grep("DateTime", tmp)], split = "=")[[1]][2])
    if (length(grep("PMT", tmp)) != 0) {
      PMT <- gsub("\"", "",
                   strsplit(tmp[grep("PMT", tmp)], split = "=")[[1]][2])
      PMT <- strsplit(PMT, split="\t")
      PMTR <- PMT[[1]][1]
      PMTG <- PMT[[1]][2]

    }
    if (length(grep("NormalizationMethod", tmp)) != 0)
      normMethod <- gsub("\"", "",
                         strsplit(tmp[grep("NormalizationMethod",tmp)],
                                  split = "=")[[1]][2])
    if (length(grep("NormalizationFactor", tmp)) != 0)
        normCoef1 <- gsub("\"", "",
                         strsplit(tmp[grep("NormalizationFactor", tmp)],
                                  split = "=")[[1]][2])
    normCoef <- gsub("\t", ",", normCoef1)

    print(paste("Reading", f))
    h <- strsplit(readLines(f, n = skip + 1), split = sep)
    h <- as.list(unlist(h[[length(h)]]))
    names(h) <- gsub("\"", "", unlist(h))
    dat <- scan(f, quiet = TRUE, what = h, sep = sep, skip = skip + 
                1, quote = quote, ...)

    block <- cbind(block, dat[["Block"]])
    column <- cbind(column, dat[["Column"]])
    row <- cbind(row, dat[["Row"]])
    name <- cbind(name, dat[["Name"]])
    id <- cbind(id, dat[["ID"]])
    Gfmedian <- cbind(Gfmedian, as.numeric(dat[["F532 Median"]]))
    Gfmean <- cbind(Gfmean, as.numeric(dat[["F532 Mean"]]))        
    Gfsd <- cbind(Gfsd, as.numeric(dat[["F532 SD"]]))
    Gbmedian <- cbind(Gbmedian, as.numeric(dat[["B532 Median"]]))
    Gbmean <- cbind(Gbmean, as.numeric(dat[["B532 Mean"]]))        
    Gbsd <- cbind(Gbsd, as.numeric(dat[["B532 SD"]]))
    Rfmedian <- cbind(Rfmedian, as.numeric(dat[["F635 Median"]]))
    Rfmean <- cbind(Rfmean, as.numeric(dat[["F635 Mean"]]))
    Rfsd <- cbind(Rfsd, as.numeric(dat[["F635 SD"]]))
    Rbmedian <- cbind(Rbmedian, as.numeric(dat[["B635 Median"]]))
    Rbmean <- cbind(Rbmean, as.numeric(dat[["B635 Mean"]]))
    Rbsd <- cbind(Rbsd, as.numeric(dat[["B635 SD"]]))        
    Gfsat <- cbind(Gfsat, as.numeric(dat[["F532 % Sat."]]))
    Rfsat <- cbind(Rfsat, as.numeric(dat[["F635 % Sat."]]))
    spotDia <- cbind(spotDia, as.numeric(dat[["Dia."]]))
    spotFArea <- cbind(spotFArea, as.numeric(dat[["F Pixels"]]))
    spotBArea <- cbind(spotBArea, as.numeric(dat[["B Pixels"]]))
    Flags <- cbind(Flags,as.numeric(dat[["Flags"]]))

    gprData <- list(File = f, Date = Date, PmtR = PMTR, PmtG = PMTG,
                    Normalization = normMethod, NormCoefficient = normCoef,
                    Name = name, ID = id,
                    Block = block, Column = column, Row = row,
                    GfMedian = Gfmedian, GfMean = Gfmean, GfSD = Gfsd,
                    GbMedian = Gbmedian, GbMean = Gbmean, GbSD = Gbsd,
                    RfMedian = Rfmedian, RfMean = Rfmean, RfSD = Rfsd,
                    RbMedian = Rbmedian, RbMean = Rbmean, RbSD = Rbsd,
                    GfSaturation = Gfsat, RfSaturation = Rfsat,
                    SpotDiameter = spotDia, spotArea = spotFArea,
                    bgArea = spotBArea, Flags = Flags)
    
    # Result
    return(gprData)
  }


###################################################
## Given a gpr file
## Computes needed statistics to assess quality
###################################################

## Call readGRP.R
## One slide only

slideQuality <- function(slideFile = NULL, path=".",
                          output = TRUE, echo=TRUE,...)
  {
    if (echo) print("SlideQuality starting")

    # Check input argument
    if (missing(slideFile) | is.null(slideFile))
      slideFile <-  dir(path, pattern = "*\\.gpr$")

    slide <- slideFile[1]
    print(slide)
    print(path)
    
    # Read data in
    gprData <- readGPR(fnames=slide, path=path)

    if (echo) print("SlideQuality 1")

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
                                
    if (echo) print("SlideQuality 3")
    
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
    if (echo) print("SlideQuality 4")
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
    if (echo) print("SlideQuality 5")
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
    if (echo) print("SlideQuality 6")

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
    if (echo) print("SlideQuality 7")
    goodMed <-  !(is.infinite(Mmedian) | is.na(Mmedian))

    pt <- cbind(Block = as.vector(as.numeric(gprData[["Block"]])),
                M = Mmedian)
    #get the mean of Mvalues by print-tip group
    printTip <-  by(pt[,2],factor(pt[,1]), mean)
    goodPt <- !(is.infinite(printTip) | is.na(printTip))
    msePtip <- MSE(printTip[goodPt], center=0)


    if (echo) print("SlideQuality 8")

    # Controls
    Control <- arrayControls(gprData)
    numE <- numNeg <- numPos <- numProb <- 0
    emp <- neg <- pos <- NULL

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

    EmptyMed <- median(emp, na.rm=TRUE)
    NegativeMed <- median(neg, na.rm=TRUE)
    PositiveMed <- median(pos, na.rm=TRUE)

    difEmptyNegative <- EmptyMed - NegativeMed
    difPositiveNegative <- PositiveMed - NegativeMed

    # Replicates
    if (echo) print("SlideQuality 9")

    gId <- gprData[["ID"]]
    Replicates <- arrayReplicates(gprData)
  
    index <- NULL
    for(r in Replicates){
      for(i in 1:length(gId)){
        if (r == gId[i]) index <- c(index,i)
      }
    }
    
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

    sortedMeasures <- c(#"File", "Date",
                        "range RF", "range GF", "range RB", "range GB",
                        "Median empty ctrl", "Median negative ctrl",
                        "Median positive ctrl",
                        "MSE by print-tip",
                        "MSE lowess",
                        "Spot radius",
                        "Percentage of flagged spots",
                        "RB variance","GB variance",
                        "RB mad", "GB mad", "RB IQR", "GB IQR",
                        "RS2N", "GS2N",
                        "MMRmad","Percentage of spots MMRmad > 0.5",
                        "MMR IQR",
                        "Difference empty/negative",
                        "Difference positive/negative",
                        "Variance of replicated spots")
        

    sortedRes <- c(#gprData[["File"]], gprData[["Date"]],
                   rangeRf, rangeGf, rangeRb, rangeGb,
                   EmptyMed, NegativeMed, PositiveMed,
                   msePtip, mseFit,
                   spotRadius,
                   percentFlag,
                   RBvar, GBvar,
                   RbMad, GbMad, RbIqr, GbIqr,
                   RS2Nmedian, GS2Nmedian,
                   MMRmad, percentSpotOverMmrLim, mmrIqr,
                   difEmptyNegative, difPositiveNegative,
                   varRepA
                   )
    
    numResult <- cbind(sortedMeasures, sortedRes)
                
    if (echo) print("SlideQuality done...")
    return(numResult)         
  }


###################################################
## Given a gpr files
## Computes needed statistics to assess quality
###################################################

## Takes all .gpr file into acccount

globalQuality <- function (fnames = NULL, path = ".",
                           output = TRUE, echo = TRUE,...)
  {

    # Check input arguments
    if (echo) print("Starting global qality")
    if (missing(path) | is.null(path))
      path <- "."

    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.gpr$")

    # Prepares results
    quality <- NULL
    tmp <- NULL

    sortedMeasures <- c(#"File", "Date",
                        "range RF", "range GF", "range RB", "range GB",
                        "Median empty ctrl", "Median negative ctrl",
                        "Median positive ctrl",
                        "MSE by print-tip",
                        "MSE lowess",
                        "Spot Radius",
                        "Percentage of flagged spots",
                        "RB variance","GB variance",
                        "RB mad", "GB mad", "RB IQR", "GB IQR",
                        "RS2N", "GS2N",
                        "MMRmad","Percentage of spots MMRmad > 0.5",
                        "MMR IQR",
                        "Difference empty/negative",
                        "Difference positive/negative",
                        "Variance of replicated spots")
    
      
    quality <- sortedMeasures

    # Call to slideQuality for each gpr file
    for (f in fnames)
      {
        if (echo) print("In the loop ")
        tmp <- slideQuality(slideFile = f, path=path)
        quality <- cbind(quality, tmp[,2])
      }
    colnames(quality) <- c("measure",fnames)
    # Results
    write.table(quality, "quality.txt",sep="\t")
                #row.names=FALSE, col.names=TRUE, sep="\t")

    return(quality)                                 
  }


###############################################
## Given a set of reference slides
## Returns mean, variance and iqr for each measure
##
## To be used to create reference table for scaling 
##############################################

qualRefTable <- function(fnames=NULL, path=".")
  {
    print("Starting qualRefTable")
    if (missing(path) | is.null(path))
      path <- "."

    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.gpr$")

    qt <- globalQuality(fnames=fnames, path=path)
    tab <- qt[,2:dim(qt)[2]]
        
    ave <- character(0)
    sigma <- character(0)
    iqr <-  character(0)
       
    for(i in 1:dim(tab)[1])
      {
        ave <-c(ave,
               mean(as.numeric(tab[i,]), na.rm=TRUE))
        sigma <- c(sigma, var(as.numeric(tab[i,]), na.rm=TRUE))
        iqr <- c(iqr, IQR(as.numeric(tab[i,]), na.rm=TRUE))
      }
      
    reftab <- cbind(mean=as.numeric(ave),
                    sigma=as.numeric(sigma),
                    iqr=as.numeric(iqr))
    return(reftab)
  }

#reftab <- qualRefTable(path="C:/Mydoc/Projects/quality/DemoA")

## Takes a matrix/vector of numbers as argument
## Scales this matrix according to values in reftab
## To use with globalQuality, remove 1st column (text)
## of globalQuality result

arrayScal <- function(numMat, scalingData=NULL)
{
  print("Starting arrayScal")
  if(missing(numMat))
    stop("Input error, matrix to scale missing")
  else
    {
      if(missing(scalingData)|is.null(scalingData))
        {
          if(!("reftab" %in% ls(1)))
          reftab <- qualRefTable(path=path)
        }
      else
        reftab <- scalingData
      #print("Dim reference table")
      #print(dim(reftab))
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
          numMat <- as.matrix(numMat)
        }
      #print(nrow)
      #print(ncol)
      #Check dimensions 
      if(nrow!=dim(reftab)[1])
        stop("Error is scaling, files of different length")
      else
        {
          tmp <- numMat
          for(j in 1:ncol)
            for(i in 1:nrow)
            {
              if(reftab[i,"sigma"]!=0)
                tmp[i,j] <- (as.numeric(tmp[i,j]) - as.numeric(reftab[i,"mean"]))/as.numeric(reftab[i,"sigma"])
              else
                {
                  print("One of the ranges = 0")
                  tmp[i,j] <- as.numeric(tmp[i,j]) - as.numeric(reftab[i,"mean"])
                }              
            }
          return(tmp)
        }
    }
}


#################################
## Plot function

# 1. Plot the boxplot
# 2. SuperImpose value for slides of interest (1 or more)

#arrayQuality and reference are results from globalQuality

qualPlot <- function(arrayQuality,reference)
  {
    print("starting plot")
    # Reference = output of globalQuality for ref slides
    scalref <- arrayScal(reference[,-1])
    scalarray <- arrayScal(arrayQuality[,-1])

    if(is.null(dim(scalref))|nrow(arrayQuality)!=nrow(reference))
      stop("Input must be a matrix resulting from globalQuality.R")
    else
      {
        print("Boxplot")
        #Boxplot of reference arrays quality measure
        nr <- nrow(scalref)
        lim <- range(as.numeric(scalarray),as.numeric(scalref),na.rm=TRUE)
        plot(0:(nr+1),0:(nr+1),ylim=lim,type="n", axes=FALSE,ylab="",xlab="")
        axis(1,1:nr,labels=reference[,1],las=2)
        axis(2)

        for(i in 1:nr)
          boxplot(as.numeric(scalref[i,]), at=i, add=TRUE)

        #Line for tested arrays
        
        nc <- ncol(scalarray)
        col <- rainbow(nc)
        print("lines")
        for(i in 1:nc)
          lines(1:nr, scalarray[,i],col=col[i])

        #legend
        leg <- colnames(arrayQuality)
        leg.txt <- leg[-1]
        legend(1,60,leg.txt,lty=1,col=col)
      }
    
  }


