#####################################################
## Set of function written to estimate
## quality of microarray slides
##
## Author: Agnes Paquet
## Modified: 04/16/2004


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

## returns a list of vector containing info from grp file column
readGPR <- function (fnames = NULL, path= ".", DEBUG=FALSE, skip = 0,
                     sep ="\t", quote= "\"",...)
  {
    if (DEBUG) print("Starting readGPR")
    if (DEBUG) print(path)
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
    if (DEBUG) print(skip)
  
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

    if (DEBUG) print(paste("Reading", f))
    h <- strsplit(readLines(f, n = skip + 1), split = sep)
    h <- as.list(unlist(h[[length(h)]]))
    names(h) <- gsub("\"", "", unlist(h))
    dat <- scan(f, quiet = TRUE, what = h, sep = sep, skip = skip + 
                1, quote = quote, ...)

    block <- cbind(block, as.numeric(dat[["Block"]]))
    column <- cbind(column, as.numeric(dat[["Column"]]))
    row <- cbind(row, as.numeric(dat[["Row"]]))
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

    gprData <- list(File = fnames[1], Date = Date, PmtR = PMTR, PmtG = PMTG,
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

## Argument: result of readGPR
## Returns: matrix of numbers

slideQuality <- function(gprData=NULL, DEBUG=FALSE,...)
  {
    if (DEBUG) print("SlideQuality starting")
    
    # Check input argument
    if (is.null(gprData))
      stop("Slide information is missing:")
    
    # Read data in

    if (DEBUG) print("SlideQuality 1")

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
    Control <- arrayControls(gprData)
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
                        "range RF", "range GF",
                        "- RB mad", "- GB mad",
                        "Median RS2N", "Median GS2N",
                        #"range RB", "range GB",
                        "- Median A for empty ctrl",
                        "- Median A for neg ctrl",
                        "Median A for positive ctrl",
                        "Pos ctl median A - Neg ctl median A",
                        "- Var replicated spots A values",
                        "- Mvalues MSE by print-tip",
                        "- MSE lowess",
                        #"Spot radius",
                        #"Percentage of flagged spots",
                        "- % flagged spots",
                        #"RB variance","GB variance",
                        #"RB IQR", "GB IQR",
                        "- Mvalues MMRmad",
                        "- % spots with Mvalues MMRmad>0.5",
                        #"Percentage of spots MMRmad > 0.5",
                        #"MMR IQR",
                        #"Difference empty/negative"
                        )
        

    sortedRes <- c(#gprData[["File"]], gprData[["Date"]],
                   rangeRf, rangeGf,
                   -RbMad, -GbMad,
                   RS2Nmedian, GS2Nmedian,
                   #rangeRb, rangeGb,
                   -EmptyMed, -NegativeMed, PositiveMed,
                   difPositiveNegative,
                   -varRepA,
                   -msePtip, -mseFit,
                   #spotRadius,
                   -percentFlag,                                      
                   -MMRmad, -percentSpotOverMmrLim,
                   #RBvar, GBvar,
                   #RbIqr, GbIqr,
                   #mmrIqr,
                   #difEmptyNegative,
                   
                   )
    
    #numResult <- cbind(sortedMeasures, sortedRes)
    numResult <- as.matrix(sortedRes)
    rownames(numResult) <- sortedMeasures
    colnames(numResult) <- gprData[["File"]]
                
    if (DEBUG) print("SlideQuality done...")
    return(numResult)         
  }


###################################################
## Given a gpr files
## Computes needed statistics to assess quality
###################################################

## Takes all .gpr file into acccount
## Plot quality boxplot and diagnostic plots
## creates html report
## if output: writes quality and normalized data to file
## return quality measure and marrayRaw object in a list

#############
## Example:
## test <- gpQuality(path="C:/Mydoc/Projects/quality/DemoA/", resdir="QualPlot")

## Reference = output from globalQuality
## ScalingTable: output from qualRef
## Must be run on the same gpr files
gpQuality <- function(fnames = NULL, path = ".",
                      organism=c("Mm", "Hs"),
                      compBoxplot = TRUE,
                      reference=NULL,
                      scalingTable=NULL,
                      output=FALSE,
                      resdir=".",
                      dev="png", #set default to be png 
                      DEBUG = FALSE,...)
{

    print("Starting gpQuality...")
    
    # Check input arguments
    if (missing(path) | is.null(path))
      path <- getwd()
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.gpr$")

    organism <- organism[1]
    opt <- list(...)
    
 ###################
 ## Setting up output device

    if (DEBUG) print("Name of output device")
    plotdef <- switch(dev,
                      "bmp" = list(dev=list(width=800, height=600, bg="white"), suffix="bmp"),
                      "jpeg" = list(dev=list(quality=100, width=800, height=600, bg="white"), suffix="jpeg"),
                      "jpg" =  list(dev=list(quality=100, width=800, height=600, bg="white"), suffix="jpeg"),
                      "postscript" = list(dev=list(paper="special", width=8, height=6, bg="white"), suffix="ps"),
                      "png" =  list(dev=list(width=800, height=600, bg="white"), suffix="png"),
                      list(dev=list(width=800, height=600,bg="white"), suffix="png"),
                    )
    if(!is.element(dev, c("bmp", "jpeg","png","postscript","jpg")))
      print("Format error, format will be set to PNG")

    if (DEBUG) print(paste("compBoxplot = ",compBoxplot, sep=""))
    
    if (compBoxplot)
      {
    
        # Prepares results
        quality <- NULL
        tmp <- NULL
        QCp <- c()
        Dp <- c()
        
        curdir <- getwd()
        if (!file.exists(resdir))
          dir.create(resdir)
        if (DEBUG) print(getwd())
        
        #Allocation of matrix for marrayraw object

        if(DEBUG) print(path)
        if(DEBUG) print(resdir)
        f <- fnames[1]
        gp <- readGPR(fnames=f, path=path)
        nrow <- length(gp[["RfMedian"]])
        ncol <- length(fnames)
        
        mlayout <- maCompLayout(as.matrix(cbind(gp[["Block"]],
                                                gp[["Row"]], gp[["Column"]])))
        tmp <- new("marrayInfo", maInfo=data.frame(gp[["Name"]],gp[["ID"]]))
        mlayout@maControls <- as.factor(maGenControls(tmp))
        rm(f, gp)
        
        Rf <- Gf <- Rb <- Gb <- weight <- matrix(0,nrow=nrow, ncol=ncol)
        filenames <- c()
        nb <- c()
        
        # Call to slideQuality for each gpr file

        for (i in 1:length(fnames))
          {
            if (DEBUG) print("In the loop ")
            f <- fnames[i]
            gp <- readGPR(fnames = f, path=path)
            restmp <- slideQuality(gp, DEBUG=DEBUG)
            scal <- arrayScal(restmp, organism=organism, scalingData=scalingTable)
            
            ###start plot
            Gf[,i] <- gp[["GfMedian"]]; Gb[,i] <-gp[["GbMedian"]] 
            Rf[,i] <- gp[["RfMedian"]]; Rb[,i] <- gp[["RbMedian"]]
            weight[,i] <- gp[["Flags"]]
            filenames <- c(filenames, gp[["File"]])
        
            #qualBoxplot
            setwd(resdir)
            plotname <- paste("qualPlot",unlist(strsplit(f, ".gpr")), dev,sep=".")
            plotdef <- c(plotdef, list(main=paste(f, ": Quantitative Diagnostic Plots")))
            do.call(dev, maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)))
            par(mar=c(3,14,2,6))
            nbtmp <- qualBoxplot(restmp, reference=reference, scalingTable=scalingTable)
            dev.off()
            setwd(curdir)

            #nb = matrix ncol=2, nrow=length(fnames)
            nb <- rbind(nb, nbtmp)
            QCp <- c(QCp, plotname)
            if (DEBUG) print("End of plot")
            if (DEBUG) print(paste("save as ",plotname))
            if (DEBUG) print("Binding results")
            quality <- cbind(quality, restmp[,1])
            meas <- rownames(restmp)
          }
        
        print("Comparative plots done")
        
        #Create marrayRaw
        colnames(Gf) <- colnames(Gb) <- colnames(Rf) <- colnames(Rb) <- filenames
        if (DEBUG) print("building mraw")
        mraw <- new("marrayRaw", maRf=Rf, maGf=Gf, maRb=Rb,
                    maGb=Gb, maNotes="", maLayout=mlayout,
                    maW=weight, maGnames=tmp)
        #maQualityPlots
        setwd(resdir)
        print("Starting maQualityPlots")
        maQualityPlots(mraw, DEBUG=DEBUG)
        
        #get diagnostic plots names
        tmpname <- sub(".gpr", "",colnames(mraw@maGf))
        pn <- paste("diagPlot", tmpname, sep=".")

        dirfiles <- dir(".")
        
        for(i in 1:length(pn))
          Dp <- c(Dp, dirfiles[grep(pn[i], dirfiles)[1]])
        
        if (DEBUG) print("After for loop")
        if (DEBUG) print(nb)
        quality2HTML(fnames=fnames,path=resdir, QCplot=QCp, DiagPlot=Dp,nbBelow=nb)
        
        colnames(quality) <- fnames
        rownames(quality) <- meas
        
        print("gpQuality done")
        
       ####### Results
        
        if (output)
          {
            print("Printing results to file")
            write.table(quality, "quality.txt",sep="\t", col.names=NA)
            colnames(mraw@maGnames@maInfo) <- c("Name", "ID")
            outputNormData(mraw)
          }

        setwd(curdir)
        return(list(mraw=mraw, quality=quality))
      }

   else {

     curdir <- getwd()
     if (!file.exists(resdir))
       dir.create(resdir)
     if (DEBUG) print(getwd())

     
     mraw <- read.GenePix(fnames, path)
     colnames(maGf(mraw)) <- fnames
     print("Starting maQualityPlots")
     setwd(resdir)
     maQualityPlots(mraw, DEBUG=DEBUG)

     if (output)
       {
         print("Printing results to file")
         outputNormData(mraw)
       }

     setwd(curdir)
     return(list(mraw=mraw))     
   }

    
  }






###############################################
## Given a set of reference slides
## Returns mean and iqr for each measure
##
## To be used to create reference table for scaling 
##############################################

qualRefTable <- function(fnames=NULL, path=".")
  {
    if (missing(path) | is.null(path))
      path <- "."

    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.gpr$")

    qt <- globalQuality(fnames=fnames, path=path)
    tab <- qt
    ave <- character(0)
    sigma <- character(0)
    iqr <-  character(0)
       
    for(i in 1:dim(tab)[1])
      {
        ave <-c(ave,
               median(as.numeric(tab[i,]), na.rm=TRUE))
        iqr <- c(iqr, IQR(as.numeric(tab[i,]), na.rm=TRUE))
      }
      
    reftab <- as.matrix(cbind(mean=as.numeric(ave),
                              iqr=as.numeric(iqr)))
    rownames(reftab) <- rownames(qt)
    return(reftab)
  }


## Takes a matrix/vector of numbers as argument
## Scales this matrix according to values in scalingData
## if scalingData is missing, use scaling matrix (Mm by default, or Hs)
## of globalQuality result

arrayScal <- function(numMat, scalingData=NULL, organism=c("Mm", "Hs"))
{
  organism=organism[1]

  if(missing(numMat))
    stop("Input error, matrix to scale missing")
  else
    {
      if(missing(scalingData)|is.null(scalingData))
        {
          if(organism == "Mm")
            {
              if(!("MmScalingTable" %in% ls(1)))
                data(MmScalingTable)
              reftab <- MmScalingTable
            }
          else
            {
              if(!("HsScalingTable" %in% ls(1)))
                data(HsScalingTable)              
              reftab <- HsScalingTable
            }
        }
      else
        reftab <- scalingData
      
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
          print(dim(reftab))
          stop("Error is scaling, files of different length")
        }
      else
        { 
          tmp <- numMat
          for(j in 1:ncol)
            for(i in 1:nrow)
            {
              if(reftab[i,"iqr"]!=0)
                tmp[i,j] <- (as.numeric(tmp[i,j]) - as.numeric(reftab[i,"mean"]))/as.numeric(reftab[i,"iqr"])
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


################################################
## HTML report
# Read in names of plots
# suppose directory contains
# same number of qualPlot and diagPlot
# based on alphabetical order

# nbBelow is a matrix: for each slide: number of measure below range of good slides and total number of measures

quality2HTML <- function(fnames=NULL, path=".", DiagPlot=NULL, QCplot=NULL, resdir=".", nbBelow=NULL)
  {

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
    datadir <- system.file("data", package="arrayQuality")
    html <- paste(readLines(file.path(datadir, "index.html")), "\n", collapse="")

    split <- unlist(strsplit(html, split="<a name=\"table\"></a>"))

    tab <- tableTag(fnames,fullQCp, fullBoxp)
    
    cat(split[1], tab,split[2], file=outputfile)     
    close(outputfile)
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

# Reads in gpr files
# returns matrix of QC measures
globalQuality <- function(fnames = NULL, path = ".",
                          organism=c("Mm", "Hs"),
                          output=FALSE,
                          resdir=".",
                          DEBUG = FALSE,...)
  {
    # Check input arguments
    if (DEBUG) print("Starting globalQuality")
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.gpr$")
    
    organism <- organism[1]
        
    # Prepares results
    quality <- NULL
    
    # Call to slideQuality for each gpr file

    for (i in 1:length(fnames))
      {
        if (DEBUG) print("In the loop ")
        f <- fnames[i]
        gp <- readGPR(fnames = f, path=path)
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


qualBoxplot <- function(arrayQuality=NULL,  reference=NULL, scalingTable=NULL, organism=c("Mm", "Hs"),...)
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
          reference <- data(HsReferenceDB)
      }

    scalref <- arrayScal(reference, scalingData=scalingTable)
    scalarray <- arrayScal(arrayQuality, scalingData=scalingTable)
    
    if(is.null(dim(scalref))|nrow(arrayQuality)!=nrow(reference))
      stop("Input must be a matrix resulting from slideQuality.R")
    else
      {
        #Boxplot of reference arrays quality measure

        nr <- nrow(scalref)

        goodLim <- matrix(0,nrow=nr,ncol=1)

        lim <- range(as.numeric(scalarray),as.numeric(scalref),na.rm=TRUE)
        plot(0:(nr+1),0:(nr+1),xlim=lim,type="n", axes=FALSE,ylab="",xlab="", xaxt="n",
             main="Array Quality Control Comparison")

        axis(2, at=1:nr, labels=rownames(arrayQuality), cex.axis=0.7, las=2)
        tmp <- c(paste(as.character(round(score)), "(", as.character(round(arrayQuality,1)), ")"), "% (value)")
        axis(4, at=1:(nr+1), labels=, tmp, cex.axis=0.8, las=2)
        axis(1, at=c(lim[1]+0.5, lim[2]-0.5), labels= c("Problematic", "Good"))

        for(i in 1:nr)
          {
            bp <- boxplot(as.numeric(scalref[i,]), at=i, add=TRUE, horizontal=TRUE, axes=FALSE)
            goodLim[i,1] <- bp$stats[1]
            
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
                if(as.numeric(scalarray[j,i]) < goodLim[j,1] ||
                   is.na(scalarray[j,i]))
                  isBelowLim[j,i] <- TRUE
              }
          }
        box()
        #legend
        leg.txt <- colnames(arrayQuality)
        legend(lim[1],nr+1,leg.txt,lty=1,col=col, cex=0.7)

        res <- c()
        for(i in 1:nc)
          res <- c(res, length(isBelowLim[isBelowLim[,i],i]))
        return(cbind(res, nr))
      }
  }

outputNormData <- function(mraw)
  {
    mnorm <- maNorm(mraw)
    write.marray(mnorm, "NormalizedData.xls")
  }
