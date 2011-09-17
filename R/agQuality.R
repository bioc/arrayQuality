###################################################
## READ Agilent file .txt file, returns a list of measures
## extracted from the Agilent file
## source("~/Projects/madman/Rpacks/arrayQuality/R/agQuality.R")
##
###################################################

## No background substraction
## No normalization
## Suppose that fnames is a .gpr file.
## Returns a list of data needed for quality
## Combines read.marrayRaw, maQualityMain, gpTools
## One slide only!!!

## returns a list of vector containing info from grp file column
readAgilent <- function (fnames = NULL, path= ".", DEBUG=FALSE, skip = 0,
                     sep ="\t", quote= "\"", controlId = c("ProbeName"), ...)
  {
    if (DEBUG) print("Starting readAgilent")
    if (DEBUG) print(path)
    # Test if input data is OK
    if (is.null(path))
      path <- "."
 
    if (missing(fnames) | is.null(fnames)) 
      fnames <- dir(path, pattern = ".*\\.txt$")
    fullfnames <- file.path(path, fnames)
    f <- fullfnames[1]
    opt <- list(...)
    
    # Search in the Agilent file where are the colums starting
    y <- readLines(fullfnames[1], n = 100)
##    skip <- grep("gMedianSignal", y)[1] - 1
    skip <- grep("gMedianSignal", y)[1]
    
    if (DEBUG) print(skip)
  
    # Read data in       
    name <- id <- c()
    block <- column <- row <- c()
    Date <- PMTR <- PMTG <- normMethod <- normCoef <- c()
    Gfmedian <- Gfmean <- Gfsd <- c()
    Gbmedian <- Gbmean <- Gbsd <- c()
    Rfmedian <- Rfmean <- Rfsd <- c()
    Rbmedian <- Rbmean <- Rbsd <- c()
    Gfsat <- Rfsat <- c()
    spotDia <- spotFArea <- spotBArea <- c()
    Flags <- c()
    
    tmp <- readLines(f, n = 40)
    if (length(grep("DateTime", tmp)) != 0) 
      Date <- unlist(strsplit(readLines(f, n=3)[3], split="\t"))[4]
    if (length(grep("PMT", tmp)) != 0) {
      PMT <- gsub("\"", "",
                   strsplit(tmp[grep("PMT", tmp)], split = "=")[[1]][2])
      PMT <- strsplit(PMT, split="\t")
      PMTR <- PMT[[1]][1]
      PMTG <- PMT[[1]][2]
    }
    
##    if (length(grep("NormalizationMethod", tmp)) != 0)
##      normMethod <- gsub("\"", "",
##                        strsplit(tmp[grep("NormalizationMethod",tmp)],
##                                  split = "=")[[1]][2])
##    if (length(grep("NormalizationFactor", tmp)) != 0)
##        normCoef1 <- gsub("\"", "",
##                         strsplit(tmp[grep("NormalizationFactor", tmp)],
##                                  split = "=")[[1]][2])
##    normCoef <- gsub("\t", ",", normCoef1)

    if (DEBUG) print(paste("Reading", f))
    h <- strsplit(readLines(f, n = skip), split = sep)
    h <- as.list(unlist(h[[length(h)]]))
    names(h) <- gsub("\"", "", unlist(h))
##    dat <- scan(f, quiet = TRUE, what = h, sep = sep, skip = skip + 1,
##                quote = quote, ...)
    dat <- scan(f, quiet = TRUE, what = h, sep = sep, skip = skip,
                quote = quote, ...)

##    block <- cbind(block, as.numeric(dat[["Block"]]))
    column <- cbind(column, as.numeric(dat[["Col"]]))
    row <- cbind(row, as.numeric(dat[["Row"]]))
    name <- cbind(name, dat[["Description"]])
    id <- cbind(id, dat[[controlId]])
    Gfmedian <- cbind(Gfmedian, as.numeric(dat[["gMedianSignal"]]))
    Gfmean <- cbind(Gfmean, as.numeric(dat[["gMeanSignal"]]))        
    Gfsd <- cbind(Gfsd, as.numeric(dat[["gPixSDev"]]))
    Gbmedian <- cbind(Gbmedian, as.numeric(dat[["gBGMedianSignal"]]))
    Gbmean <- cbind(Gbmean, as.numeric(dat[["gBGMeanSignal"]]))        
    Gbsd <- cbind(Gbsd, as.numeric(dat[["gBGPixSDev"]]))
    Rfmedian <- cbind(Rfmedian, as.numeric(dat[["rMedianSignal"]]))
    Rfmean <- cbind(Rfmean, as.numeric(dat[["rMeanSignal"]]))
    Rfsd <- cbind(Rfsd, as.numeric(dat[["rPixSDev"]]))
    Rbmedian <- cbind(Rbmedian, as.numeric(dat[["rBGMedianSignal"]]))
    Rbmean <- cbind(Rbmean, as.numeric(dat[["rBGMeanSignal"]]))
    Rbsd <- cbind(Rbsd, as.numeric(dat[["rBGPixSDev"]]))        
    Gfsat <- cbind(Gfsat, as.numeric(dat[["gNumSatPix"]]))
    Rfsat <- cbind(Rfsat, as.numeric(dat[["rNumSatPix"]]))
##    spotDia <- cbind(spotDia, as.numeric(dat[["Dia."]]))  Not used
    spotFArea <- cbind(spotFArea, as.numeric(dat[["gNumPix"]]))
    spotBArea <- cbind(spotBArea, as.numeric(dat[["gBGNumPix"]]))
    Flags <- cbind(Flags, as.numeric(as.logical(dat[["gIsFeatNonUnifOL"]]) & as.logical(dat[["rIsFeatNonUnifOL"]]))) 
    gprData <- list(File = fnames[1], Date = Date, PmtR = PMTR, PmtG = PMTG,
                    Normalization = NA, NormCoefficient = NA,
                    Name = name, ID = id,
                    Block = 1, Column = column, Row = row,
                    GfMedian = Gfmedian, GfMean = Gfmean, GfSD = Gfsd,
                    GbMedian = Gbmedian, GbMean = Gbmean, GbSD = Gbsd,
                    RfMedian = Rfmedian, RfMean = Rfmean, RfSD = Rfsd,
                    RbMedian = Rbmedian, RbMean = Rbmean, RbSD = Rbsd,
                    GfSaturation = Gfsat, RfSaturation = Rfsat,
                    SpotDiameter = NA, spotArea = spotFArea,
                    bgArea = spotBArea, Flags = Flags)
    
    # Result
    return(gprData)
  }


###########################################################################
##
##  Move controlCode and maGenControls from marrayTools to marrayClasses
##  May 7, 2003
##  Sept 21, 2004
## TO USE:
## 1. Create Reference slides
##    reference <- globalQuality(fnames, inputsource = "readAgilent")
## 2. Run : some alternative
##        -- agQuality(fnames, reference = reference)  ## Use what you have created
##        -- agQuality(fnames, reference = reference, compBoxplot=FALSE)
##                                    ## Only generate qualitative plots
##        -- agQuality(fnames)  ## Use exisiting reference info
###########################################################################

agcontrolCode <-
structure(c("\\(+\\)*", "\\(-\\)*", "Positive", "Negative"),.Dim = c(2, 2), .Dimnames = list(c("1", "2"), c("Pattern", "Name")))

agQuality <- function(fnames = NULL, path = ".",
                      organism = c("Mm", "Hs"),
                      compBoxplot = TRUE,
                      reference = NULL,
                      controlMatrix = agcontrolCode,
                      controlId = c("ProbeName"),
                      output = FALSE,
                      resdir = ".",
                      dev = "png", #set default to be png 
                      DEBUG = FALSE,...)
{

    print("Starting agQuality...")
    
    # Check input arguments
    if (missing(path) | is.null(path))  path <- "."
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = ".*\\.txt$")

    organism <- organism[1]
    controlId <- controlId[1]
    if (DEBUG) print(controlId)
    if (DEBUG) print(controlMatrix)
    opt <- list(...)

    #get Normalization method if any
    #norm.defs <- maDotsDefaults(opt, list(norm="p"))
    #defs <- list(norm="p")
    #norm.defs <- maDotsMatch(maDotsDefaults(opt, defs), formals(args("maNorm")))

 ###################
 ## Setting up output device

    if (DEBUG) print("Name of output device")
    plotdef <- switch(dev,
                      "bmp" = list(dev=list(width=800, height=600, bg="white"), suffix="bmp"),
                      "jpeg" = list(dev=list(quality=100, width=800, height=600, bg="white"), suffix="jpeg"),
                      #"jpg" =  list(dev=list(quality=100, width=800, height=600, bg="white"), suffix="jpeg"),
                      #"postscript" = list(dev=list(paper="special", width=8, height=6, bg="white"), suffix="ps"),
                      "postscript" = list(dev=list( bg="white"), suffix="ps"),
                      "png" =  list(dev=list(width=800, height=600, bg="white"), suffix="png"),
                      list(dev=list(width=800, height=600,bg="white"), suffix="png"))
    if(!is.element(dev, c("bmp", "jpeg","png","postscript","jpg")))
      {
        print("Format error, format will be set to PNG")
        dev = "png"
      }

      # was print("Format error, format will be set to PNG")

    if (DEBUG) print(paste("compBoxplot = ", compBoxplot, sep=""))

    
##############################################################
    ## COMPBOXPLOT = TRUE
##############################################################
    
    if (compBoxplot)
      {

        if(DEBUG) print("Starting compBoxplot")
        
        # Prepares results
        quality <- c()  
        tmp <- c()
        QCp <- c()
        Dp <- c()
        
        curdir <- getwd()
        if (!file.exists(resdir))
          dir.create(resdir)
        if (DEBUG) print(getwd())
        
        #Allocation of matrix for marrayraw object

        if (DEBUG) print(path)
        if (DEBUG) print(resdir)
#        f <- fnames[1]

        ## Need to call read.Agilent here to take care of
        ## layout problem
        
        if (DEBUG) print("call read.Agilent")
        gp <- read.Agilent(fnames = fnames[1], path = path, DEBUG=DEBUG)
        numrow <- nrow(gp@maRf)
        numcol <- length(fnames)
        
        Rf <- Gf <- Rb <- Gb <- weight <- matrix(0,nrow=numrow, ncol=numcol)
        mlayout <- maLayout(gp)
        gnames <- maGnames(gp)
        mlayout@maControls <- as.factor(maGenControls(gnames, controlcode = controlMatrix,
                                                      id=controlId))
 
        nb <- filenames <- c()
        
        # Call to slideQuality for each Agilent txt file

        if (DEBUG) print("Starting loop for boxplot")
        for (i in 1:length(fnames))
          {
            if (DEBUG) print(i)
            if (DEBUG) print("In the loop ")
            f <- fnames[i]
            gp <- readAgilent(fnames = f, path=path)
            ##Warning: Here, controlId must be = to a name in readAgilent
            ##ProbeName IS NOT in the list (ID is!)
            restmp <- slideQuality(gp, controlMatrix = controlMatrix, #controlId = controlId,
                                   DEBUG=DEBUG)
                        
            ###start plot
            Gf[,i] <- gp[["GfMedian"]]; Gb[,i] <-gp[["GbMedian"]] 
            Rf[,i] <- gp[["RfMedian"]]; Rb[,i] <- gp[["RbMedian"]]
            weight[,i] <- gp[["Flags"]]
            filenames <- c(filenames, gp[["File"]])
        
            #qualBoxplot
            if (DEBUG) print("Ploting")
            setwd(resdir)
            plotname <- paste("qualPlot",unlist(strsplit(f, ".txt")), plotdef$suffix,sep=".")
            plotdef <- c(plotdef, list(main=paste(f, ": Quantitative Diagnostic Plots")))

            ##do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)),
             ##                        formals(args(dev))))
            
            if(plotdef$suffix != "ps")
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)), formals(args(dev))))
            else
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)), formals(args(dev))))

            
            par(mar=c(3,14,2,6))
            nbtmp <- qualBoxplot(restmp, reference=reference, organism=organism, DEBUG=DEBUG)
            dev.off()
            setwd(curdir)

            if (DEBUG) print("Done ploting")
            
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
        colnames(Gf) <- colnames(Gb) <- colnames(Rf) <- colnames(Rb) <- colnames(weight) <- filenames
        if (DEBUG) print("building mraw")
        mraw <- new("marrayRaw", maRf=Rf, maGf=Gf, maRb=Rb,
                    maGb=Gb, maNotes="", maLayout=mlayout,
                    maW=weight, maGnames=gnames)

        ##maQualityPlots
        setwd(resdir)
        print("Starting maQualityPlots")

        defs <- list(norm="l")
        norm.defs <- maDotsDefaults(opt, defs)     
        
        do.call(maQualityPlots, c(list(mrawObj=mraw, controlId=controlId,
                                         DEBUG=DEBUG, dev=dev),
                                    norm.defs))
                
        #get diagnostic plots names
        tmpname <- sub(".txt", "",colnames(mraw@maGf))
        pn <- paste("diagPlot", tmpname, sep=".")

        dirfiles <- dir(".")
        
        for(i in 1:length(pn))
          Dp <- c(Dp, dirfiles[grep(pn[i], dirfiles)[1]])
        
        if (DEBUG) print("After for loop")
        if (DEBUG) print(nb)
        quality2HTML(fnames=fnames,path=resdir, QCplot=QCp, DiagPlot=Dp,nbBelow=nb)
        
        colnames(quality) <- fnames
        rownames(quality) <- meas
        
        print("agQuality done")
        
       ####### Results
        
        if (output)
          {
            print("Printing results to file")
            write.table(quality, "quality.txt",sep="\t", col.names=NA)
            do.call(outputNormData, c(list(mraw=mraw, val=c("maM", "maA")), norm.defs))            
          }

        setwd(curdir)
        return(list(mraw=mraw, quality=quality))
      }

   else {

     ##############################################################
     ## COMPBOXPLOT = FALSE
     ##############################################################
     
     curdir <- getwd()
     if (!file.exists(resdir))
       dir.create(resdir)
     if (DEBUG) print(getwd())

     
     mraw <- read.Agilent(fnames, path)
     colnames(maGf(mraw)) <- fnames
     print("Starting maQualityPlots")
     defs <- list(norm="l")
     norm.defs <- maDotsDefaults(opt, defs)     
     setwd(resdir)
     do.call(maQualityPlots, c(list(mrawObj=mraw, controlId=controlId, DEBUG=DEBUG, dev=dev),
                                 norm.defs))

     print("agQuality done")
     
     if (output)
       {
         print("Printing results to file")
         do.call(outputNormData, c(list(mraw=mraw, val=c("maM", "maA")), norm.defs))
       }

     setwd(curdir)
     return(list(mraw=mraw))     
   }
  }



