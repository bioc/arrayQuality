###################################################
## READ  Spot file .txt file, returns a list of measures
## extracted from the Agilent file
## source("~/Projects/madman/Rpacks/arrayQuality/R/spotQuality.R")
##
###################################################

## No background substraction
## No normalization
## Suppose that fnames is a .spot file.
## Returns a list of data needed for quality
## Combines read.marrayRaw, maQualityMain, spotTools
## One slide only!!!

## returns a list of vector containing info from grp file column
readSpot <- function (fnames = NULL, path= ".", galfile = NULL, DEBUG=FALSE, skip = 0,
                     sep ="\t", quote= "\"", controlId = "ID", ...)
  {
    if (DEBUG) print("Starting readSpot")
    if (DEBUG) print(path)
    controlId <- controlId[1]
    # Test if input data is OK
##    if (is.null(path))  path <- "."

    if(!is.null(galfile))
      galf <- file.path(path, galfile)
    else
      {
        galf <- dir(pattern = "*\\.gal$")[1]
        if(is.na(galf))
          stop("Missing galfile information.")
        else
          cat("User did not specify a galfile name, reading in ", galf, " from current directory. \n")
      }

    opt <- list(...)
    defs <- list(galfile=galf, path = path, sep=sep, skip=skip,
                  quote = quote) 
    read.Galfile.defs <- maDotsMatch(maDotsDefaults(opt, defs),
                                     formals(read.Galfile))
    if (DEBUG) print("Calling read.Galfile")
    #read.Galfile.defs <- maDotsDefaults(opt, defs)
    galdata <- do.call("read.Galfile", read.Galfile.defs)
    
    #galdata <- read.Galfile(galf)

    if (missing(fnames) | is.null(fnames)) 
      fnames <- dir(path, pattern = "*\\.spot$")
    fullfnames <- file.path(path, fnames)
    f <- fullfnames[1]
    opt <- list(...)
    
    # Search in the Spot file where are the colums starting
    if (DEBUG) print("Estimating number of lines to skip")
    y <- readLines(fullfnames[1], n = 100)
    skip <- grep("Gmedian", y)[1] - 1
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

    if (DEBUG) print(paste("Reading", f))
    h <- strsplit(readLines(f, n = skip + 1), split = sep)
    h <- as.list(unlist(h[[length(h)]]))
    names(h) <- gsub("\"", "", unlist(h))
#    dat <- scan(f, quiet = TRUE, what = h, sep = sep, skip = skip + 
 #               1, quote = quote, ...)
    scan.defs <- list(file=f,quiet = TRUE, what = h, sep = sep, skip = skip + 
                1, quote = quote)
    scan.args <- maDotsMatch(maDotsMatch(opt, scan.defs),
                             formals(args("scan")))
    
    dat <- do.call("scan", scan.args[-grep("nmax", names(scan.args))])
    gc <- as.numeric(dat[["grid.r"]])
    gr <- as.numeric(dat[["grid.c"]])
    numcol <- max(gc)
    block <- cbind(block,  (gr-1) * numcol + gc)
    column <- cbind(column, as.numeric(dat[["spot.c"]]))
    row <- cbind(row, as.numeric(dat[["spot.r"]]))
    name <- cbind(name, galdata$gnames@maInfo[["Name"]])
    id <- cbind(id, galdata$gnames@maInfo[[controlId]])
    if (DEBUG)
      {
        print("cchecking id")
        print(class(id))
        print(length(id))
      }
    Gfmedian <- cbind(Gfmedian, as.numeric(dat[["Gmedian"]]))
    Gfmean <- cbind(Gfmean, as.numeric(dat[["Gmean"]]))        
    Gfsd <- cbind(Gfsd, as.numeric(dat[["GIQR"]]))
    Gbmedian <- cbind(Gbmedian, as.numeric(dat[["bgGmed"]]))
    Gbmean <- cbind(Gbmean, as.numeric(dat[["bgGmean"]]))        
    Gbsd <- cbind(Gbsd, as.numeric(dat[["bgGSD"]]))
    Rfmedian <- cbind(Rfmedian, as.numeric(dat[["Rmedian"]]))
    Rfmean <- cbind(Rfmean, as.numeric(dat[["Rmean"]]))
    Rfsd <- cbind(Rfsd, as.numeric(dat[["RIQR"]]))
    Rbmedian <- cbind(Rbmedian, as.numeric(dat[["bgRmed"]]))
    Rbmean <- cbind(Rbmean, as.numeric(dat[["bgRmean"]]))
    Rbsd <- cbind(Rbsd, as.numeric(dat[["bgRSD"]]))        
##    Gfsat <- cbind(Gfsat, as.numeric(dat[["gNumSatPix"]]))
##    Rfsat <- cbind(Rfsat, as.numeric(dat[["rNumSatPix"]]))
##    spotDia <- cbind(spotDia, as.numeric(dat[["Dia."]]))  Not used
    spotFArea <- cbind(spotFArea, as.numeric(dat[["area"]]))
##    spotBArea <- cbind(spotBArea, as.numeric(dat[["gBGNumPix"]]))
    Flags <- cbind(Flags, as.numeric(dat[["circularity"]]) < 0.5)
    gprData <- list(File = fnames[1], Date = NA, PmtR = NA, PmtG = NA,
                    Normalization = NA, NormCoefficient = NA,
                    Name = name, ID = id,
                    Block = block, Column = column, Row = row,
                    GfMedian = Gfmedian, GfMean = Gfmean, GfSD = Gfsd,
                    GbMedian = Gbmedian, GbMean = Gbmean, GbSD = Gbsd,
                    RfMedian = Rfmedian, RfMean = Rfmean, RfSD = Rfsd,
                    RbMedian = Rbmedian, RbMean = Rbmean, RbSD = Rbsd,
                    GfSaturation = NA, RfSaturation = NA,
                    SpotDiameter = NA, spotArea = spotFArea,
                    bgArea = NA, Flags = Flags)
    
    # Result
    return(gprData)
  }



###################################################
## Given a spot files
## Computes needed statistics to assess quality
###################################################

## Takes all .spot file into acccount
## Plot quality boxplot and diagnostic plots
## creates html report
## if output: writes quality and normalized data to file
## return quality measure and marrayRaw object in a list

#############
## Example:
## test <- gpQuality(path="C:/Mydoc/Projects/quality/DemoA/", resdir="QualPlot")

## Reference = output from globalQuality
## ScalingTable: output from qualRef
## Must be run on the same spot files

###########################################################################
##
##  Move controlCode and maGenControls from marrayTools to marrayClasses
##  May 7, 2003
##  Sept 21, 2004
## TO USE:
## 1. Create Reference slides
##    reference <- globalQuality(fnames, inputsource = "readSpot")
## 2. Run : some alternative
##        -- spotQuality(fnames, reference = reference)  ## Use what you have created
##        -- spotQuality(fnames, reference = reference, compBoxplot=FALSE)
##                                    ## Only generate qualitative plots
##        -- spotQuality(fnames)  ## Use exisiting reference info
###########################################################################

spotQuality <- function(fnames = NULL, path = ".", galfile = NULL,
                      organism=c("Mm", "Hs"),
                      compBoxplot = TRUE,
                      reference=NULL,
                      #scalingTable=NULL,
                      controlMatrix = controlCode,
                      controlId = c("ID"),
                      output=FALSE,
                      resdir=".",
                      dev="png", #set default to be png 
                      DEBUG = FALSE,...)
{

    print("Starting spotQuality...")
    
    # Check input arguments
    if (missing(path) | is.null(path))  path <- "."
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = "*\\.spot$")

    organism <- organism[1]
    controlId <- controlId[1]
    if (DEBUG) print(controlId)
    if (DEBUG) print(controlMatrix)
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

    ##############################################################
    ## COMPBOXPLOT = TRUE
    ##############################################################
    
    if (compBoxplot)
      {
        if(DEBUG) print("Starting compBoxplot")

        ## Prepares results
        quality <- NULL
        tmp <- NULL
        QCp <- c()
        Dp <- c()
        
        curdir <- getwd()
        if (!file.exists(resdir))
          dir.create(resdir)
        if (DEBUG) print(getwd())
        
        #Allocation of matrix for marrayraw object
        if(!is.null(galfile))
          galf <- galfile
        else
          {
            galf <- dir(path=path, pattern = "*\\.gal$")[1]
            if(is.na(galf))
              stop("Missing galfile information.")
            else
              cat("No specific galfile names, reading in ", galf, " from ", path,  "directory. \n")
          }
        
        read.Galfile.defs <- maDotsMatch(maDotsDefaults(opt, list(galfile=galf)), formals(read.Galfile))
        galdata <- do.call("read.Galfile", read.Galfile.defs)
 
        #galdata <- read.Galfile(galf, path=path)
        
        if(DEBUG) print(path)
        if(DEBUG) print(resdir)
        
        f <- fnames[1]
        gp <- readSpot(fnames=f, path=path, galfile=galfile, controlId=controlId,...)
        numrow <- length(gp[["GfMedian"]])
        numcol <- length(fnames)

        
        if(DEBUG) print("creating layout")
        mlayout <- galdata$layout
        gnames <- galdata$gnames
        ##read controls
        mlayout@maControls <- as.factor(maGenControls(gnames, controlcode = controlMatrix,
                                                      id=controlId))
        
        Rf <- Gf <- Rb <- Gb <- weight <- matrix(0,nrow=numrow, ncol=numcol)
        filenames <- c()
        nb <- c()
        rm(f, gp)
        
        # Call to slideQuality for each spot file
        if (DEBUG) print("Before loop")
        for (i in 1:length(fnames))
          {
            if (DEBUG) print("In the loop ")
            f <- fnames[i]
            gp <- readSpot(fnames = f, path=path, galfile=galfile,controlId = controlId, ...)
            restmp <- slideQuality(gp, controlMatrix = controlMatrix, DEBUG=DEBUG) #controlId = "ID"
                        
            ###start plot
            Gf[,i] <- gp[["GfMedian"]]; Gb[,i] <-gp[["GbMedian"]] 
            Rf[,i] <- gp[["RfMedian"]]; Rb[,i] <- gp[["RbMedian"]]
            weight[,i] <- gp[["Flags"]]
            filenames <- c(filenames, gp[["File"]])
        
            #qualBoxplot
            setwd(resdir)
            plotname <- paste("qualPlot",unlist(strsplit(f, ".spot")), dev,sep=".")
            plotdef <- c(plotdef, list(main=paste(f, ": Quantitative Diagnostic Plots")))
            
            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)),
                                     formals(args(dev))))
            par(mar=c(3,14,2,6))
            nbtmp <- qualBoxplot(restmp, reference=reference, organism=organism)
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
                    maW=weight, maGnames=gnames)
        #maQualityPlots
        setwd(resdir)
        print("Starting maQualityPlots")

        defs <- list(norm="p")
        norm.defs <- maDotsDefaults(opt, defs)     

        do.call("maQualityPlots", c(list(mrawObj=mraw, controlId=controlId, DEBUG=DEBUG),
                                    norm.defs))
        
        #get diagnostic plots names
        tmpname <- sub(".spot", "",colnames(mraw@maGf))
        pn <- paste("diagPlot", tmpname, sep=".")

        dirfiles <- dir(".")
        
        for(i in 1:length(pn))
          Dp <- c(Dp, dirfiles[grep(pn[i], dirfiles)[1]])
        
        if (DEBUG) print("After for loop")
        if (DEBUG) print(nb)
        quality2HTML(fnames=fnames,path=resdir, QCplot=QCp, DiagPlot=Dp,nbBelow=nb)
        
        colnames(quality) <- fnames
        rownames(quality) <- meas
        
        print("SpotQuality done")
        
       ####### Results
        
        if (output)
          {
            print("Printing results to file")
            write.table(quality, "quality.txt",sep="\t", col.names=NA)
            #colnames(mraw@maGnames@maInfo) <- c("Name", "ID")
            do.call("outputNormData", c(list(mraw=mraw), norm.defs))
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

      if(!is.null(galfile))
        galf <- galfile
      else
        {
          galf <- dir(path=path, pattern = "*\\.gal$")[1]
          if(is.na(galf))
            stop("Missing galfile information.")
          else
            cat("No specific galfile names, reading in ", galf, " from ", path,  "directory. \n")
        }
        read.Galfile.defs <- maDotsMatch(maDotsDefaults(opt, list(galfile=galf)), formals(read.Galfile))
        galdata <- do.call("read.Galfile", read.Galfile.defs)
      
      #galdata <- read.Galfile(galf, path=path)
      
      if(DEBUG) print(path)
      if(DEBUG) print(resdir)
      
      if(DEBUG) print("creating layout")
      mlayout <- galdata$layout
      gnames <- galdata$gnames
      mlayout@maControls <- as.factor(maGenControls(gnames, controlcode = controlMatrix,
                                                    id=controlId))
     
      mraw <- read.Spot(fnames, path, layout=mlayout, gnames=gnames)
      colnames(maGf(mraw)) <- fnames
     
      print("Starting maQualityPlots")
      defs <- list(norm="p")
      norm.defs <- maDotsDefaults(opt, defs)     

      setwd(resdir)
      do.call("maQualityPlots", c(list(mrawObj=mraw, controlId=controlId, DEBUG=DEBUG),
                                  norm.defs))
      
      if (output)
        {
          print("Printing results to file")
         # colnames(mraw@maGnames@maInfo) <- c("Name", "ID")
          do.call("outputNormData", c(list(mraw=mraw), norm.defs))
        }

      setwd(curdir)
      return(list(mraw=mraw))     
    }
  }

