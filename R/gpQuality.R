#####################################################
## Set of function written to estimate
## quality of microarray slides
##
## Author: Agnes Paquet
## Modified: 04/16/2004
##           09/21/2004
##           09/24/2004 
## source("~/Projects/madman/Rpacks/arrayQuality/R/readGPR.R")

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
      fnames <- dir(path, pattern = ".*\\.gpr$")
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
      {
        normCoef1 <- gsub("\"", "",
                          strsplit(tmp[grep("NormalizationFactor", tmp)],
                                   split = "=")[[1]][2])
        normCoef <- gsub("\t", ",", normCoef1)
      }

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
## Must be run on the same gpr files
gpQuality <- function(fnames = NULL, path = ".",
                      organism=c("Mm", "Hs"),
                      compBoxplot = TRUE,
                      reference=NULL,
                      ##scalingTable=NULL,
                      controlMatrix = controlCode,
                      controlId = c("ID", "Name"),
                      output=FALSE,
                      resdir=".",
                      dev="png", #set default to be png 
                      DEBUG = FALSE,...)
{

    print("Starting gpQuality...")
    
    # Check input arguments
    if (missing(path) | is.null(path))
      path <- "."
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = ".*\\.gpr$")

    organism <- organism[1]
    controlId <- controlId[1]
    if (DEBUG) print(controlId)

    opt <- list(...)

    if (DEBUG) print(controlMatrix)
        
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
    
        # Prepares results
        quality <- NULL
        #tmp <- NULL
        QCp <- c()
        Dp <- c()
        
        curdir <- getwd()
        if (!file.exists(resdir))
          dir.create(resdir)
        if (DEBUG) print(getwd())
        
        #Allocation of matrix for marrayraw object

        if(DEBUG) print(path)
        if(DEBUG) print(resdir)

        #f <- fnames[1]
        #gp <- readGPR(fnames=f, path=path)
        #nrow <- length(gp[["RfMedian"]])
        #ncol <- length(fnames)
        #if(DEBUG) print("creating layout")
        #mlayout <- maCompLayout(as.matrix(cbind(gp[["Block"]],
        #                                        gp[["Row"]], gp[["Column"]])))
        #tmp <- new("marrayInfo", maInfo=data.frame(gp[["Name"]],gp[["ID"]]))
        #mlayout@maControls <- as.factor(maGenControls(tmp, controlcode = controlMatrix,
         #                                             id=controlId)) 
        #rm(f, gp)

        ## Read layout and allocation of matrix for marrayRaw

        if (DEBUG) print("call read.Galfile")
        galf <- read.Galfile(galfile = fnames[1], path = path)
        mlayout <- galf$layout

        gnames <- galf$gnames
        mlayout@maControls <- as.factor(maGenControls(gnames, controlcode = controlMatrix,
                                                      id=controlId))
        galOrder <- galf$neworder
        nb <- filenames <- c()
        numrow <- maNspots(galf$layout)
        numcol <- length(fnames)
        Rf <- Gf <- Rb <- Gb <- weight <- matrix(0,nrow=numrow, ncol=numcol)
        ##numrow <- nrow(gp@maRf)
        ##numcol <- length(fnames)
        
        ## Call to slideQuality for each gpr file

        for (i in 1:length(fnames))
          {
            if (DEBUG) print("In the loop ")
            f <- fnames[i]
            gp <- readGPR(fnames = f, path=path)
            restmp <- slideQuality(gp, controlMatrix = controlMatrix, DEBUG=DEBUG)
            
            ###start plot
            Gf[,i] <- gp[["GfMedian"]]; Gb[,i] <-gp[["GbMedian"]] 
            Rf[,i] <- gp[["RfMedian"]]; Rb[,i] <- gp[["RbMedian"]]
            weight[,i] <- gp[["Flags"]]
            filenames <- c(filenames, gp[["File"]])
        
            #qualBoxplot
            setwd(resdir)
            plotname <- paste("qualPlot",unlist(strsplit(f, ".gpr")), dev,sep=".")
            plotdef <- c(plotdef, list(main=paste(f, ": Quantitative Diagnostic Plots")))

            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)),
                                     formals(args(dev))))

            ## was: do.call(dev, maDotsDefaults(opt, c(list(filename=plotname), plotdef$dev)))
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
        colnames(Gf) <- colnames(Gb) <- colnames(Rf) <- colnames(Rb) <- colnames(weight) <- filenames
        if (DEBUG) print("building mraw")

        ## Warning: all vectors must be reordered according to
        mraw <- new("marrayRaw", maRf=Rf, maGf=Gf, maRb=Rb,
                    maGb=Gb, maNotes="", maLayout=mlayout,
                    maW=weight, maGnames=gnames)

        mraw <- mraw[galf$neworder,]
        mraw@maLayout@maSub <- mlayout@maSub
 
        #maQualityPlots
        setwd(resdir)
        print("Starting maQualityPlots")
        
        defs <- list(norm="p")
        norm.defs <- maDotsDefaults(opt, defs)     
        
        do.call("maQualityPlots", c(list(mrawObj=mraw, controlId=controlId, DEBUG=DEBUG),
                                    norm.defs))
        
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

     if (DEBUG) print("In compBoxplot=FALSE")
     curdir <- getwd()
     if (!file.exists(resdir))
       dir.create(resdir)
     if (DEBUG) print(getwd())

     
     mraw <- read.GenePix(fnames, path)
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
         #colnames(mraw@maGnames@maInfo) <- c("Name", "ID")
         do.call("outputNormData", c(list(mraw=mraw), norm.defs))
       }

     setwd(curdir)
     return(list(mraw=mraw))     
   }
  }


