##############################################################
## Wrapper function for MEEBO QC
## Date: 07/03/2006
## Author: Agnes
##############################################################

## To do:
## 1) Read in file + layout
## 2) Normalization
## 3) Filtering
## 4) Plot: Mismatch, tiling, BE, spike-ins

##source("C:/MyDoc/Projects/Rcode/arrayQualityTmp/R/meeboQuality.R")
##source("C:/MyDoc/Projects/madman/Rpacks/arrayQuality/R/meeboQuality.R")
##meeboQuality(spikeTypeFile="DopingTypeFile2.txt",DEBUG=TRUE, norm="n")

##load("C:/MyDoc/Projects/Statistics/MEEBO/Ash/MEEBOtiling.RData")

meeboQuality <- function(fnames = NULL, path = ".",
                         galfile=NULL,
                         source="genepix.median",
                         other.columns = c("Flags"),
                         controlMatrix = MeeboSpotTypes, 
                         controlId = c("ID", "Name"),
                         DOPING=TRUE,
                         meeboSetQC=TRUE,
                         SpotTypeFile=NULL,
                         SpikeTypeFile=NULL,
                         cy3col="CY3.ng._MjDC_V1.7",
                         cy5col="CY5.ng._MjDC_V1.7",
                         id="SeqID",
                         namecol=c("Symbol","Name"),
                         annot=NULL,
                         bgMethod="none",
                         normMethod="p",
                         ##FILTER=FALSE,
                         diagnosticPlot=TRUE,
                         output=TRUE,
                         resdir = ".", 
                         dev = "png",
                         organism="Mm",
                         DEBUG = FALSE, ...)
  {
    print("Staring meeboQuality")
    require(MEEBOdata)
    ## Check input arguments
    if (missing(path) | is.null(path))
      path <- "."
    
    if (missing(fnames) | is.null(fnames))
      fnames <- dir(path, pattern = ".*\\.gpr$")

    if(!DOPING)
      print("Spike-in QC plots won't be generated")
    else
      {
        if (missing(SpikeTypeFile) | is.null(SpikeTypeFile))
          {
            print("Missing doping control information")
            print("Spike-in QC plots will not be generated")
            DOPING <- FALSE
          }
      }
       
    if (DEBUG) print(paste("DOPING=",DOPING,sep=""))
    if(!meeboSetQC)
      {
        print("MEEBO set quality control plots will not be generated")
      }
    if(DEBUG) print(paste("meeboSetQC = ",meeboSetQC,sep=""))
    if(DEBUG) print(paste("diagnosticPlot = ",diagnosticPlot,sep=""))
    
    controlId <- controlId[1]
    if (DEBUG) print(controlId) 
    opt <- list(...)
 
    ## Load required R objects
    if(is.null(annot))
      {
        if(!("MEEBOset" %in% ls(1)))
          data(MEEBOset)
        annot=MEEBOset
      }

    if(DOPING)
      {
        if(DEBUG) print("Parsing Spike Type File")
        spike.defs <- list(file=SpikeTypeFile,path=path,cy5col=cy5col,cy3col=cy3col)
        spike.args <- maDotsMatch(maDotsMatch(opt, spike.defs), formals(args("readSpikeTypes")))
        SpikeList <- do.call("readSpikeTypes",spike.args)
        spike.args <- maDotsMatch(maDotsMatch(opt, spike.defs), formals(args("readSpotTypes")))
        MeeboSpikeTypes <- do.call(readSpotTypes,spike.args)
        if(length(namecol)>1)
          namecol=namecol[1]

      }

    ## Reading data in
    if (!is.null(galfile))
      {
        if(DEBUG) print("Reading Galfile ... ")
        gal.defs <- list(galfile = galfile, path=path, fill=TRUE)        
        gal.args <- maDotsMatch(maDotsMatch(opt, gal.defs), formals(args("readGAL")))
        gal.args <- maDotsMatch(maDotsDefaults(opt, gal.defs), formals(args("readGAL")))
        gal <- do.call("readGAL", gal.args)

        gpr.defs <- list(files=fnames, path=path,source=source,other.columns=other.columns)
        gpr.args <- maDotsMatch(maDotsMatch(opt,gpr.defs),formals(args(read.maimages)))
        
        RG <- do.call("read.maimages",gpr.args)
        RG$genes <- gal
        RG$printer <- getLayout(RG$genes)
      }
    else
      {
        if (DEBUG) print("Reading .gpr ...")
        gpr.defs <- list(files=fnames, path=path,source=source,other.columns=other.columns)
        gpr.args <- maDotsMatch(maDotsMatch(opt,gpr.defs),formals(args(read.maimages)))     
        RG <- do.call("read.maimages",gpr.args)
        RG$printer <- getLayout(RG$genes)
      }

    ## Set control status
    if(DEBUG) print("Setting up control status")
    if(!is.null(SpotTypeFile))
      {
        spot.defs <- list(file=SpotTypeFile)
        spot.args <- maDotsMatch(maDotsMatch(opt,spot.defs),formals(args(read.maimages)))     
        controlMatrix <- do.call("readSpotTypes",spot.args)
      }
    RG$genes$Status <- controlStatus(controlMatrix,RG,regexpcol=controlId,verbose=FALSE)
    if(DEBUG) print(table(RG$genes$Status))

    ## Don't forget to set up rownames==ID
    rownames(RG$R) <- rownames(RG$G) <- RG$genes[,controlId]

    ##############################################################
    ## Pre-processing: Bg subtraction, normalization

    ## Bg subtraction
    bg.defs <- list(method=bgMethod, RG=RG)
    bg.args <- maDotsMatch(maDotsMatch(opt,bg.defs),formals(args(backgroundCorrect)))
    RGsub <- do.call("backgroundCorrect",bg.args)
    rownames(RGsub$R) <- rownames(RGsub$G) <- RG$genes[,controlId]

    ## Normalization
    norm.defs <- list(method=normMethod, object=RGsub)
    norm.args <- maDotsMatch(opt,norm.defs)
    MA <- do.call("normalizeWithinArrays",norm.args)
    ## need to set up rownames for MA$M and MA$A
    rownames(MA$A) <- rownames(MA$M) <- MA$genes[,controlId]
    
    ## Set up output device
    if (DEBUG) print("Setting up default output device")
    plotdef <- switch(dev,
                      "bmp" = list(dev=list(width=800, height=600, bg="white"), suffix="bmp"),
                      "jpeg" = list(dev=list(quality=100, width=800, height=600, bg="white"), suffix="jpeg"),
                      "postscript" = list(dev=list( bg="white"), suffix="ps"),
                      "png" =  list(dev=list(width=800, height=600, bg="white"), suffix="png"),
                      list(dev=list(width=800, height=600,bg="white"), suffix="png"),
                    )
    
    if(!is.element(dev, c("bmp", "jpeg","png","postscript")))
      {
        print("Format error, format will be set to PNG")
        dev = "png"
      }

    ## Switch do result directory if needed
    if (DEBUG) print("Switching to result directory")
    if (!file.exists(resdir))
      dir.create(resdir)
    curdir <- getwd()
    setwd(resdir)

    ###############################################################
    ## Generate diagnostic plots
    if(diagnosticPlot)
      {
        print("Starting meeboQualityPlots")
        qp.defs <-list(mrawObj=RG,controlId=controlId,DEBUG=DEBUG,
                       dev=dev,norm=normMethod, organism=organism)
        qp.args <- maDotsDefaults(opt, qp.defs)
        do.call("meeboQualityPlots", qp.args)
      }

    ###############################################################
    ## Write normalized data to file

    if (output)
      {
        print("Printing results to file")
        ##colnames(mraw@maGnames@maInfo) <- c("Name", "ID")
        ## do.call("outputNormData", c(list(mraw=RG, val=val), norm.defs))
        tmp <- c()
        cnames <- c()
        for(i in 1:ncol(MA$A))
          {
            tmp <- cbind(tmp, MA$M[,i],MA$A[,i])
            cnames <- c(cnames, paste("M",colnames(MA$M)[i],sep="_"),
                        paste("A", colnames(MA$A)[i],sep="_"))
          }
        colnames(tmp) <- cnames
        write.table(cbind(Block=MA$genes$Block,
                          Column=MA$genes$Column, Row=MA$genes$Row,
                          Name=MA$genes$Name, ID=MA$genes$ID,
                          tmp),file="NormalizedData.xls",
                    sep="\t",quote=FALSE,row.names=FALSE)
        rm(cnames,tmp)
      }
    
    ##########################################################
    ## Specific individual control plots
    if(meeboSetQC)
      {
        if(DEBUG) print("Starting specific control plots...")
    
        ## Loading controls
        if(DEBUG) print("Loading mismatch control and binding energy info")

        if(!("MEEBOtilingres" %in% ls(1)))
          data(MEEBOtilingres)
        
        if(!("MEEBOctrl" %in% ls(1)))
          data(MEEBOctrl)
                
        for(i in 1:dim(MA)[2])
          {
            if (DEBUG) print(i)
            
            ## Mismatch plot
            if(DEBUG) print("qpMisMatchPlot")
            plotname <- paste("Mismatch",colnames(MA$A)[i],plotdef$suffix,sep=".")
            if(plotdef$suffix != "ps")
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname,width=1600,height=1200),
                                                             plotdef$dev)),
                                       formals(args(dev))))
            else
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                       formals(args(dev))))
            par(mfrow=c(3,4), cex=1.5, cex.axis=0.8, cex.main=1.2,mar=c(2,2,2,2),oma=c(3,3,3,3))
            qpMisMatchPlot(MA[,i],ctrllist=MEEBOctrl,organism=organism,DEBUG=DEBUG)
            title(paste("Mismatch plot for", colnames(MA$A)[i],sep=" "),outer=TRUE,cex=1.2)
            dev.off()
            
            ## Binding Energy ------------------------ 
            
            if(DEBUG) print("qpBEPlot")
            plotname <- paste("BindingEnergy",colnames(MA$A)[i],plotdef$suffix,sep=".")
            if(plotdef$suffix != "ps")
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname),
                                                             plotdef$dev)),
                                       formals(args(dev))))
            else
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                       formals(args(dev))))
            par(mar=c(5,5,5,5))
            qpBEplot.linear(RGsub[,i], ctrllist=MEEBOctrl,DEBUG=DEBUG)
            dev.off()                           
            ##-------------------------------------------------------------
            
            ## Titling controls
            if(DEBUG) print("qpTiling")
            plotname <- paste("Tiling",colnames(MA$A)[i],plotdef$suffix,sep=".")
            if(plotdef$suffix != "ps")
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname,width=1500,height=1200),
                                                             plotdef$dev)),formals(args(dev))))
            else
              do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                       formals(args(dev))))
            par(mfcol=c(3,4), cex=1.3, cex.axis=1, mar=c(3,4,3,2),oma=c(2,2,4,2))
            qpTiling(RGsub[,i], tilingres=MEEBOtilingres,organism=organism,DEBUG=DEBUG)
            title(paste("Raw signal vs. 3' distance:",colnames(MA$A)[i],sep=" "),
                  cex=2,outer=TRUE)
            dev.off()
            
            if(DEBUG) print("End of Meebo set QC")
          }
      }

        ##-------------------------------------------------------------       
        ## Spike-ins
        
        if(DOPING)
          {
            for(i in 1:dim(MA)[2])
              {
                
                if(DEBUG) print("Starting QC plots for doping controls")
                ## Cy3 vs Cy5 raw signal
                if(DEBUG) print("Spike.Cy5vsCy3")
                spike.defs <- list(rawobj=RGsub[,i],spikesTypes=MeeboSpikeTypes,id=id,annot=annot,
                                   gnames=RGsub$genes[,controlId],cy5col=cy5col,cy3col=cy3col)
                spike.args <- maDotsDefaults(maDotsMatch(opt, spike.defs),formals(args(Spike.Cy5vsCy3)))

                plotname <- paste("Spike.Cy5vsCy3",colnames(MA$A)[i],plotdef$suffix,sep=".")
                if(plotdef$suffix != "ps")
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=plotname,width=800,height=800),
                                                                 plotdef$dev)),formals(args(dev))))
                else
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                           formals(args(dev))))
                
                if(DEBUG) print(names(spike.args))
                par(mar=c(5,5,5,5))
                do.call("Spike.Cy5vsCy3",spike.args)
                dev.off()
                
                ## MMplot obserevd vs. expected
                spike.defs <- list(mval=MA$M[,i],spikeList=SpikeList,id=id,annot=annot,
                                   gnames=MA$genes[,controlId],namecol=namecol,cy5col=cy5col,cy3col=cy3col)
                spike.args <- maDotsMatch(maDotsDefaults(opt,spike.defs),formals(args(Spike.MMplot)))
                
                ## Will generate 2 plots/spike type
                ntypes <- length(SpikeList)
                
                ## Spike.MMplot
                if(DEBUG) print("Spike.MMplot")
                plotname <- paste("Spike.MMplot",colnames(MA$A)[i],plotdef$suffix,sep=".")
                if(plotdef$suffix != "ps")
                  do.call(dev, maDotsMatch(maDotsDefaults(opt,c(list(filename=plotname,width=800,height=500*ntypes),
                                                                plotdef$dev)),formals(args(dev))))
                else
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                           formals(args(dev))))
                par(mfrow=c(ntypes,2),cex.main=2,cex.lab=1.5,oma=c(7,3,3,3))
                do.call("Spike.MMplot",spike.args)
                title(paste("Observed and expected ratios for",colnames(MA$A)[i],sep=" "),cex=2.5,outer=TRUE)
                dev.off()
                
                ##Observed ratio/expected ratio scattered
                
                if(DEBUG) print("Spike.MM.Scatter")
                plotname <- paste("Spike.MM.Scatter",colnames(MA$A)[i],plotdef$suffix,sep=".")
                if(plotdef$suffix != "ps")
                  do.call(dev, maDotsMatch(maDotsDefaults(opt,c(list(filename=plotname,width=600,height=500*ntypes),
                                                                plotdef$dev)),formals(args(dev))))
                else
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                           formals(args(dev))))
                par(mfrow=c(ntypes,1),oma=c(1,3,3,3),mar=c(7,3,3,3))
                spike.args <-  maDotsMatch(maDotsDefaults(opt,spike.defs),formals(args(Spike.MM.Scatter)))
                do.call("Spike.MM.Scatter",spike.args)
                title(paste("Observed and expected ratios for",colnames(MA$A)[i],sep=" "),cex=2.5,outer=TRUE)
                dev.off()
                
                ## Spike-in sensitivity
                
                if(DEBUG) print("Spike Sensitivity")            
                spike.defs <- list(rawobj=RGsub[,i],spikeList=SpikeList,id=id,annot=annot,ctllist=MEEBOctrl,
                                   gnames=RGsub$genes[,controlId],namecol=namecol,cy5col=cy5col,cy3col=cy3col)
                spike.args <- maDotsMatch(maDotsDefaults(opt,spike.defs),formals(args(Spike.Sensitivity)))
                
                plotname <- paste("Spike.Sensitivity",colnames(MA$A)[i],plotdef$suffix,sep=".")
                if(plotdef$suffix != "ps")
                  do.call(dev, maDotsMatch(maDotsDefaults(opt,c(list(filename=plotname,width=800,height=400*ntypes),
                                                                plotdef$dev)),formals(args(dev))))
                else
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                           formals(args(dev))))
                
                par(mfrow=c(ntypes,2),oma=c(3,3,3,3),cex.lab=1.2,cex.axis=1.2)
                do.call("Spike.Sensitivity",spike.args)
                title(paste("Sensitivity", colnames(MA$A)[i],sep=" "),outer=TRUE,cex=2.5)
                dev.off()
                
                if(DEBUG) print("Spike Sensitivity, individual spikes plots")            
                
                ##Spike.Individual.Sensitivity
                
                plotname <- paste("Spike.Sensitivity.Indi",colnames(MA$A)[i],plotdef$suffix,sep=".")
                if(plotdef$suffix != "ps")
                  do.call(dev, maDotsMatch(maDotsDefaults(opt,c(list(filename=plotname,width=1000,height=400*ntypes),
                                                                plotdef$dev)),formals(args(dev))))
                else
                  do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=plotname), plotdef$dev)),
                                           formals(args(dev))))
                par(mfrow=c(ntypes,2),oma=c(3,3,3,3),cex.lab=1.2,cex.axis=1.2)
                spike.defs <- list(rawobj=RGsub[,i],spikeList=SpikeList,id=id,annot=annot,ctllist=MEEBOctrl,
                                   gnames=RGsub$genes[,controlId],namecol=namecol,cy5col=cy5col,cy3col=cy3col)
                spike.args <- maDotsMatch(maDotsDefaults(opt,spike.defs),
                                          formals(args(Spike.Individual.Sensitivity)))
                do.call("Spike.Individual.Sensitivity",spike.args)
                title(paste("Sensitivity", colnames(MA$A)[i],sep=" "),outer=TRUE,cex=2.5)
                dev.off()
              }
          }

     ## Go back to previous directory
    
    if(DEBUG) print("End of meeboQuality")
    setwd(curdir)
    return(MA)
  }


