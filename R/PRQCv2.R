###################################################### 
## Print-run QC
## March 12, 2004
##
## Two functions: PRvQCHyb
####################################################### 

## PRvQCHyb(fnames="9mm.186.gpr", prname="9Mm", DEBUG=TRUE)
## PRvQCHyb(path="9Mm", prname="9Mm", DEBUG=TRUE)


PRvQCHyb<-  function(fnames,
                     path=".",
                     dev = "png",  #set default to be png
                     DEBUG=FALSE,
                     prargs=NULL,
                     samepr=TRUE,
                     prname="xMm",
                     save = TRUE,
                     col,
                     ...) 
{
  
  ## Setting defaults
  require(marray)
  require(limma)

  print("Starting PRvQCHyb")
  
  defs <- list(fill = TRUE, quote = "\"", check.names=FALSE, as.is=TRUE,
               name.Gf = "F532 Median", name.Rf = "F635 Median",
               name.Gb="B532 Median", name.Rb="B635 Median",
               name.W="Flags", comment.char="", sep="\t", header=TRUE)
  opt <- list(...)
  args <- maDotsMatch(opt, defs)

  ## Getting File names
 if (DEBUG) print("Getting File Names")
  if (missing(fnames))
    fnames <- file.path(path, dir(path, pattern = "*\\.gpr$"))
  
  if (is.null(path)) {
    fnames <- fnames
    resdir <- paste(prname, "PRQC", sep="")
  }else {
    resdir <- file.path(path, paste(prname,"PRQC", sep=""))    
  }

  if(!file.exists(resdir))
    dir.create(resdir)
  if(DEBUG) print(resdir)
  
  for(f in fnames)
    {
      
      ###################
      ## Set up output name
      if (DEBUG) print("Name the output file")
      tmp <- unlist(strsplit(f, "/"))
      fnopath <- tmp[length(tmp)]  ## filename without the path attached
      tmp <- unlist(strsplit(fnopath, "\\."))
      fstart <- paste(tmp[-length(tmp)], collapse=".")

      ###################
      ## marrayRaw
      read.args <- maDotsMatch(args, formals(args("read.GenePix")))
      read.args$fnames <- f
      #read.args$path <- NULL
      read.args$path <- path
      if(DEBUG) cat("Reading", read.args$file, "...\n")
      bgraw <- do.call("read.GenePix", read.args)
      GInfo <- maGeneTable(bgraw)
      mraw <- bgraw; mraw@maGb <- mraw@maRb <- matrix(0,0,0)

      ###################
      ## load existing DE genes

      if(length(grep("Mm",  prname)) == 1)
        {
          data(MmDEGenes)
          DEGenes <- MmDEGenes;
        }
             
      else
        {
          if(length(grep("Hs",  prname)) == 1)
            {
              stop("No previous DE results for Human yet")
              #data(HsDEGenes)
              #DEGenes <- HsDEGenes;
            }
          else
            {
              stop("No previous DE results for Human yet")
              print("No previous DE results")
            }

        }

      
      ###################
      ## Setting up output device
      ###################
      if (DEBUG) print("Name of output device")
      plotdef <- switch(dev,
                        "bmp" = list(dev=list(width=1600, height=1000, bg="white"), suffix="bmp"),
                        "jpeg" = list(dev=list(quality=100, width=1600, height=1000, bg="white"), suffix="jpeg"),
                        "jpg" =  list(dev=list(quality=100, width=1600, height=1000, bg="white"), suffix="jpeg"),
                        "postscript" = list(dev=list(paper="special", width=16, height=10, bg="white"), suffix="ps"),
                        "png" =  list(dev=list(width=1600, height=1000, bg="white"), suffix="png"),
                        list(dev=list(width=1600, height=1000,bg="white"), suffix="png"),
                    )
      if(!is.element(dev, c("bmp", "jpeg","png","postscript","jpg")))
        print("Format error, format will be set to PNG")

      ## Set up output name
      fname <- paste("QCHyb", fstart,  plotdef$suffix, sep=".")
      plotdef <- c(plotdef, list(main=paste(fname, ": QC Hyb")))
            
      ###################
      ## Plot
      if(DEBUG) print("start layout")
      if(save)  do.call(dev, maDotsDefaults(opt, c(list(filename=file.path(resdir, fname)), plotdef$dev)) ) 

 
      ## Layout 
      layout(matrix(c(11,1,1,2,11,3,3,7,11,4,4,7,11,5,5,7,11,6,6,7,11,8,9,10), 4, 6),
             height=c(0.5, 2, 2, 4), width = c(12,5.5,2,5.5,2,7)) 

      ## setting controls cols
      ifelse(missing(col), colcode<- setCtlCol(mraw) , colcode <- col)
      
      ##colcode <- unique(as.integer(maControls(mraw))+1)
      ##ctlcode <- levels(maControls(mraw))
      ##names(colcode) <- ctlcode

      ## 1) MA-plot (Before Normalization)
      if(DEBUG) print("start 1")
      qpMAPlots(mraw, addp=TRUE, main="MAPlot: Raw M", ...)
      addLines(mraw)

      ## 2) Boxplot split by Print-tip
      if(DEBUG) print("start 2")
      par(mar=c(5,5,5,3))
      boxplot(mraw, yvar="maM", ylab="M values", main="Boxplot:M by print-tip groups")

      ## 3 & 4) maM (No Normalization)
      if(DEBUG) print("start 3 & 4")
      qpImage(mraw, xvar="maM", main="Spatial plot: Rank(M-raw)")

      ## 5 & 6) maA 
      if(DEBUG) print("start 5 & 6")
      qpImage(mraw, xvar="maA", main="Spatial plot: A")

      ## 7) Compare old and new
      if(DEBUG) print("start 7")
      par(mar=c(5,5,5,3))
      mnorm <- maNorm(mraw)
      rownames(mnorm@maM) <- as.vector(maGeneTable(mnorm)[,"ID"])
      slope <- lm(mnorm@maM[rownames(DEGenes),1] ~ DEGenes[,"Median"])$coef[2]
      if(slope > 0){
        plot(DEGenes[,"Median"],mnorm@maM[rownames(DEGenes),1],xlab="Past average M", ylab="current array", type="n",
             main="Known DE genes: comparison to past average M")
        text(DEGenes[,"Median"],mnorm@maM[rownames(DEGenes),1],DEGenes[,"Ref"])
        }
      else{
        plot(DEGenes[,"Median"],-mnorm@maM[rownames(DEGenes),1],xlab="Past average M", ylab="current array", type="n",
             main="Known DE genes: comparison to past average M")
        text(DEGenes[,"Median"],-mnorm@maM[rownames(DEGenes),1],DEGenes[,"Ref"])
      }
      abline(0, 1, lty=2, lwd=3, col="red")
      abline(-1, 1, lty=2, lwd=3, col="blue")
      abline(1, 1, lty=2, lwd=3, col="blue")

      ## 8 & 9  Red and Green Signal to Noise (background corrected)
      if(DEBUG) print("start 9 & 10")
      qpS2N(bgraw, channel="red", colcode=colcode)
      qpS2N(bgraw, channel="green", colcode=colcode)

      ## 10 Dot Plot
      if(DEBUG) print("start 10")
      if(length(maControls(mraw))!=0)
        {
          qpDotPlots(mraw, x="maA",  col=colcode)
          title( main="Dot plot: Controls A", cex=0.7)
        }
      
      
      ## 11 Title
      if(DEBUG) print("start 11")
      layout(1)
      par(mar=c(2,2,4,2))
      mtext(plotdef$main, line=3)
      mrawheader <- readGPRHeaders(file.path(path,f))
      mtext(paste("Date: ",  mrawheader$DateTime, " :: PMT", mrawheader$PMTGain), line=2, cex=0.8)
      if(DEBUG) print("Done...")
      
      ## Finishing
      if (save == TRUE) {
        cat(paste("save as", fname, "\n"))
        dev.off()
      }
    } ## end for
}## end function


