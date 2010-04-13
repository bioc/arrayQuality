#########################################################
## Set of functions written to assess quality
## of HEEBO set arrays
## Author: Agnes Paquet
## Date:   06/22/2006
## Using: R2.3
##########################################################

## 1) Get proper controlCode for HEEBO
## 2) Modify qpS2N
## 3) Modify qpDotPlots
## 4) Wrap MMplot, BE plot and tilingPlots into one

## Set up controlCode for Heebo
controlCodeHeebo <-
  structure(c("hHC", "hCP", "hCN", "mCD",
            "Human hHC", "Positive", "Negative", "Doping"),
            .Dim = c(4, 2),
            .Dimnames = list(c("1", "2", "3", "4"),
              c("Pattern", "Name")))

HeeboSpotTypes <-
  structure(c("probes","Human hHC", "Positive", "Negative","Doping",
              "*","hHC*", "hCP*", "hCN*", "hCD*",
            rep("*",5),
            "black","green","red","blue","yellow"),
            .Dim = c(5, 4),
            .Dimnames = list(c("1", "2", "3", "4","5"),
              c("SpotType","ID","Name","Color")))


qpDotPlotsEEBO <- function(mdata,  xvar="maA", id="SeqID", colcode=1, nrep=3, pch=18, meeboAnnot=MEEBOset,...)
  {
    newdata <- eval(call(xvar, mdata))
    xlim <- range(newdata, na.rm=TRUE)
    Cindex <- maControls(mdata) != "probes"

    ## Get ordered sequence id
    matchIds <- match(maGeneTable(mdata)[,"ID"], meeboAnnot[,"uniqID"])
    
    Ctl <- cbind(maInfo(maGnames(mdata)), SeqID=meeboAnnot[matchIds,id],maControls(mdata),row.names=NULL)

    IDindex <- grep(id, colnames(Ctl))  ## Set ID columns
    y <- split(Ctl, Ctl[,ncol(Ctl)])  ## The last column of Ctl is the control status
    
    if(length(y[(names(y) != "probes") & (names(y) != "Doping") & (names(y) != "Mouse mMC") & (names(y) != "Human hHC")])!= 0)
      ## check that there exist control spots
      {
        ## There are control spots
        exty <- lapply(y[(names(y) != "probes") & (names(y) != "Doping")& (names(y) != "Mouse mMC") & (names(y) != "Human hHC")],
                       function(x){
                         ##ext <- split(x, x[, IDindex])
                         ext <- split(x,x[,IDindex],drop=TRUE)
                         extid <- lapply(ext, function(xx){as.integer(row.names(xx))})
                         extid[lapply(extid, length) > nrep]
                       })
        exty <- exty[lapply(exty, length) != 0]
        ## plot replicated controls and boxplot for other types (negative/empty) if any
        numplots <- 0
        plotnames <- c()
        if(("Negative" %in% names(y)) & !("Negative" %in% names(exty)))
          {
            numplots <- numplots+1
            plotnames <- c(plotnames,"Negative")
          }
        if(("Empty" %in% names(y)) & !("Empty" %in% names(exty)))
          {
            numplots <- numplots+1
            plotnames <- c(plotnames, "Empty")
          }
        if(("EMPTY" %in% names(y)) & !("EMPTY" %in% names(exty)))
          {
            numplots <- numplots+1
            plotnames <- c(plotnames, "EMPTY")
          }
        
        ylim <- c(1, (sum(unlist(lapply(exty, length)))+numplots))

        par(mar=c(4,7,3,2), cex=1)  ## A wide left side to allow for gene names
        plot(1,1, type="n", xlim=xlim, ylim=ylim, axes=FALSE, xlab=xvar, ylab="",...)
        ii <- 1
        for(i in 1:length(exty))
          for(j in 1:length(exty[[i]]))
            {
              ind <- exty[[i]][[j]]
              boxplot(newdata[ind], add=TRUE,horizontal=TRUE,
                      at=ii,col=colcode[names(exty)[i]], axes=FALSE)
              ii <- ii + 1
            }
        labs <- c()
        if(numplots>0)
          {
            for(k in 1:numplots)
              {
                boxplot(newdata[as.integer(row.names(y[[plotnames[k]]]))], add=TRUE,horizontal=TRUE,
                        at=ii,col=colcode[plotnames[k]],axes=FALSE)
                ii <- ii+1
                labs <- c(labs,paste(plotnames[k], "(n=",length(row.names(y[[plotnames[k]]])), ") ",sep=""))
              }
          }
        axis(1,at=seq(round(xlim[1]),round(xlim[2]),by=2))

        labstmp <- unlist(lapply(exty,names))
        ## Replace sequence ids by symbols if necessary
        labstmp[grep("SQ", labstmp)] <- 
          as.character(meeboAnnot[match(labstmp[grep("SQ", labstmp)], as.character(meeboAnnot[,"SeqID"])),"Symbol"])        
        
        labsnum <- unlist(lapply(exty, lapply, length))
        labs <- c(paste(labstmp, " (n=", labsnum, ") ",sep=""), labs)
      
        axis(2, at=1:ylim[2], labels=labs, las=2, cex.axis=0.8)
        box()

      } else
    {
         ## There are NO control spots
      plot(1, 1, axes=FALSE, xlab="", ylab="", type="n")
      text(1, 1, "No Control Genes")
      box()
    }
    return()
  }

qpS2Neebo <- function(mdata, channel=c("red", "green"), colcode=1, organism=c("Mm","Hs"),controlId="ID",...)
  {
    if(missing(organism))
      organism <- "Mm"
    else
      organism <- organism[1]
    print(organism)
    
    if(organism == "Mm")
      {
        print("in organism == Mm")
        ##coreCollection <- grep("mMC", rownames(maRf(mdata))) #24958
        coreCollection <- grep("mMC", maGeneTable(mdata)[,controlId])
        if(length(coreCollection) == 0)
          stop("Organism argument doesn't match your array")
      }
    else
      {
        if(organism == "Hs")
          {
            print("in organism == Hs")
            coreCollection <- grep("hHC", rownames(maRf(mdata)))
            print(length(coreCollection))
            if (length(coreCollection) == 0)
              {
                print("coreCollection = 0")
                stop("Organism argument doesn't match your array")
              }
          }
        else
          stop("The input organism does not match MEEBO or HEEBO oligo sets")
      }
    

    ##Same as qpS2N, but using only core collection probes for histogram
    
    if(channel == "red"){
      if(length(maRb(mdata))!=0)
        S2N <- as.vector(log.na(maRf(mdata),2) - log.na(maRb(mdata),2))
      else
        S2N <- as.vector(log.na(maRf(mdata),2))}

    if(channel == "green"){
      if(length(maGb(mdata))!=0)
        S2N <- as.vector(log.na(maGf(mdata),2) - log.na(maGb(mdata),2))
      else
        S2N <- as.vector(log.na(maGf(mdata),2))}

     lab <- paste("S2N : ","mean:", round(mean(rm.na(S2N[coreCollection])), 2)," : ", "var:",
                 round(var(rm.na(S2N[coreCollection])), 2))
    hist(rm.na(S2N[coreCollection]), main=lab, col=channel, nclass=50, freq=FALSE, ylim=c(0,1.1), xlab="Results for Core collection");

    if(length(maControls(mdata))!=0)
      {
        tmp <- split(S2N, maControls(mdata))
        bw <- sd(rm.na(S2N)/4)
        tmp2 <- lapply(tmp, function(x){density(rm.na(x), bw=bw)})
        for(i in 1:length(tmp2))
          {
            lines(tmp2[[i]], lwd=2, col=colcode[names(tmp2)[i]])
          }
        xrange <- range(rm.na(S2N))
        xcood <- xrange[1] + (xrange[2]-xrange[1]) * 0.7
        legend(xcood, 1, names(tmp), lty=1, lwd=2, col=colcode[names(tmp)], cex=1.2)
      }
  }

 
heeboQualityPlots <-  function(mrawObj, headerInfo="",
                               save = TRUE,
                               dev = "png",  #set default to be png
                               col=NULL, badspotfunction=NULL,
                               controlId=c("ID", "Name"),
                               seqId="SeqID",
                               organism="Hs",
                               DEBUG=FALSE, ...)
{
  require(HEEBOdata)
  
  if (DEBUG) print("function starting")
  controlId <- controlId[1]

  if(!("HEEBOset" %in% ls(1)))
    data(HEEBOset)

  #Convert RGList to mraw if needed

  if (setequal(class(mrawObj), "RGList"))
      {
        mrawTmp <- as(mrawObj, "marrayRaw")
        maControls(mrawTmp) <- mrawObj$genes$Status
        rownames(maRf(mrawTmp)) <- rownames(maGf(mrawTmp)) <- rownames(mrawObj$R)
        rownames(maRb(mrawTmp)) <- rownames(maGb(mrawTmp)) <- rownames(mrawObj$R)
        mrawObj <- mrawTmp
        rm(mrawTmp)
      }

  if(DEBUG) print(dim(mrawObj))
  for(i in 1:dim(mrawObj)[2])
    {
      mraw <- mrawObj[,i]
      opt <- list(...)

      ## re-evaluate W
      if (DEBUG) print("Re-evaluate Weight")
      if(missing(badspotfunction)||is.null(badspotfunction))
        {
          tmp <- do.call(gpFlagWt, list(mraw@maW))
          mraw@maW <- tmp 
        }
      else
        if(!is.null(badspotfunction))
          mraw@maW <- do.call(badspotfunction, list(mraw@maW))

      ## setting controls cols
       ifelse(missing(col), colcode<- setCtlCol(mraw) , colcode <- col)

      ## To be fixed: using controlCode to generate maControls, not MeeboSPotTypes
      ## gives character, but colcode is integer
      #if(DEBUG) print("Set up color")
      #if(missing(col))
      #  {
      #    if(length(grep("Color",colnames(MeeboSpotTypes)))>0)
      #       colcode <- MeeboSpotTypes[,"Color"]
      #    else
      #      colcode <- setCtlCol(mraw)
      #  }
      #else
      #  colcode <- col
      
      if(DEBUG) cat("check Control color code", colcode, "\n")

      ## Set up no backgroud data and normalization
      if (DEBUG) print("Set up data and normalization")
      nbgraw <- mraw;
      if(length(mraw@maGb) != 0)
        nbgraw@maGb <- nbgraw@maRb <- matrix(0,0,0)

      if (DEBUG) print("Reading normalization parameters")
      if (DEBUG) print(opt)
      defs <- list(norm="p")
      norm.defs <- maDotsMatch(maDotsDefaults(opt, defs), formals(args("maNorm")))

      if (DEBUG) cat("Using normalization method:  ", norm.defs$norm, "\n")

      mnorm <- do.call(maNorm, c(list(nbgraw), norm.defs))

      ## Set up output name
      if (DEBUG) print("Name the output file")
      tmp <- unlist(strsplit(colnames(mraw@maGf), "/"))
      tmp3 <- unlist(strsplit(tmp[length(tmp)], "\\."))
      tmp2 <-  sub("/", "", tmp3)

      if(length(tmp2) > 1)
        fstart <- paste(tmp2[-length(tmp2)], collapse=".")
      else
        fstart <- tmp2

      if (DEBUG) print(fstart)


  ###################
  ## Setting up output device
  ###################
      if (DEBUG) print("Name of output device")
      plotdef <- switch(dev,
                        "bmp" = list(dev=list(width=1600, height=1200, bg="white"), suffix="bmp"),
                        "jpeg" = list(dev=list(quality=100, width=1600, height=1200, bg="white"), suffix="jpeg"),
                        ##"postscript" = list(dev=list(paper="special", width=32.0, height=24.0, bg="white"), suffix="ps"),
                        "postscript" = list(dev=list( bg="white"), suffix="ps"),
                        "png" =  list(dev=list(width=1600, height=1200, bg="white"), suffix="png"),
                        list(dev=list(width=1600, height=1200,bg="white"), suffix="png"),
                        )
      if(!is.element(dev, c("bmp", "jpeg","png","postscript")))
              {
        print("Format error, format will be set to PNG")
        dev = "png"
      }

      fname <- paste("diagPlot", fstart,  plotdef$suffix, sep=".")
      if (DEBUG) print(fname)
      plotdef <- c(plotdef, list(main=paste(fname, ": Qualitative Diagnostic Plots")))

      ###################
      ## Plot
      ###################
      ## Match args and calls function
      ## Start device and layout
      if(DEBUG) print("start layout")

      if(save)
        {
          if(plotdef$suffix != "ps")
            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=fname), plotdef$dev)), formals(args(dev))))
          else
            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=fname), plotdef$dev)), formals(args(dev))))
        }

      
      ## Feb 23, 2005: JYH
      ## Modified width from 1.5 to 2 for the color bar.  This stops the error
      ## Figure margins too large when dev="postscript".  It was not a problem
      ## for any other dev setting. 
      
      ##      layout(matrix(c(14, 1,2,2, 14,0,3,3, 14,4,6,6, 14, 5, 7, 7, 14, 8, 10, 11,
      ##                      14, 9, 10, 11, 14, 12, 13, 13), 4, 7),
      ##             height=c(1, 10, 5, 5), width = c(11, 2, 5, 1.5 ,5, 1.5, 7))
      layout(matrix(c(14, 1,2,2, 14,0,3,3, 14,4,6,6, 14, 5, 7, 7, 14, 8, 10, 11,
                      14, 9, 10, 11, 14, 12, 13, 13), 4, 7),
             height=c(1, 10, 5, 5), width = c(11, 2, 5, 2 ,5, 2, 7))

      
      ## 1) Split MA-plot (Before Normalization)
      if(DEBUG) print("start 1")
      qpMAPlots(nbgraw, addp=TRUE, main="MA-Plot :: raw", ...)
      if(!is.null(opt$norm))
        {
          if(opt$norm == "p")
            addLines(nbgraw)
        }
      else
        addLines(nbgraw)

      ## 2) HEXbin MA-plot (After Normalization)
      if(DEBUG) print("start 2")
      qpHexbin(mnorm, main="MA-Plot :: Norm")

      ## 3 & 4) maM (Before Normalizaton)
      if(DEBUG) print("start 3, 4")
      qpImage(nbgraw, xvar="maM", main="Spatial: Rank(M-Raw)")

      ## 5 & 6) maM (After Normalization)
      if(DEBUG) print("start 5 & 6")
      ov.sub <- as.vector(maW(mnorm)) < 0
      qpImage(mnorm, xvar="maM", main="Spatial: Rank(M-Norm)", overlay=ov.sub)

      ## 7 & 8) maA
      if(DEBUG) print("start 7 & 8")
      qpImage(nbgraw, xvar="maA", main="Spatial: A")

      ## 9 & 10  Red and Green Signal to Noise (background corrected)
      if(DEBUG) print("start 9 & 10")
      if(DEBUG)
        print(paste("organism = ",organism,sep=""))
      qpS2Neebo(mraw, channel="red", colcode=colcode, organism=organism,controlId=controlId)
      qpS2Neebo(mraw, channel="green", colcode=colcode, organism=organism, controlId=controlId)

      ## 11 maM Dot plot -- Replaced by controls boxplot
      if(DEBUG) print("start 11")
      if(length(maControls(mnorm))!=0)
        qpBoxplotMeebo(mraw,xvar="maA", col=colcode, main="Control A", cex.main=0.8, id=controlId,meeboAnnot=HEEBOset)

      ## 12 maM Dot plot
      if(DEBUG) print("start 12")
      if(length(maControls(mraw))!=0)
        qpDotPlotsEEBO(mraw, xvar="maM", col=colcode, main="Control M", cex.main=0.8, id=seqId, meeboAnnot=HEEBOset)
       #title(main= "Controls A")

      if(DEBUG) print("start 13")
      ## 13
      layout(1)
      par(mar=c(2,2,4,2))
      mtext(plotdef$main, line=3)
      mtext(headerInfo, line=2, cex = 0.7)
      ##mtext(paste("Call:", maNormCall(mnorm)[3]), line=1, cex = 0.7)
      mtext(paste("Normalization: ",norm.defs$norm), line=1, cex = 0.7)
      
      ## Finishing
      if(DEBUG) cat("Done...")
      if (save == TRUE) {
        cat(paste("save as", fname, "\n"))
        dev.off()
      }
    }
}

