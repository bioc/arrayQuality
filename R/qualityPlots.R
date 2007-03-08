##########################################################3
## Quality control plots
##
## Author: Jean Yang, Agnes Paquet
## Modified Sept, 14, 2004
## source("~/Projects/madman/Rpacks/arrayQuality/R/qualityPlots.R")



##############################
## rm.na is in package marray


########################
## [Internal Function]
## Set control colors
########################

setCtlCol <- function(mdata)
  {
    colcode <- c(2:(length(unique(as.integer(maControls(mdata))))+1))
    tmp <- levels(maControls(mdata));  ctlcode <- c()
    if("Positive" %in% tmp) ctlcode <- c(ctlcode, "Positive")
    if("probes" %in% tmp) ctlcode <- c(ctlcode, "probes")
    if("Negative" %in% tmp) ctlcode <- c(ctlcode, "Negative")
    if("Empty" %in% tmp) ctlcode <- c(ctlcode, "Empty")
    names(colcode) <- c(ctlcode, setdiff(tmp, ctlcode))
    return(colcode)
  }

########################
## [Internal Function]
## Gene Pix reweighting
########################
gpFlagWt <- function(x)
{
  res <- rep(NA, length(x))
  res[x < -60] <- -100
  res[x == -50] <- 0
  res[x > -50] <- 100
  return(as.matrix(res))
}

##########################################################
## [Internal Function]
## Dot plot:  not really for general use [Internal Function]
##########################################################
qpDotPlots <- function(mdata,  xvar="maA", id="ID", colcode=1, nrep=3, pch=18, ...)
  {
    newdata <- eval(call(xvar, mdata))
    xlim <- range(newdata, na.rm=TRUE)
    Cindex <- maControls(mdata) != "probes"
    #Ctl <- cbind(maInfo(maGnames(mdata)), maControls(mdata))
    ## combined control status and name

    Ctl <- cbind(maInfo(maGnames(mdata)), maControls(mdata), row.names=NULL)
    IDindex <- grep(id, colnames(Ctl))  ## Set ID columns
    y <- split(Ctl, Ctl[,ncol(Ctl)])  ## The last column of Ctl is the control status

    if(length(y[names(y) != "probes"]) != 0)  ## check that there exist control spots
      {
        ## There are control spots
        exty <- lapply(y[names(y) != "probes"], function(x){
          ext <- split(x, x[, IDindex])
          extid <- lapply(ext, function(xx){as.integer(row.names(xx))})
          extid[lapply(extid, length) > nrep]
        })
        exty <- exty[lapply(exty, length) != 0]
        ylim <- c(1, sum(unlist(lapply(exty, length))))

        par(mar=c(4,7,3,2), cex=1)  ## A wide left side to allow for gene names
        plot(1,1, type="n", xlim=xlim, ylim=ylim, axes=FALSE, xlab=xvar, ylab="",...)
        ii <- 1
        for(i in 1:length(exty))
          for(j in 1:length(exty[[i]]))
            {
              ind <- exty[[i]][[j]]
              points(newdata[ind], rep(ii, length(newdata[ind])), pch=pch, col=colcode[names(exty)[i]])
              points(median(newdata[ind], na.rm=TRUE), ii, pch=18, col="black")
              ii <- ii + 1
            }
        axis(1)
        lab <- paste(unlist(lapply(exty, names)), " (n=",unlist(lapply(exty, lapply, length)), ") ", sep="")
        axis(2, at=1:ylim[2], labels=lab, las=2, cex.axis=0.6)
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


##########################################################
## [Internal Function]
## Hexbin Plots
##########################################################

## This is the old function
## qpHexbin <-  function(mdata, main="", ...)
##{
##  par(mar=c(5,4,3,2))
##  y <- maM(mdata)
##  x <- maA(mdata)
##  nalim <- !is.na(x) & !is.na(y)
##  bin <- hexbin(x[nalim],y[nalim])
##  col <- BTY
##  maxcnt  <-  max(bin$cnt)
##  colorcut <-  seq(0, 1, length = min(17,maxcnt))
##
##  yrange <- c(min(maM(mdata), na.rm=TRUE), max(maM(mdata), na.rm=TRUE) + 1)
##  plot(bin$xbnds, bin$ybnds, type = "n", xlab="A", ylab="M", ylim=yrange, main=main)
##  hexagons(bin, colramp = col, colorcut = colorcut, maxcnt = maxcnt)

##  par(mar=c(2,1,2,1))
##  plot(seq(-0, 4, length=11),  -1:9, type="n", axes=FALSE, xlab="", ylab="")
##  hex.legend(legend=4, ysize=8, lcex=1, inner=1,  maxcnt = maxcnt,
##             colorcut= colorcut, colramp = col)
##}

## New version of qpHexbin from Paul Murrell
qpHexbin <- function (mdata, main = "", ...)
{
  y <- maM(mdata)
  x <- maA(mdata)
  nalim <- !is.na(x) & !is.na(y)
  bin <- hexbin(x[nalim], y[nalim])
  col <- BTY
  maxcnt <- max(bin@count)
  colorcut <- seq(0, 1, length = min(17, maxcnt))
  yrange <- c(min(maM(mdata), na.rm = TRUE), max(maM(mdata),
                                na.rm = TRUE) + 1)
  ## move to next traditional graphics plot region
  frame()
  ## Set up grid viewports that correspond to the
  ## current traditional figure region
  require(gridBase)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure)
  pushViewport(viewport(gp=gpar(cex=0.6)))
  ## draw a complete hexbin plot with legend
  plot(bin, xlab = "A", ylab = "M",
       ## ylim = yrange,
       main = main, colramp = col, colorcut = colorcut, maxcnt = maxcnt,
       legend = 1,
       lcex = 1,
       ## ysize = 8,
       ## inner = 1)
       newpage=FALSE)
  upViewport(3)
  ## move to next traditional graphics plot region
  ## i.e., skip over region that used to be used for hexbin legend
  frame()
}


##########################################################
## [Internal Function]
## MA-Plots with bad spots highlighted
##########################################################

qpMAPlots <-  function(mdata, addp=TRUE, main="", ...)
{
  opt <- list(...)
  par(mar=c(5,4,3,2), cex=0.6)
  y <- max(maM(mdata), na.rm=TRUE) + 2
  x <- min(maA(mdata), na.rm=TRUE)
  defs <- maDotsDefaults(opt, maDefaultPar(mdata, x="maA", y="maM", z="maPrintTip"))
  plot(mdata, lines.func=NULL, legend.func=NULL, col="gray",
       ylim=c(min(maM(mdata), na.rm=TRUE),y),  main=main) ## took xlim=c(0,16) out
  abline(h=0,col="blue",lwd=2.5, lty=2)
  legend.func <- do.call("maLegendLines", defs$def.legend)
  legend.func(x, y)
  ## HARD CODE red represents bad spots
  if(addp) points(mdata, subset = maW(mdata) < 0, col="red", pch=18, cex=0.6)
  par(mar=c(5,4,4,2) + 0.1)
}


##########################################################
## [Internal Function]
## Print-tip split loess plots
## Modified from Gordon's code
##########################################################

qpPTLoess <- function(mdata, span=0.4) {
  y <- maM(mdata)
  x <- maA(mdata)
  tmp <- coplot(y~x|factor(maGridCol(mdata)) * factor(maNgr(mdata) - maGridRow(mdata) +1),
                xlab=c("A","Tip Column"),  ylab=c("M","Tip Row"),  pch=".",
                span=span,  show.given=FALSE,  panel=panel.smooth)
  invisible()
}


##########################################################
## [Internal Function]
## Spatial plots
##########################################################

qpImage <- function(mdata, xvar="maM", main="", overlay=NULL, ...)
{
  par(mar=c(2,3,5,2))
  #convert marrayRaw into marrayNorm, needed for assignment
  if(xvar == "maM")
    {
      tmpNorm <- mdata
      if(!is(mdata, "marrayNorm"))
        tmpNorm <- as(mdata, "marrayNorm")
      tmpNorm@maM <- as.matrix(rank(as.numeric(eval(call(xvar, mdata)))))
      tmp <- image(tmpNorm, xvar="maM", main=main, bar=FALSE, overlay=overlay, colorinfo = TRUE, ...)
    }
  else
    tmp <- image(mdata, xvar=xvar, main=main, bar=FALSE, zlim=c(0,16), overlay=overlay, colorinfo=TRUE, ...)
  par(mar=c(2,1, 5,4))
  maColorBar(tmp$x.bar, horizontal = FALSE, col = tmp$x.col,  main = "")
  par(mar=c(5,4,4,2) + 0.1) ## set back to default
}


##########################################################
## [Internal Function]
## Signal to Noise
##########################################################

qpS2N <- function(mdata, channel=c("red", "green"), colcode=1, ...)
  {

    if(channel == "red"){
      if(length(maRb(mdata))!=0)
        S2N <- as.vector(log(maRf(mdata),2) - log(maRb(mdata),2))
      else
        S2N <- as.vector(log(maRf(mdata),2))}

    if(channel == "green"){
      if(length(maGb(mdata))!=0)
        S2N <- as.vector(log(maGf(mdata),2) - log(maGb(mdata),2))
      else
        S2N <- as.vector(log(maGf(mdata),2))}

    lab <- paste("S2N : ","mean:", round(mean(rm.na(S2N)), 2)," : ", "var:", round(var(rm.na(S2N)), 2))
    hist(rm.na(S2N), main=lab, col=channel, nclass=50, freq=FALSE, ylim=c(0,1.1));

    if(length(maControls(mdata))!=0)
      {
        tmp <- split(S2N, maControls(mdata))
        bw <- sd(rm.na(S2N)/4)
        tmp2 <- lapply(tmp, function(x){density(rm.na(x), bw=bw)})
        for(i in 1:length(tmp2))
          lines(tmp2[[i]], lwd=2, col=colcode[names(tmp2)[i]])
        xrange <- range(rm.na(S2N))
        xcood <- xrange[1] + (xrange[2]-xrange[1]) * 0.7
        legend(xcood, 1, names(tmp), lty=1, lwd=2, col=colcode[names(tmp)], cex=0.8)
      }
  }


## mraw == marrayRaw object
## mrawheader == Heading information
maQualityPlots <-  function(mrawObj, headerInfo="",
                            save = TRUE,
                            dev = "png",  #set default to be png
                            col=NULL, badspotfunction=NULL,
                            controlId=c("ID", "Name"),
                            DEBUG=FALSE, ...)
{
  require(hexbin)
  if (DEBUG) print("function starting")
  controlId <- controlId[1]

  #Convert RGList to mraw if needed

  if (setequal(class(mrawObj), "RGList"))
      {
        mrawTmp <- as(mrawObj, "marrayRaw")
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
          tmp <- do.call("gpFlagWt", list(mraw@maW))
          mraw@maW <- tmp 
        }
      else
        if(!is.null(badspotfunction))
          mraw@maW <- do.call(badspotfunction, list(mraw@maW))

      ## setting controls cols
      ifelse(missing(col), colcode<- setCtlCol(mraw) , colcode <- col)
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

      mnorm <- do.call("maNorm", c(list(nbgraw), norm.defs))

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
#<<<<<<< .mine
        #do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=fname), plotdef$dev)), formals(args(dev))))
##      do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=fname), plotdef$dev)), formals(args(dev))))
#=======
        {
          if(plotdef$suffix != "ps")
            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(filename=fname), plotdef$dev)), formals(args(dev))))
          else
            do.call(dev, maDotsMatch(maDotsDefaults(opt, c(list(file=fname), plotdef$dev)), formals(args(dev))))
        }
#>>>>>>> .r11253
      
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
      qpS2N(mraw, channel="red", colcode=colcode)
      qpS2N(mraw, channel="green", colcode=colcode)

      ## 11 maM Dot plot
      if(DEBUG) print("start 11")
      if(length(maControls(mnorm))!=0)
        qpDotPlots(mnorm, xvar="maM", col=colcode,
                   main="Control normalized M", cex.main=0.8, id=controlId)
#      title(main= "Controls normalized M")

      ## 12 maM Dot plot
      if(DEBUG) print("start 12")
      if(length(maControls(mraw))!=0)
        qpDotPlots(mraw, xvar="maA", col=colcode, main="Control A", cex.main=0.8, id=controlId)
       #title(main= "Controls A")

      if(DEBUG) print("start 13")
      ## 13
      layout(1)
      par(mar=c(2,2,4,2))
      mtext(plotdef$main, line=3)
      mtext(headerInfo, line=2, cex = 0.7)
      mtext(paste("Call:", maNormCall(mnorm)[3]), line=1, cex = 0.7)

      ## Finishing
      if(DEBUG) cat("Done...")
      if (save == TRUE) {
        cat(paste("save as", fname, "\n"))
        dev.off()
      }
    }
}

