#################################################
## Separate control plot functions for MEEBO set
## Author: Agnes
## Date modified: 10/31/2005
## from ctlPlots.R. Replaced marray accessors to limma accessors
#################################################

## source("C:/MyDoc/Projects/Statistics/MEEBO/Ash/Oct31/qpMeeboPlotsOct31_limma.R")
## source("C:/MyDoc/Projects/Statistics/MEEBO/01Nov5/qpMeeboPlotsOct31_limma.R")
## load("C:/MyDoc/Projects/Statistics/MEEBO/Ash/MEEBOctrl.RData")

## Mismatch controls plots
## A values vs. %mismatch
## Includes anchored MM, distributed MM and PM
## Work on marrayNorm object containing only 1 array
## Need to set up rownames of A and M to oligo id before
## rownames(mnorm$M) <- rownames(maA(mnorm)) <-
##   make.names(as.vector(maGeneTable(mraw)[,"ID"]),unique=TRUE)

qpMisMatchPlot <- function(mnorm, ctrllist=MEEBOctrl, DEBUG=FALSE)
  {
    if (DEBUG) print("In qpMisMatchPlot")
    ## Set up layout and par
    ##op <- par()
    ##par(mfrow=c(5,2), cex=1.5, cex.axis=0.8, cex.main=1.2,mar=c(2,3,2,3),oma=c(3,3,3,3))
    
    ##Set up y-axis range (samef or all plots)
    yrange <- c(2,16)  ##!!!Needs to be think about more

    ## Get all "normal"  genes for axis set-up
    mMC <- grep("mMC", rownames(mnorm$A)) #24958
    
    #MEEBO set includes 10 unique PM and coresponding MM
    if (DEBUG) print(paste("length ctrllist = ", length(ctrllist), sep=""))
    
    for(i in 1:10)
      {
        if (DEBUG) print(i)
        
        ## Takes control matrix corresponding to PM i
        ctl <- ctrllist[[i]]
        pmid <- as.character(ctl[1, "uniqID"])

        ## Get plot title
        if(as.character(ctl[1, "Symbol"]) != "")
          mtitle <- as.character(ctl[1, "Symbol"])
        else
          mtitle <- as.character(ctl[1, "Probe_Name"])

        idremove <- function(x,y)
          {
            return(is.na(x) | is.na(y) | !is.finite(x) | !is.finite(y))
          }

         ## Plot results for Anchored MM
        
        mma <- grep("a", as.character(ctl[,"MismatchType"]))
  
        plot(as.numeric(ctl[mma,"NumberMM"])*100/70,
             mnorm$A[as.vector(ctl[mma, "uniqID"]),1],
             col="blue",pch=18,
             xlab="% MisMatch", ylab="A value",
             ylim=yrange, xlim=c(-10,105),
             main=mtitle, yaxt="n", axes=FALSE
             )

        ## Plot results for Distributed MM
        mmd <- grep("d", as.character(ctl[,"MismatchType"]))
        points(as.numeric(ctl[mmd,"NumberMM"])*100/70,
               mnorm$A[as.vector(ctl[mmd, "uniqID"]),1],
               col="green",pch=18)

        ## Boxplot for WT
        wt <- grep("WT", as.character(ctl[,"MismatchType"]))
        boxplot(mnorm$A[as.vector(ctl[wt, "uniqID"]),1], width=15, boxwex=3,
                at=-10, add=TRUE, na.rm=TRUE, axes=FALSE, col="red", outline=FALSE,
                medpch=NA, medlty="blank")

        ## Boxplot for negative controls
#        boxplot(mnorm$A[ctrllist[[11]],1], width=10, boxwex=1.2,outline=FALSE,
#                at=100, add=TRUE, na.rm=TRUE, axes=FALSE, col="red",
#                medpch=NA, medlty="blank")
        boxplot(mnorm$A[ctrllist[[11]],1],outline=FALSE,width=15, boxwex=3,
                at=100, add=TRUE, na.rm=TRUE, axes=FALSE, col="red",
                medpch=NA, medlty="blank")
        
 
        goodid <- !idremove(as.numeric(ctl[mma,"NumberMM"])*100/70,
                            mnorm$A[as.vector(ctl[mma, "uniqID"]),1])      
        lines(lowess(as.numeric(ctl[mma,"NumberMM"])[goodid]*100/70,
                     mnorm$A[as.vector(ctl[mma, "uniqID"]),1][goodid]), col="blue")

        goodid <- !idremove(as.numeric(ctl[mmd,"NumberMM"])*100/70,
                            mnorm$A[as.vector(ctl[mmd, "uniqID"]),1])      
        lines(lowess(as.numeric(ctl[mmd,"NumberMM"])[goodid]*100/70,
                     mnorm$A[as.vector(ctl[mmd, "uniqID"]),1][goodid]), col="green")
        
        abline(h=median(mnorm$A[ctrllist[[11]],1], na.rm=TRUE),
               col="black", lty=2)
        abline(h=quantile(mnorm$A[c(mMC),1],c(0.25,0.75,0.9),na.rm=TRUE), col="red", lty=2)
        axis(1, at=seq(0,90,10), label=seq(0,90,10))
        axis(1, at=c(-10,100), label=c("WT", "Neg ctl"))
        axis(2, at=seq(yrange[1], yrange[2],2), las=2)
        axis(4, at=quantile(mnorm$A[c(mMC),1], c(0.25,0.75,0.9,1), na.rm=TRUE),
             label=c(0.25,0.75,0.9,1), las= 2)
        box()
      }
    ##Plot legend in 11th plot
    if(DEBUG) print("Plotting legend")
    plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
    leg.txt <- c("Anchored MM", "Distributed MM",
                 "Median neg. ctl", "% A")
    legend(1,9,leg.txt, col=c("blue", "green", "black", "red"),
           lty=c(1,1,2,2),lwd=3,bty="n")
  }


## Binding Energy boxplot
## Create bins for histogram in BEbin

## Remove controls with low signal:
## set limit to be 75% of negative controls?

## Last element of MEEBOctrl contains Neg Ctrl ids
## One array at a time!!
## Work on marrayNorm object containing only 1 array
## Need to set up rownames of A and M to oligo id before
## Same, but using linear R+G instead of A

BEbin.linear<- function(mraw, ctrllist, DEBUG=FALSE)
  {
    ## take first slide only
    if(dim(mraw)[2]>1)
      mraw <- mraw[,1]

    Anc <- c()
    Dist <- c()
    WT <- c()
    Rc <- mraw$R
    Rc[Rc < 0] <- 0
    Gc <- mraw$G
    Gc[Gc < 0] <- 0

    Int <- mraw$R +mraw$G
    rownames(Int) <- mraw$genes[,"ID"]

    #mMC <- grep("mMC", rownames(mnorm$A)) #24958
    #Alim <- quantile(mnorm$A[mMC], 0.75, na.rm=TRUE)
    Alim <- quantile(Int, 0.75, na.rm=TRUE)
    ##Alim <- quantile(mnorm$A[ctrllist[[11]],], 0.75, na.rm=TRUE)

    if (DEBUG) print(paste("Alim = ", Alim, sep=""))
          
    for(i in 1:(length(ctrllist)-1))
      {
        ctlmat <- ctrllist[[i]]
        mma <- grep("a", as.character(ctlmat[,"MismatchType"]))
        mmd <- grep("d", as.character(ctlmat[,"MismatchType"]))
        wt <- grep("WT", as.character(ctlmat[,"MismatchType"]))
        wtInt <- Int[as.vector(ctlmat[wt,"uniqID"]),]

        ##Normalization : divide by median A for WT
        ##normFactor <- median(mnorm$A[as.vector(ctlmat[wt,"uniqID"]),1], na.rm=TRUE)
        normFactor <- median(as.numeric(wtInt), na.rm=TRUE)
        if (DEBUG) print(paste("Normfactor = ", normFactor,sep=""))

        if( normFactor > Alim)
          {            
            Anc <- rbind(Anc,
                         cbind(Int[as.vector(ctlmat[mma,"uniqID"]),1]/normFactor,
                               round(as.numeric(ctlmat[mma, "BE_PerfectMatch"])/5)*5))
            Dist <- rbind(Dist,
                          cbind(Int[as.vector(ctlmat[mmd,"uniqID"]),1]/normFactor,
                                round(as.numeric(ctlmat[mmd, "BE_PerfectMatch"])/5)*5))
            WT <- rbind(WT,
                        cbind(Int[as.vector(ctlmat[wt,"uniqID"]),1]/normFactor,
                              round(as.numeric(ctlmat[wt, "BE_PerfectMatch"])/5)*5))
          }
      }
    return(list(mmA=Anc, mmD=Dist, WT=WT)) 
  }

## Wrapper function
qpBEplot.linear <- function(mraw, ctrllist=MEEBOctrl,DEBUG=FALSE)
  {
    test <- BEbin.linear(mraw=mraw, ctrllist=ctrllist, DEBUG=DEBUG)
    tmp <- split(test$WT[,1], -test$WT[,2]+2)
    tmp2 <- split(test$mmA[,1], -test$mmA[,2])
    tmp3 <- split(test$mmD[,1], -test$mmD[,2]+1)

    yrange <- range(c(test$WT[,1], test$mmA[,1], test$mmD[,1]), na.rm=TRUE)
    plot(seq(-50,120,1),seq(-50,120,1),
         ylim=yrange,
         type="n", axes=FALSE,
         xlab="Binding energy", ylab="Normalized signal intensity",
         main=paste("Normalized A vs. Binding Energy",
           colnames(mraw$R), sep=" "))
    
    boxplot(tmp, add=TRUE, at=as.numeric(names(tmp)),
            col="red",axes=FALSE, outline=FALSE)
    boxplot(tmp2, add=TRUE, at=as.numeric(names(tmp2)),
            col="blue", axes=FALSE, outline=FALSE)
    boxplot(tmp3, add=TRUE, at=as.numeric(names(tmp3)),
            col="green", axes=FALSE, outline=FALSE)
    abline(h=1, lty=2)
    points(names(tmp),lapply(tmp, median, na.rm=TRUE), pch=15, cex=0.8)
    points(names(tmp2),lapply(tmp2, median, na.rm=TRUE), pch=15, cex=0.8)
    points(names(tmp3),lapply(tmp3, median, na.rm=TRUE), pch=15, cex=0.8)
    lines(lowess(-test$mmA[,2], test$mmA[,1]), col="blue", lwd=2)
    lines(lowess(-test$mmD[,2]+1, test$mmD[,1]), col="darkgreen", lwd=2)
    
    axis(1, at=seq(-50,120,10), label=seq(50,-120,-10))
    axis(2, las=2)
    box()
    legend(-50, yrange[2], c("Wild-Type", "Anchored", "Distributed", "Median"),
           col=c("red", "blue", "green", "black"), pch=15)
  }


## 3' distance plost (tiling controls)
## 11 controls stored in tilingres object (MEEBOTiling.RData)
## Work on marrayRaw object
## Need to assign oligoIds to rownames for maRf and maGf
## rownames(maRf(raw)) <- rownames(maGf(raw)) <-
##   make.names(as.vector(maGeneTable(raw)[,"ID"]),unique=TRUE)

qpTiling <- function(raw, tilingres=tilingres, DEBUG=FALSE)
  {
 
    ##Set up y-axis range (samef or all plots)
    ##yrange <- c(2,16)  ##!!!Needs to be think about more
    
    for(i in 1:11)
      {
        tmp <- tilingres[[i]]
        ##r<-log(ifelse(r>0, r, NA),2)
        Gf <- log.na(raw$G,2)[names(tmp),1]
        Rf <- log.na(raw$R,2)[names(tmp),1]
        #Gf <- maLG(raw)[names(tmp),1]
        #Rf <- maLR(raw)[names(tmp),1]
        distance <- tmp
        mMC <- grep("mMC", rownames(raw$G)) #24958
        mMCA <- (log.na(raw$G,2)[mMC,1] +  log.na(raw$R,2)[mMC,1])/2
        
        plot(-as.numeric(distance), Gf,col="darkgreen",pch=18,
             xlab="3' distance", ylab="Log Intensity",
             main=names(tilingres[i]), ylim=c(2, 16),
             yaxt="n", axes=FALSE)
        points(-as.numeric(distance), Rf,col="red",pch=18)
        
        AnegCt <- (log.na(raw$G,2)[MEEBOctrl[[11]],1] +  log.na(raw$R,2)[MEEBOctrl[[11]],1])/2
        abline(h=median(AnegCt, na.rm=TRUE),
               col="black", lty=2,lwd=2)
        
        axis(1, at=unique(sort(-as.numeric(distance))),
             label=rev(unique(sort(as.numeric(distance)))))
        ##axis(4, at=quantile(mMC, seq(0.5,1, 0.1), na.rm=TRUE),
                                        ##label=seq(0.5,1,0.1), las= 2)
        axis(2, at=seq(2,16,2), las=2) ##to fix, only on left side
        box()
        
        idremove <- function(x,y)
          {
            return(is.na(x) | is.na(y) | !is.finite(x) | !is.finite(y))
          }
        
        goodid <- !idremove(-as.numeric(distance), Gf)
        lines(lowess(-as.numeric(distance)[goodid], Gf[goodid]),
              col="darkgreen", lwd=2)
        goodid <- !idremove(-as.numeric(distance), Rf)
        lines(lowess(-as.numeric(distance)[goodid], Rf[goodid]),
              col="red", lwd=2)
      }
    plot(-as.numeric(distance), Gf,  ylim=c(4, 16),type="n",
         yaxt="n", axes=FALSE, ylab="", xlab="")
    legend(min(-as.numeric(distance)), 12, c("Cy5", "Cy3", "Median Neg ctl"),
           col=c("red", "darkgreen", "black"), lty=c(1,1,2), lwd=3,
           bty="n")
   }


## Example of use for the QC plots

#png("BEplots.png", width=2500, height=600)
#par(mfrow=c(1,5))
#for(i in 1:5)
#  {
#    print(i)
    #plot(1:5, 1:5)
#    qpBEplot(mnorm[,i], MEEBOctrl=MEEBOctrl,DEBUG=FALSE)
#  }
#dev.off()

#png("Tiling_meebo-P1-258.png", width=800, height=2000)
#par(mfcol=c(3,4), cex=1.5, cex.axis=1.2, mar=c(4,4,3,3))
#qpTiling(mraw[,5], tilingres=tilingres)
#dev.off()

##png("Tiling_meebo-14MmnoF5.png", width=1500, height=1200)
##par(mfcol=c(3,4), cex=1.3, cex.axis=1, mar=c(3,4,3,2),oma=c(2,2,4,2))
##qpTiling(RGsub, tilingres=tilingres)
##dev.off()

#png("MMplots.png", width=1200, height=2000)
#par(mfrow=c(5,2), cex=1.5, cex.axis=0.8, mar=c(3,5,3,5))
#qpMMPlots()
#dev.off()
