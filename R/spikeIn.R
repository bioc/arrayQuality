#################################################
## Separate control plot functions for MEEBO set
## Author: Agnes and Yuanyuan
## Date modified: 11/02/2005
## Spike-ins plots -- Using limma functions
#################################################

##source("C:/MyDoc/Projects/madman/Rpacks/arrayQuality/R/SpikeIn.R")

###############################################################
## LOAD REQUIRED PACKAGES and OBJECTS

library(marray)
library(RColorBrewer)
##load("MEEBOset.RData")


###############################################################
## Spike parsing functions
###############################################################

## input: text file describing the spike-ins (names/cy3-cy5/seqID/Type)
## output: list of matrix containing info for each spike type

readSpikeTypes <- function(file="DopingTypeFile2.txt",path=NULL,cy5col="MassCy5",cy3col="MassCy3",...)
  {
    spikeTypes <- readSpotTypes(file=file,path=path,...)
    spikeTypes <- spikeTypes[!is.na(spikeTypes[,cy5col]),]
    if(length(grep("Type", as.character(colnames(spikeTypes))))>0)
      spikelist <- split(spikeTypes,spikeTypes[,"Type"])
    else
      stop("Wrong SpikeTypeFile")
    return(spikelist)
  }

## Examples:
#spikes <- readSpikeTypes() ##List of matrix by  different types of spikes
#spikes <- readSpikeTypes(file="StanfordDCV1.7complete.txt",
#                         cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7")

## -------------------------------------------------------------------------
## getSpikeIds
## Input: list of matrix, one for for each type of spike
## Output: list of list:
## For each spike-type: look at all individual spikes,
## if on array, find all replicates (using seqID) and return their MEEBO ids
## 
## Some spikes are not used- Remove from search

getSpikeIds <- function(spikeList, id="SeqID", annot=MEEBOset, namecol=c("Symbol","Name"))
  {
    SpikeIds <- c()
    if(length(namecol)>2)
      namecol=namecol[1]
    ## find all replicated sequences/probe_Name
    for(type in 1:length(spikeList))
      {
        usedSeq <- spikeList[[type]][!is.na(spikeList[[type]][,id]),id]
        TmpIds <- lapply(spikeList[[type]][,id],
                         function(x){as.character(annot[grep(as.character(x),as.character(annot[,id])),"uniqID"])})
        if(nchar(spikeList[[type]][1,namecol])>0)
          names(TmpIds) <- spikeList[[type]][,namecol]
        else
          names(TmpIds) <- spikeList[[type]][,"Probe_Name"]
        SpikeIds <- c(SpikeIds,list(TmpIds))
      }
    names(SpikeIds) <- names(spikeList)
    return(SpikeIds)
  }

## Example:
##test <- getSpikeIds(spikeList=spikes)

## ---------------------------------------------------------------------------
## getSpikeIndex
## Input: Input: list of matrix, one for for each type of spike
## Output: list of list:
## For each spike-type: look at all individual spikes,
## if on array, find all corresponding Meebo ids
## return index of each replicate in all-array index

getSpikeIndex <- function(spikeList,id="SeqID",annot=MEEBOset,gnames, namecol=c("Symbol","Name"))
  {
    if(length(namecol)>1)
      namecol=namecol[1]
    SpikeIds <- getSpikeIds(spikeList=spikeList,id=id,annot=annot, namecol=namecol)
    SpikeIndex <- lapply(SpikeIds,function(x){lapply(x,function(y){match(y,gnames)})})
  }

## Example:
##test2<- getSpikeIndex(spikes,gnames=maGeneTable(mnorm)[,"ID"])

## ---------------------------------------------------------------------------
## Get all spikes at once, no distinctions
## Input: spike type file
## Output: list of all spike-ins Meebo ids for each spike
readAllSpikes <- function(file="DopingTypeFileTest.txt",path=NULL,cy5col="MassCy5",
                          cy3col="MassCy3",id="SeqID",annot=MEEBOset,...)
  {
    ## read everything into a matrix
    spikesTypes <- readSpotTypes(file=file,path=path,...)
    ## remove missing data NA or ""
    spikesTypes <- spikesTypes[!is.na(spikesTypes[,cy5col]),]
    ## get all replicates on array using SeqId
    TmpIds <- sapply(spikesTypes[,id],
                     function(x){as.character(annot[grep(as.character(x),as.character(annot[,id])),"uniqID"])})
    ##
    return(TmpIds)
  }

## spikeids <- readAllSpikes()


###############################################################
## Individual Spike-type plotting functions
###############################################################

## Spike.MMplot: 
## Input: results of readSpikeTypes
## output: expected ratio vs. Observed ratios, log scale

Spike.MMplot <- function(mval,spikeList, id="SeqID", gnames=RG$genes[,"ID"],
                         annot=MEEBOset,cy5col="MassCy5",cy3col="MassCy3",namecol="Symbol")
  {
    SpikeIndex <- getSpikeIndex(spikeList,id=id,annot=annot,gnames=gnames,namecol=namecol)
##    print(SpikeIndex[[1]])
##    print(length(spikeList))
    for(type in 1:length(spikeList))
      {
 ##       print("here")
        expected <- log.na(as.numeric(spikeList[[type]][,cy5col]) / as.numeric(spikeList[[type]][,cy3col]),2)
        max <- max(expected,na.rm=TRUE)
        min <- min(expected,na.rm=TRUE)
        ##set the color
        ##separate MJ and Ambion/Strata
        
        if (length(SpikeIndex[[type]])>16)
          {
           print("More than 16 spikes") 
           colvect <- "blue"
           ##Plot dots
           plot((min-1):(max+1),(min-1):(max+1),type="n",xlab="expected ratios", ylab="observed ratios",
                main=spikeList[[type]][1,"Type"])
           for(i in 1:length(SpikeIndex[[type]]))
             {
               yy <- mval[SpikeIndex[[type]][[i]]]
               points(rep(expected[i],length(yy)),yy,col=colvect, pch=19,cex=1)
             }
           abline(0,1,lty=2)
           ##plot legend
           plot(rep(1,length(SpikeIndex[[type]])),1:length(SpikeIndex[[type]]),
                type="n", xlab="", ylab="", axes=FALSE, xlim=c(1,1))
         }
        else
          {
            colvect <- rep(brewer.pal(ceiling(length(SpikeIndex[[type]])/2),"Set1"),2)
            ##plot
            plot((min-1):(max+1),(min-1):(max+1),type="n",xlab="expected ratios", ylab="observed ratios",
                 main=spikeList[[type]][1,"Type"])
            for(i in 1:length(SpikeIndex[[type]]))
              {
                yy <- mval[SpikeIndex[[type]][[i]]]
                points(rep(expected[i],length(yy)),yy,col="black", pch=LETTERS[i],cex=1.5)
                if(length(grep("EC",as.character(spikeList[[type]][i,"Probe_Name"])))>0)
                  points(expected[i], median(yy,na.rm=TRUE), col=colvect[i], pch=19,cex=3)
                else
                  {
                    points(expected[i], median(yy,na.rm=TRUE), col=colvect[i], pch=17,cex=3)
                    ##print(median(yy,na.rm=TRUE))
                  }
              }
            
            abline(0,1,lty=2)
            op <- par(mar=c(5,2,5,10))
            plot(rep(1,length(SpikeIndex[[type]])),1:length(SpikeIndex[[type]]),
                 type="n", xlab="", ylab="", axes=FALSE, xlim=c(1,1))
            for(i in 1:length(SpikeIndex[[type]]))
              {
                points(1, i,col=colvect[i],
                       pch=ifelse(length(grep("EC",spikeList[[type]][i,"Probe_Name"]))>0,19,17),cex=3)
                axis(4,at=1:length(SpikeIndex[[type]]),
                     paste(LETTERS[1:length(SpikeIndex[[type]])], as.character(names(SpikeIndex[[type]])),sep=":"),
                     las=2,cex.axis=1.5,tick=FALSE,line=-3)
              }
            par(op)
          }
      }
  }

## Example
#png("testMMplot_MJs.png",width=1000,height=1000)
#par(mfrow=c(2,2),cex.main=2,cex.lab=1.5,oma=c(4,4,4,4))
#Spike.MMplot(Mval, spikes[1],cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7")
#Spike.MMplot(Mval, spikes[2],cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7")
#dev.off()

## --------------------------------------------------------------------------------

## Spike.Sensitivity
## concentration (mass) vs raw signal for each channel
## Input: output from read.SpikeTypes
## Output: boxplot of log raw signal / log concentration
## Bg sub or not depending on input

Spike.Sensitivity <- function(rawobj,spikeList, id="SeqID", gnames=RG$genes[,"ID"], annot=MEEBOset,cy5col="MassCy5",cy3col="MassCy3",namecol=c("Name","Symbol"))
  {
    if(length(namecol)>2)
      namecol <- namecol[1]
    
    SpikeIndex <- getSpikeIndex(spikeList,id=id,annot=annot,gnames=gnames,namecol=namecol)
    for(type in 1:length(spikeList))
      {
        expectedR <- log.na(as.numeric(spikeList[[type]][,cy5col]),2)
        expectedG <- log.na(as.numeric(spikeList[[type]][,cy3col]),2)
        max <- max(c(expectedR,expectedG),na.rm=TRUE)
        min <- min(c(expectedR,expectedG),na.rm=TRUE)
        Cy5 <- c()
        Cy3 <- c()
        massCy5 <- c()
        massCy3 <- c()
        for(i in 1:length(SpikeIndex[[type]]))
          {
            ##yy <- maLR(rawobj)[SpikeIndex[[type]][[i]]]
            yy <- log.na(rawobj$R[SpikeIndex[[type]][[i]]],2)
            Cy5 <- c(Cy5,yy)
            massCy5 <- c(massCy5,rep(expectedR[i],length(yy)))
            yy <- log.na(rawobj$G[SpikeIndex[[type]][[i]]],2)
            Cy3 <- c(Cy3,yy)
            massCy3 <- c(massCy3,rep(expectedG[i],length(yy)))
          }
        resCy5<- split(data.frame(cbind(massCy5,Cy5)),factor(massCy5))
        resCy3<- split(data.frame(cbind(massCy3,Cy3)),factor(massCy3))
        plot((min-1):(max+1),(min-1):(max+1),type="n",xlab="Log2 mass",
             ylab="Signal intensity",ylim=c(0,16),xlim=c(0,10),
             main=paste(spikeList[[type]][1,"Type"],"Cy5",sep=" "))
        for(i in 1:length(resCy5))
          {
            boxplot(resCy5[[i]][,2], at=resCy5[[i]][1,1],col="red",add=TRUE)
          }
        plot((min-1):(max+1),(min-1):(max+1),type="n",xlab="Log2 mass",
             ylab="Signal intensity",ylim=c(0,16),xlim=c(0,10),
             main=paste(spikeList[[type]][1,"Type"],"Cy3",sep=" "))
        for(i in 1:length(resCy3))
          {
            boxplot(resCy3[[i]][,2], at=resCy3[[i]][1,1],col="green",add=TRUE)
          } 
      }
  }

## Example
#png("Sensitivity.png", width=800,height=1000)
#par(mfrow=c(2,2))
#Spike.Sensitivity(nbg,spikes)
#dev.off()

## --------------------------------------------------------------------
## Spike.MM.Scatter
# MMplot: scatter plot + expected ratio superimposed
## Input: results from readSpikeTypes
## Output: scatter plot of observed ratio and expected ratio for each spike

Spike.MM.Scatter <- function(mval,spikeList, id="SeqID", gnames=RG$genes[,"ID"], annot=MEEBOset,
                             cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7",namecol=NULL)
  {
    SpikeIndex <- getSpikeIndex(spikeList,id=id,annot=annot,gnames=gnames,namecol=namecol)
    yrange <- range(mval[unlist(SpikeIndex)],na.rm=TRUE)
    for(type in 1:length(SpikeIndex))
      {
        expected <- log.na(as.numeric(spikeList[[type]][,cy5col]) / as.numeric(spikeList[[type]][,cy3col]),2)
        ##set the color
        #if (length(SpikeIndex[[type]])>16)
        #  colvect <- "blue"
        #else
        #  colvect <- rep(brewer.pal(ceiling(length(SpikeIndex[[type]])/2),"Set1"),2)
        ##plot
        yrange <- range(yrange,expected,na.rm=TRUE)
        plot(1:length(SpikeIndex[[type]]),1:length(SpikeIndex[[type]]),
             type="n",xlab="", ylab="Observed ratios",ylim=yrange,axes=FALSE,
             main=spikeList[[type]][1,"Type"])
        for(i in 1:length(SpikeIndex[[type]]))
          {
            yy <- mval[SpikeIndex[[type]][[i]]]
            points(rep(i,length(yy)),yy,col="black", pch=18,cex=1.5)
            points(i,expected[i],pch=19,cex=2,col="red")
          }
        abline(h=c(-2,-1,0,1,2),lty=2)
        axis(1,at=1:length(SpikeIndex[[type]]),labels=names(SpikeIndex[[type]]),las=2)
        axis(2,las=2);box()
      }
  }

## Example
#png("MM.scatter.png", width=600,height=800)
#par(mar=c(7,4,4,4),mfrow=c(2,1),oma=c(2,2,2,2))
#Spike.MM.Scatter(Mval,spikes)
#dev.off()

##-------------------------------------------------------------------------
## Cy3 vs Cy5 raw signal, linear scale
## Bg sub depends on input
## All spikes
## Input: spikeTypes: results of readSpotTypes
Spike.Cy5vsCy3 <- function(rawobj,spikesTypes, id="SeqID",
                           annot=MEEBOset,gnames=RG$genes[,"ID"],
                           cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7",...)
  {
    ## read everything into a matrix
##    spikesTypes <- readSpotTypes(file=file,path=path)
    ## remove missing data NA or ""
    spikesTypes <- spikesTypes[!is.na(spikesTypes[,cy5col]),]
    ## get all replicates on array using SeqId
    TmpIds <- sapply(spikesTypes[,id],
                     function(x){as.character(annot[grep(as.character(x),as.character(annot[,id])),"uniqID"])})
    ## Get the colors
   ## expectedRatio <- spikesTypes[,cy5col]/spikesTypes[,cy3col]
    expectedRatio <- as.numeric(spikesTypes[,cy5col])/as.numeric(spikesTypes[,cy3col])
    
    ratioSpiked <- sort(unique(expectedRatio))
    colvect <- rainbow(length(ratioSpiked))
    indexcol <- match(expectedRatio,ratioSpiked)

    ## Get the plot
    axisrange <- range(log.na(c(rawobj$R,rawobj$G),2),na.rm=TRUE)
    axisrange[2] <- axisrange[2]+1.5
    
    plot(1:length(TmpIds),1:length(TmpIds),type="n",xlim=axisrange,
         ylim=axisrange,xlab="Log2 Cy5 signal", ylab="Log2 Cy3 signal",
         main="Comparison of Cy5 and Cy3 raw signal intensity")      
    for(i in 1:length(TmpIds))
      {
        index <- match(TmpIds[[i]],as.character(gnames))
        points(log.na(rawobj$R[index,],2),log.na(rawobj$G[index,],2),
               ##pch=18,col=1)
        pch=18,col=colvect[indexcol[i]])

      }
    legend(axisrange[2]-2.2,axisrange[2]-0.5,"Expected log2 ratios",
           bty="n",cex=0.7)
    legend(axisrange[2]-1,axisrange[2]-1,pch=18,col=colvect,bty="o",
           as.character(log2(ratioSpiked)))
    abline(0,1,lty=2)
  }

##----------------------------------------------------------------     
## Spike.Sensitivity indivudual boxplots
## concentration (mass) vs raw signal for each channel
## Input: output from read.SpikeTypes
## Output: boxplot of log raw signal ordered by input mass
## Bg sub or not depending on input

##Example
#png("testIndiBplot.png",width=1000,height=1000)
#par(mfrow=c(2,2), oma=c(3,3,3,3),mar=c(5,4,5,4))
#Spike.Individual.Sensitivity(RGsub,spikeList,namecol="Symbol")
#dev.off()

Spike.Individual.Sensitivity <- function(rawobj,spikeList=NULL, id="SeqID", gnames=RG$genes[,"ID"], annot=MEEBOset,cy5col="MassCy5",cy3col="MassCy3",namecol=c("Name","Symbol"),ctllist=MEEBOctrl,DEBUG=FALSE)
  {
    if(DEBUG) print("Starting Spike.Individual.Sensitivity")
    if(length(namecol)>2)
      namecol <- namecol[1]
 #   Aval <- (log.na(rawobj$R,2) + log.na(rawobj$G,2))/2
    
    spikeIndex <- getSpikeIndex(spikeList,id=id,annot=annot,gnames=gnames,namecol=namecol)
    negCtlCy5 <- log.na(rawobj$R[ctllist[[11]],],2)
    negCtlCy3 <- log.na(rawobj$G[ctllist[[11]],],2)
    for(type in 1:length(spikeIndex))
      {
        orderCy5 <- sort(as.numeric(spikeList[[type]][,cy5col]),index.return=TRUE)$ix
        orderCy3 <- sort(as.numeric(spikeList[[type]][,cy3col]),index.return=TRUE)$ix

        plot(1:(length(spikeIndex[[type]])+1),1:(length(spikeIndex[[type]])+1),
             type="n",ylim=c(0,16),xlab="",axes=FALSE,
             ylab="Log2 signal",main=paste(spikeList[[type]][1,"Type"]," : log2(Cy5) by increasing concentration",sep=" "))
        boxplot(negCtlCy5,add=TRUE,col="blue",at=1,axes=FALSE,na.rm=TRUE)
        for(j in 1:length(orderCy5))
          {
            boxplot(log.na(rawobj$R[spikeIndex[[type]][[orderCy5[j]]],],2),add=TRUE,col="red",
                    at=(j+1),axes=FALSE,na.rm=TRUE)
          }
        axis(1,at=1:(length(orderCy5)+1),labels=c("Neg Ctrl",spikeList[[type]][orderCy5,namecol]),las=2)
        axis(2,las=2, at=seq(0,16,2));box()
#        axis(4, at=quantile(Aval, c(0.25,0.5,0.75,0.9,1), na.rm=TRUE),
#             label=c(0.25,0.5,0.75,0.9,1), las= 2)
#        abline(h=quantile(Aval, c(0.25,0.5,0.75,0.9,1), na.rm=TRUE),lwd=2,lty=2)

        plot(1:(length(spikeIndex[[type]])+1),1:(length(spikeIndex[[type]])+1),
             type="n",ylim=c(0,16),xlab="",axes=FALSE,
             ylab="Log2 signal",main=paste(spikeList[[type]][1,"Type"]," : log2(Cy3) by increasing concentration",sep=" "))
        boxplot(negCtlCy3,add=TRUE,col="blue",at=1,axes=FALSE)
        for(j in 1:length(orderCy3))
          {
            boxplot(log.na(rawobj$G[spikeIndex[[type]][[orderCy3[j]]],],2),add=TRUE,col="green",
                    at=(j+1),axes=FALSE)
          }
        axis(1,at=1:(length(orderCy3)+1),labels=c("Neg Ctrl",spikeList[[type]][orderCy3,namecol]),las=2)
        axis(2,las=2,at=seq(0,16,2));box()
      }
  }
        
 
##Spike.Individual.Sensitivity(RGsub,spikeList[1],namecol="Symbol")
