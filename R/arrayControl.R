arrayControls <- function (data = NULL, SFGHControlcode = controlCode, id = "ID") 
{
  if(missing(data)) stop("Missing input data")

  gId <- as.vector(data[[id]])
  
  Control <- rep("probes", length(gId))
  for (i in 1:nrow(SFGHControlcode)) {
    position <- grep(as.vector(SFGHControlcode[i, "Pattern"]), 
                     gId)
    if (length(position) > 0) 
      Control[position] <- as.vector(SFGHControlcode[i, "Name"])
    }
  
  return(Control)
}

arrayReplicates <- function(data=NULL, id="ID", cut=5)
  {
    if(missing(data)) stop("Missing input data")
    if(!is.element(id, names(data)))
      stop("Data type is not correct")

    geneList <- as.vector(data[["ID"]])
    y <- table(geneList)
    Replicates <- names(y)[y > cut]
    return(Replicates)
  }



replicatesAvariance <- function(slide = NULL, id="ID")
  {
    if(missing(slide)) stop("Missing input data")
    gprData <- readGPR(slide)

    if(!is.element(id, names(gprData)))
      stop("Can't identify the genes names")
    
    Replicates <- arrayReplicates(gprData)
    gId <- gprData[["ID"]]
    
    index <- NULL
    for(r in Replicates){
      for(i in 1:length(gId)){
        if (r == gId[i]) index <- c(index,i)
      }
    }
    repA <- (log.na(gprData[["RfMedian"]][index],2) +
                log.na(gprData[["GfMedian"]][index],2))/2
    varRepA <- var(repA, na.rm=TRUE)
    return(varRepA)
  }
