\name{heeboQualityPlots}
\alias{heeboQualityPlots}

\title{Qualitative diagnostic plot for general hybridization, specific
  to HEEBO set.}

\description{
  This function generates diagnostic plots for a qualitative assessment of
  slide quality.  
}
\usage{
heeboQualityPlots(mrawObj, headerInfo = "", save = TRUE, dev = "png",
col = NULL, badspotfunction = NULL, controlId = c("ID", "Name"), seqId =
"SeqID", organism = "Hs", DEBUG = FALSE, ...)
}

\arguments{
  \item{mrawObj}{\code{marrayRaw}  or \code{RGList} object representing
    the slides to be tested.}

  \item{headerInfo}{Text to be used as header in the diagnostic plot.}

  \item{save}{Logical. If 'TRUE', the plot is saved to a file.}

  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}

  \item{col}{Vector of colors to use to describe different
    categories of spots.}

  \item{badspotfunction}{Function to use for bad spots on the array.}

  \item{controlId}{Character string. Name of the column of the gpr
    file used to define controls (or provides the name of the probes).}

  \item{seqId}{A "character" string naming the column of the MEEBO
    annotation to use to retrieve replicated sequences.}

  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Hs". It is used to retrieve the corresponding Oligo set annotations.}
  
  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{additional arguments}
}

\details{
The right panels of the diagnostic plot contain boxplots of the various
controls in the set, like positive or negative controls, as set in the
marrayRaw of RGlist object. Replicated controls are recognized based on
their HEEBO sequence id.
}

\author{Agnes Paquet}

\seealso{\code{\link{heeboQuality}}, \code{\link{gpQuality}},
  \code{\link{maQualityPlots}},\code{\link{HeeboSpotTypes}},
  \code{\link{controlCodeHeebo}}}

\examples{
if (interactive())
{
require(HEEBOdata)
datadir <- system.file("Heebo", package="HEEBOdata")
gal <- readGAL("hoc.gal", path=datadir)
RG <-
read.maimages(files=c("63421.gpr"), path=datadir, source="genepix.median",other.columns="Flags")
RG$genes <- gal
RG$printer <- getLayout(RG$genes)
RG$genes$Status <- controlStatus(HeeboSpotTypes,RG,verbose=TRUE)
rownames(RG$R) <- rownames(RG$G) <- RG$genes[,"ID"]
heeboQualityPlots(RG)
}
}
\keyword{hplot}

