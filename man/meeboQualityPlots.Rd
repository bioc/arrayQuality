\name{meeboQualityPlots}
\alias{meeboQualityPlots}

\title{Qualitative diagnostic plot for general hybridization, specific
  to MEEBO set.}

\description{
  This function generates diagnostic plots for a qualitative assessment of
  slide quality.  
}
\usage{
meeboQualityPlots(mrawObj, headerInfo = "", save = TRUE, dev = "png",
col = NULL, badspotfunction = NULL, controlId = c("ID", "Name"), seqId =
"SeqID", organism = "Mm", DEBUG = FALSE, ...)
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
    set to "Mm". It is used to retrieve the corresponding Oligo set annotations.}
  
  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{additional arguments}
}

\details{
The right panels of the diagnostic plot contain boxplots of the various
controls in the set, like positive or negative controls, as set in the
marrayRaw of RGlist object. Replicated controls are recognized based on
their MEEBO sequence id.
}

\author{Agnes Paquet}

\seealso{\code{\link{meeboQuality}}, \code{\link{gpQuality}},
  \code{\link{maQualityPlots}},\code{\link{MeeboSpotTypes}},
  \code{\link{controlCodeMeebo}}}

\examples{

if (interactive())
{
require(MEEBO)
datadir <- system.file("Meebo", package="MEEBO")
mraw <- read.GenePix(path=datadir)
maControls(mraw) <- maGenControls(maGnames(mraw),id="ID",controlcode=controlCodeMeebo)
rownames(maRf(mraw)) <- rownames(maRb(mraw)) <- maGeneTable(mraw)[,"ID"]
rownames(maGf(mraw)) <- rownames(maGb(mraw)) <- maGeneTable(mraw)[,"ID"]
meeboQualityPlots(mraw)
}
}
\keyword{hplot}

