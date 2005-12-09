\name{readSpikeTypes}
\alias{readSpikeTypes}

\title{Read Spike Types File}

\description{
  Read a table containing information about the doping control mixture used in
  the hybridization.
}

\usage{
readSpikeTypes(file = "DopingTypeFile2.txt", path = NULL, cy5col = "MassCy5", cy3col = "MassCy3", ...)
}

\arguments{
  \item{file}{A "character" string giving the name of the file specifying
    the doping controls used.}

  \item{path}{A "character" string giving the directory containing the
    file. By default this is set to the current working directory
    (".").}
  
  \item{cy5col}{A "character" string describing the name of the column
    corresponding to the controls labelled with Cy5.}
  
  \item{cy3col}{A "character" string describing the name of the column
    corresponding to the controls labelled with Cy3.}
  
  \item{\dots}{Additional arguments passed to "readSpotTypes"}
}

\details{
   The file is a text file with rows corresponding to doping controls
   and columns describing various experimental conditions. It must
   contain an oligo sequence identifier for each control, the spike type
   (e.g. Ambion or MJ) and the mass of each oligo spiked in each
   channel. By default, this function assumes that the mass unit are the same.
}

\value{
  A list of n matrices, each matrix containing information about a unique
  type of spiked controls.  
}

\author{Agnes Paquet}

\seealso{\code{\link{readSpotTypes}}}

\examples{
datadir <- system.file("Meebo", package="arrayQuality")
if (interactive())
{
spikes <-
readSpikeTypes(file="StanfordDCV1.7complete.txt",path=datadir,cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7")
}
}
\keyword{data}% at least one, from doc/KEYWORDS
