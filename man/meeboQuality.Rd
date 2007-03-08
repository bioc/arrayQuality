\name{meeboQuality}
\alias{meeboQuality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostics plots designed for MEEBO set controls}
\description{
  This component generates several exploratory plots customized to the
  various types of controls provided in the MEEBO set. All results are
  saved as an image. Tested on GenePix format only (11-18-2005). For
  more details about the plots, please refer to the HTML description.
}
\usage{
meeboQuality(fnames = NULL, path = ".", galfile = NULL, source
="genepix.median", other.columns = c("Flags"),controlMatrix=MeeboSpotTypes,controlId = c("ID", "Name"),
DOPING = TRUE,meeboSetQC = TRUE, SpotTypeFile = NULL, SpikeTypeFile =
NULL, cy3col = "CY3.ng._MjDC_V1.7", cy5col = "CY5.ng._MjDC_V1.7", id =
"SeqID", namecol = c("Symbol", "Name"), annot = NULL, bgMethod = "none", normMethod =
"p", diagnosticPlot = TRUE, output = TRUE, resdir = ".", dev = "png", organism = "Mm", DEBUG = FALSE,
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnames}{A "character" string naming the input files.}

  \item{path}{A "character" string representing the data directory.  By
    default this is set to the current working directory (".").}
  
  \item{galfile}{A "character" string naming the file descrining the
    layout of the array. If missing, meeboQuality will read the layour
    from the gpr file.}
  
  \item{source}{A "character" string specifing the image analysis program
    which produced the output files. See ?read.maimages in package limma
    for more details.}
  
  \item{other.columns}{See ?read.maimages in package limma for more
    details.}

  \item{controlMatrix}{A character matrix of n by 2 or more
    columns. One column should contain a few regular expression of
    spotted probe sequences and another column should contain the
    corresponding control status. By default, controlMatrix is set to
    MeeboSpotTypes.}
  
  \item{controlId}{Character string. Name of the column of the gpr
    file used to define controls.}
  
  \item{DOPING}{Logical. If 'TRUE', doping controls quality plots are
    generated.}
  
  \item{meeboSetQC}{Logical. If 'TRUE', mismatch and tiling controls
    quality plots are generated.}
  
  \item{SpotTypeFile}{A "character" string representing the name of the
    file containing spot type description for the array.}
  
  \item{SpikeTypeFile}{A "character" string representing the name of the
    file containing doping control information. See HTML description for
    more details.}
  
  \item{cy3col}{A "character" string representing the name of the column
    of the SpiketypeFile containing the quantity of each control spiked
    in the Cy3 channel.}
  
  \item{cy5col}{A "character" string representing the name of the column
    of the SpiketypeFile containing the quantity of each control spiked
    in the Cy5 channel.}
  
  \item{id}{A "character" string describing which column of the MEEBO
    annotation should be used to retrieve replicated oligos,
    e.g. "SeqID".}
  
  \item{namecol}{A "character" string describing which column of the
    SpiketypeFile should be used in the legend.}
  
  \item{annot}{A "character" string describing which R object should be
    used to look-up probes annotations. By default, it is set to
    MEEBOset.}
  
  \item{bgMethod}{Character string specifying which background
    correction method to use. See ?backgroundCorrect in package limma
    for more details.}
  
  \item{normMethod}{Character string specifying which normalization
    method should be used. See ?normalizeWithinArrays in package limma
    for more details.}
  
  \item{diagnosticPlot}{Logical. If 'TRUE', a quality diagnostic plot
    will be generated.}

  \item{output}{Logical. If 'TRUE', normalized M values and A values corresponding to
    the input GenePix files and additionnal quality measures are printed to a file.}

  \item{resdir}{A "character" string representing the directory where
    the results will be saved. By default, this is set to the current
    working directory ("."). }
  
  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}

  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding Oligo set annotations.}

  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{Additional arguments}
}

\details{
meeboQuality returns 3 types of quality control plots, specifically
designed for the various controls offered by the MEEBO set.
To assess the global performance of an hybridization, users can generate
a diagnostic plot summarizing several graphs and statictics by setting'diagnosticPlot=TRUE'. Then, the performance of the MEEBO set can beanalyzed by looking specifically at the mismatch and the tiling controls
('meeboSetQC=TRUE').

Finally, we are also providing several exploratory tools to assess the performance of the doping-controls
('DOPING=TRUE');these plots should be used only if a spike-in mixture
was added to the hybridization solution.
}

\value{
  meeboquality will produce several graphs, saved in an image file
  format. Please refer to the HTML description for more details. The
  function will also return the MAList object describing your tested slides.
  }

\author{Agnes Paquet}

\seealso{ \code{\link{meeboQualityPlots}}, \code{\link{gpQuality}}}

\examples{
if (interactive())
{
require(MEEBOdata)
datadir <- system.file("Meebo", package="MEEBOdata")
MA <- meeboQuality(fnames="RDI108_n.gpr",path=datadir,SpikeTypeFile="StanfordDCV1.7complete.txt",cy5col="CY5.ng._MjDC_V1.7",cy3col="CY3.ng._MjDC_V1.7",diagnosticPlot=TRUE,
DOPING=TRUE, meeboSetQC=TRUE, namecol="Name", resdir="MeeboQC", DEBUG=TRUE)
}
}

\keyword{dynamic}

