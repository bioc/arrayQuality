\name{heeboQuality}
\alias{heeboQuality}

\title{Diagnostics plots designed for HEEBO set controls}
\description{
  This component generates several exploratory plots customized to the
  various types of controls provided in the HEEBO set. All results are
  saved as an image. Tested on GenePix format only (06-29-2006). For
  more details about the plots, please refer to the HTML description.
}
\usage{
heeboQuality(fnames = NULL, path = ".", galfile = NULL, source
="genepix.median", other.columns =
c("Flags"), controlMatrix=HeeboSpotTypes,controlId = c("ID", "Name"), 
DOPING = TRUE, heeboSetQC = TRUE, SpotTypeFile = NULL, SpikeTypeFile =
NULL, cy3col = "Cy3_ng", cy5col = "Cy5_ng", id = "SeqID", namecol =
c("Symbol", "Name"), annot = NULL, bgMethod = "none", normMethod = "p",
diagnosticPlot = TRUE, output = TRUE, resdir = ".", dev = "png",
organism = "Hs", DEBUG = FALSE, ...)
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
    HeeboSpotTypes.}
  
  \item{controlId}{Character string. Name of the column of the gpr
    file (or gal file) used to define controls.}
  
  \item{DOPING}{Logical. If 'TRUE', doping controls quality plots are
    generated.}
  
  \item{heeboSetQC}{Logical. If 'TRUE', mismatch and tiling controls
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
    used to look-up probes annotations. By default, it will be set to
    HEEBOset.}
  
  \item{bgMethod}{Character string specifying which background
    correction method to use. See ?backgroundCorrect in package limma
    for more details.}
  
  \item{normMethod}{Character string specifying which normalization
    method should be used. See ?normalizeWithinArrays in package limma
    for more details.}
  
  \item{diagnosticPlot}{Logical. If 'TRUE', a quality diagnostic plot
    will be generated.}

  \item{output}{Logical. If 'TRUE', normalized M values and A values
    corresponding to the input GenePix files and additionnal quality
    measures are printed to a file.}

  \item{resdir}{A "character" string representing the directory where
    the results will be saved. By default, this is set to the current
    working directory ("."). }
  
  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}

  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Hs". It is used to retrieve the corresponding Oligo set
    annotations.}

  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{Additional arguments}
}

\details{
heeboQuality returns 3 types of quality control plots, specifically
designed for the various controls offered by the HEEBO set.
To assess the global performance of an hybridization, users can generate
a diagnostic plot summarizing several graphs and statictics by setting'diagnosticPlot=TRUE'. Then, the performance of the HEEBO set can beanalyzed by looking specifically at the mismatch and the tiling controls
('heeboSetQC=TRUE').

Finally, we are also providing several exploratory tools to assess the performance of the doping-controls
('DOPING=TRUE');these plots should be used only if a spike-in mixture
was added to the hybridization solution.
}

\value{
  heeboquality will produce several graphs, saved in an image file
  format. Please refer to the HTML description for more details. The
  function will also return the MAList object describing your tested slides.
  }

\author{Agnes Paquet}

\seealso{ \code{\link{heeboQualityPlots}}, \code{\link{gpQuality}},
  \code{\link{meeboQuality}}}

\examples{
if (interactive())
{
require(HEEBOdata)
datadir <- system.file("Heebo", package="HEEBOdata")
MA <-
heeboQuality(fnames="63421.gpr",galfile="hoc.gal",path=datadir,SpikeTypeFile="DCV2.0June06.txt",cy5col="Cy5_ng",cy3col="Cy3_ng",diagnosticPlot=TRUE, 
DOPING=TRUE, heeboSetQC=TRUE, namecol="Name", resdir="HeeboQC", DEBUG=TRUE)
}
}

\keyword{dynamic}

