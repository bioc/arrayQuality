\name{gpQuality}
\alias{gpQuality}

\title{Diagnostic plots and comparative boxplots for general hybridization}
 
\description{
  This component provides qualitative diagnostic plots and quantitative
  measures for assessing general hybridization quality. All results are
  displayed in a HTML report. GenePix format only.
}

\usage{

gpQuality(fnames = NULL, path = ".", organism = c("Mm", "Hs"),
compBoxplot = TRUE, reference = NULL, controlMatrix = controlCode,
controlId = c("ID", "Name"), output = FALSE, resdir = ".", dev= "png", DEBUG = FALSE,...)
}

\arguments{
  
  \item{fnames}{A "character" string naming the input files.}
  
  \item{path}{A "character" string representing the data directory.  By
    default this is set to the current working directory (".").}
  
  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding reference
    tables.}

  \item{compBoxplot}{Logical. If set to 'FALSE', only qualitative
    diagnostic plots will be plotted. gpQuality ouput will be limited to
    a diagnostic plot by gpr file and a marrayRaw object.}

  \item{reference}{A matrix resulting from globalQuality, to be used as
    reference table to compare slides. If 'NULL', the default table
    corresponding to organism will be used. See details for more information.}
  
 
  \item{controlMatrix}{A character matrix of n by 2 columns. The first
    column contains a few regular expression of spotted probe sequences
    and the second column contains the corresponding control status. By
    default, controlMatrix is set to controlCode.}

  \item{controlId}{Character string. Name of the column of the gpr
    file used to define controls.}

  \item{output}{Logical. If 'TRUE', normalized M values corresponding to
    the input GenePix files and quality measures are printed to a file.}

  \item{resdir}{A "character" string representing the directory where
    the results will be saved. By default, this is set to the current
    working directory ("."). }

  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}
  
  \item{DEBUG}{If 'TRUE', debug statements are printed.}

  \item{\dots}{additional arguments}
}


\details{
gpQuality returns 2 plots for each GenePix files passed as argument. The
first one is a qualitative diagnostic plot, a quick visual way to assess
slide quality. The second one is a comparative boxplot: each quality
control measure is compared to the range obtained for a database of 'good'
slides used as reference. You can use your own set of references created
using \code{globalQuality} and \code{qualRefTable} passed
in the arguments "reference" and "scalingTable", or use
the reference QC values stored in the datasets \code{MmReferenceDB} and
\code{MmScalingTable}. All results and quality scores are gathered in a
HTML report. For more details about the QC measures and the plots, please refer to the online manual.
}

\value{
  A list of 2 elements:

  \item{mraw}{A marrayRaw object created from the input files.}
  \item{quality}{A matrix containing Quality Control measures for all slides.}
}

\author{Agnes Paquet}

\seealso{\code{\link{globalQuality}}, \code{\link{qualBoxplot}},
  \code{\link{scaleRefTable}}
}

\examples{
  datadir <- system.file("gprQCData", package="arrayQuality")
  results <- gpQuality(fnames="9Mm137.gpr", path=datadir, organism="Mm")
}

\keyword{hplot}


