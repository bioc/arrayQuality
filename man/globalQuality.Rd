\name{globalQuality}
\alias{globalQuality}

\title{Quality measures for general hybridization.}

\description{
  This function provides Quality Control measures for GenePix, Spot and
  Agilent format files. It is used to create a table of measures to be used as
  reference in gpQuality, spotQuality or agQuality.
}

\usage{
globalQuality(fnames = NULL, path = ".", organism = c("Mm", "Hs"),
output = FALSE, resdir = ".", DEBUG = FALSE, inputsource = "readGPR",...)
}

\arguments{
  \item{fnames}{A "character" string naming the input files.}
  
  \item{path}{a "character" string representing the data directory.  By
    default this is set to the current working directory (".").}
  
  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding reference tables.}

  \item{output}{Logical. If 'TRUE', the quality measures are printed to
    a file.
  }

  \item{resdir}{A "character" string representing the directory where
    the results will be saved. By default, this is set to the current
    working directory ("."). }

  \item{DEBUG}{If 'TRUE', debug statements are printed.}

  \item{inputsource}{A "character" string providing the name of the function to
    use to read the input files. It should be inputsource = "readGPR"
    for GenePix format files, inputsource = "readSpot" for Spot files, or
    inputsource = "readAgilent" for Agilent format. By default,
    'inputsource' is set to "readGPR".}
  
  \item{\dots}{additional arguments}
}

\value{
  A matrix of Quality Control measures, each column representing a
  different input slide. 
}

\author{Agnes Paquet}

\seealso{\code{\link{gpQuality}}, \code{\link{slideQuality}},
  \code{\link{MmReferenceDB}}, \code{\link{readGPR}},
  \code{\link{readSpot}}, \code{\link{readAgilent}}
}

\examples{
datadir <- system.file("data", package="arrayQuality")
reference <- globalQuality(fnames="9Mm137.gpr", path=datadir, organism="Mm")
}

\keyword{programming}






