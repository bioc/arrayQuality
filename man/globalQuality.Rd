\name{globalQuality}
\alias{globalQuality}

\title{Quality measures for general hybridization.}

\description{
  This function provides Quality Control measures for GenePix
  format files. It is used to create a table of measures to be used as
  reference in gpQuality.
}

\usage{
globalQuality(fnames = NULL, path = ".", organism = c("Mm", "Hs"),
output = FALSE, resdir = ".", DEBUG = FALSE,...)
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

  \item{\dots}{additional arguments}
}

\value{
  A matrix of Quality Control measures, each column representing a
  different input slide. 
}

\author{Agnes Paquet}

\seealso{\code{\link{gpQuality}}, \code{\link{slideQuality}},
  \code{\link{MmReferenceDB}}
}

\examples{
datadir <- system.file("data", package="arrayQuality")
reference <- globalQuality(fnames="9Mm137.gpr", path=datadir, organism="Mm")
}

\keyword{programming}






