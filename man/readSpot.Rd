\name{readSpot}
\alias{readSpot}

\title{Extraction of measures from Spot (.spot) files.}

\description{This component reads a Spot file (.spot) and returns
  columns used for quality control.
}

\usage{
readSpot(fnames = NULL, path= ".", galfile=NULL,DEBUG=FALSE, skip = 0, sep ="\t",quote="\"",controlId="ID", ...)
}

\arguments{
  \item{fnames}{A "character" string naming the input file.}
  
  \item{path}{a "character" string representing the data directory.  By
    default this is set to the current working directory (".").}

  \item{galfile}{ A "character" string naming the galfile associated
    with the input file. If galfile = NULL, readSpot will use the first
    .gal file found in the working directory.}
  
  \item{DEBUG}{If 'TRUE', debug statements are printed.}

  \item{skip}{Number of lines to skip in the gpr files.}

  \item{sep}{A "character" string defining the type of separation for
    the columns in the gpr files.}

  \item{quote}{A "character" string defining the type of quote in the
    gpr files.}
  
  \item{controlId}{Character string. Name of the column of the Spot
    file used to define controls.}

  \item{\dots}{additional arguments.}
}

\value{A list of vectors containing information extracted from the
  Agilent file}

\author{Agnes Paquet}

\seealso{
  \code{\link{slideQuality}}, \code{\link{agQuality}}, \code{\link{globalQuality}}
}

\examples{
}

\keyword{programming}

