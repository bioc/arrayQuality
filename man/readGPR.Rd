\name{readGPR}
\alias{readGPR}

\title{Reading GenePix gpr file}

\description{This component reads a GenePix file (.gpr) and returns
  columns used for quality control.
}

\usage{
readGPR(fnames = NULL, path= ".", DEBUG=FALSE, skip = 0, sep ="\t",
quote="", ...)
}

\arguments{
  \item{fnames}{A "character" string naming the input file.}
  
  \item{path}{a "character" string representing the data directory.  By
    default this is set to the current working directory (".").}
  
  \item{DEBUG}{If 'TRUE', debug statements are printed.}

  \item{skip}{Number of lines to skip in the gpr files.}

  \item{sep}{A "character" string defining the type of separation for
    the columns in the gpr files.}

  \item{quote}{A "character" string defining the type of quote in the
    gpr files.}

  \item{\dots}{additional arguments.}
}

\value{A list of vectors containing information extracted from the
  GenePix file}

\author{Agnes Paquet}

\seealso{
  \code{\link{slideQuality}}, \code{\link{gpQuality}}, \code{\link{globalQuality}}
}

\examples{
#datadir <- system.file("data", package="arrayQuality")
#gprData <- readGPR(fnames="9Mm137.gpr", path=datadir)
}

\keyword{programming}

