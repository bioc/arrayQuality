\name{slideQuality}
\alias{slideQuality}

\title{Quality Control statistics for general hybridization}

\description{This component uses data extracted from GenePix file to
  provide quality control statistics.}

\usage{
slideQuality(gprData = NULL, DEBUG = FALSE,...)
}

\arguments{
  \item{gprData}{A list of vector, results from \code{readGPR},
    containing information extracted from the gpr file.}

  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{additional arguments}
}

\value{
  A matrix of numbers.}

\author{Agnes Paquet}

\seealso{\code{\link{gpQuality}}, \code{\link{globalQuality}},
  \code{\link{readGPR}}
}

\examples{
#datadir <- system.file("data", package="arrayQuality")
#gprdata <- readGPR(fnames="9Mm137.gpr", path=datadir)
#results <- slideQuality(gprdata)
}

\keyword{programming}


