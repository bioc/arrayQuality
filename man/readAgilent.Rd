\name{readAgilent}
\alias{readAgilent}

\title{Extraction of measures from Agilent (.txt) files.}

\description{This component reads an Agilent file (.txt) and returns
  columns used for quality control.
}

\usage{
readAgilent(fnames = NULL, path= ".", DEBUG=FALSE, skip = 0, sep ="\t",
quote="\"",controlId=c("ProbeName"), ...)
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
  
  \item{controlId}{Character string. Name of the column of the 
    file used to define controls.}

  \item{\dots}{additional arguments.}
}

\value{A list of vectors containing information extracted from the
  Agilent file}

\author{Agnes Paquet}

\seealso{
  \code{\link{slideQuality}}, \code{\link{agQuality}}, \code{\link{globalQuality}}
}


\keyword{programming}

