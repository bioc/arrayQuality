\name{readcontrolCode}
\alias{readcontrolCode}

\title{Control status information}
 
\description{
  This component reads spot type information from a tab delimited text
  file to a matrix.}

\usage{
readcontrolCode(file = "SpotTypes.txt", path = NULL, sep = "\t", check.names = FALSE, controlId=c("ID", "Name"), ...) 
}

\arguments{
  
  \item{file}{Character string giving the name of the file specifying the
    spot types.}
  
  \item{path}{Character string giving the directory containing the file.
    Can be omitted if the file is in the current working directory.}
  
  \item{sep}{the field separator character.}
  
  \item{check.names}{Logical, if 'FALSE' column names will not be converted to
    valid variable names, for example spaces in column names will not be left as is.}

  \item{controlId}{Character string. Name of the column of the gpr
    file used to define controls.}

  \item{\dots}{additional arguments}
}

\value{
  A 2 column matrix named controlCode.
  }

\author{Jean Yee Hwa Yang, Agnes Paquet}

\seealso{\code{\link{readSpotTypes}}}

\keyword{hplot}


