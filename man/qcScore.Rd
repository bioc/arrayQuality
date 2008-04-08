\name{qcScore}
\alias{qcScore}

\title{Quality Control score for general hybridization}

\description{This function returns, for each quality measure, the
  number of qc measures of the tested slides which are below the
  reference slides boundaries. 
}

\usage{qcScore(arrayQuality,reference)
}

\arguments{

  \item{arrayQuality}{Matrix of quality results from \code{slideQuality}
    or \code{globalQuality}.}

  \item{reference}{A matrix resulting from globalQuality, to be used as
    reference table to compare slides. If 'NULL', the default table
    corresponding to organism will be used. See details for more
    information.}
  
}

\value{qcScore returns the number of qc measures under the
  lower limits of the reference values.}

\author{Agnes Paquet}

\keyword{programming}

