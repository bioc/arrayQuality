\name{qualityScore}
\alias{qualityScore}

\title{Quality Control score for general hybridization}

\description{This function returns, for each quality measure, the
  percentage of reference slides measures which are below the tested
  slide values. For more details on this score, refer to the online
  manual.
}

\usage{qualityScore(slidequality, organism = c("Mm", "Hs"), reference =
NULL)
}

\arguments{
  \item{slidequality}{A quality matrix from slideQuality.}

  \item{reference}{A matrix resulting from globalQuality, to be used as
    reference table to compare slides. If 'NULL', the default table
    corresponding to organism will be used. See details for more
    information.}
  
  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding reference tables.}
}

\value{A matrix of percentages.}

\author{Agnes Paquet}

\keyword{programming}

