\name{scaleRefTable}
\alias{scaleRefTable}

\title{General hybridization quality scaling}

\description{This function helps you scale quality measures in order
  to compare them on the same plot. It is used on reference slides to
  create a look-up table, which will be used to scale other slides.}

\usage{scaleRefTable(reference=NULL, organism=c("Mm", "Hs"))}
}

\arguments{

  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding reference
    tables.}

  \item{reference}{A matrix resulting from globalQuality, to be used as
    reference table to compare slides. If 'NULL', the default table
    corresponding to organism will be used.}
}

\value{A matrix containing median and iqr for each quality measure for
  tested slides.}

\author{Agnes Paquet}

\seealso{\code{\link{gpQuality}}, \code{\link{globalQuality}},
  \code{\link{qualBoxplot}}}
}

\keyword{programming}
