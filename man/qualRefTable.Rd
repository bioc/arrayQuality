\name{qualRefTable}
\alias{qualRefTable}

\title{General hybridization quality score scaling table}

\description{This function helps you scale quality measures in order
  to compare them on the same plot. It is used on reference slides to
  create a look-up table, which will be used to scale other slides.}

\usage{qualRefTable(fnames = NULL, path = ".")
}

\arguments{

  \item{fnames}{A "character" string naming the input files.}
  
  \item{path}{A "character" string representing the data directory. By
    default this is set to the current working directory (".").}
}

\value{A matrix containing the mean and iqr for tested slides.}

\author{Agnes Paquet}

\seealso{\code{\link{gpQuality}}, \code{\link{globalQuality}},
  \code{\link{qualBoxplot}}, \code{\link{MmScalingTable}}
}

\keyword{programming}
