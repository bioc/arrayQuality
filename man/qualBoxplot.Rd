\name{qualBoxplot}
\alias{qualBoxplot}

\title{Comparative boxplot for general hybridization Quality Control}

\description{This functions allows you to graphically compare your slide
  quality measures to results obtained for a database of 'good quality'
  slides.}

\usage{
qualBoxplot(arrayQuality = NULL, reference = NULL, scalingTable = NULL,
organism = c("Mm", "Hs"), ...)
}

\arguments{
  \item{arrayQuality}{Matrix of quality results from \code{slideQuality}
    or \code{globalQuality}.}

  \item{reference}{A matrix resulting from globalQuality, to be used as
    reference table to compare slides. If 'NULL', the default table
    corresponding to "organism" will be used. See details for more information.}
  
  \item{scalingTable}{A matrix resulting from qualRefTable, used to
    scale quality measures before plotting them. If 'NULL',the default table
    corresponding to organism will be used. See details for more information.}
  \item{organism}{A "character" string naming the organism genome
    printed on the array, either "Mm" or "Hs". By default, organism is
    set to "Mm". It is used to retrieve the corresponding reference
    tables.}

  \item{\dots}{additional arguments}

}

\details{
  You can use your own set of references created
  using \code{globalQuality} and \code{qualRefTable} passed
  in the arguments "reference" and "scalingTable", or use
  the reference QC values stored in the datasets \code{MmReferenceDB} and
  \code{MmScalingTable}.
}
 
\value{
  Returns a score vector containing, for each column in "arrayQuality",
  the number of quality measures below the range of "reference". The
  last element of the score vector is the total number of quality
  measures tested.
}

\author{Agnes Paquet}

\seealso{\code{\link{globalQuality}}, \code{\link{qualBoxplot}}, \code{\link{qualRefTable}}}

\examples{
#datadir <- system.file("data", package="arrayQuality")
#gprData <- readGPR(fnames="9Mm137.gpr", path=datadir)
#arrayQuality <- slideQuality(gprData, organism="Mm")
#qualBoxplot(arrayQuality)
}

\keyword{hplot}
