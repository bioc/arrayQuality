\name{maQualityPlots}
\alias{maQualityPlots}

\title{Qualitative diagnostic plot for general hybridization}

\description{
  This function generates diagnostic plots for a qualitative assessment of
  slide quality.  
}

\usage{
maQualityPlots(mrawObj, headerInfo = "", save = TRUE, dev = "png", col,
badspotfunction, DEBUG = FALSE, ...)
}

\arguments{
  \item{mrawObj}{\code{marrayRaw} object representing the slides to be
    tested.}

  \item{headerInfo}{Text to be used as header in the diagnostic plot.}

  \item{save}{Logical. If 'TRUE', the plot is saved to a file.}

  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}

  \item{col}{Vector of colors to use to describe different
    categories of spots.}

  \item{badspotfunction}{Function to use for bad spots on the array.}

  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  
  \item{\dots}{additional arguments}
}

\author{Jean Yang}

\examples{
# Example uses swirl dataset...

}

\keyword{hplot}

