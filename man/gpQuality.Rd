\name{gpQuality}
\alias{gpQuality}
\alias{slideQuality}
\alias{readGPR}
%\alias{qualBoxplot}
%\alias{maQualityPlots}
%\alias{quality2HTML}

\title{Quality score and diagnostic plots for general hybridization}
 
\description{
  This component provides qualitative diagnostic plots and quantitative
  measures for assessing general hybridization quality. GenePix files
  format only.
}

\usage{
readGPR(fnames = NULL, path= ".", DEBUG=TRUE, skip = 0,
                     sep ="\t", quote="",...)

slideQuality(gprData=NULL,output = TRUE, DEBUG=TRUE,...)

gpQuality(fnames = NULL, path = ".", organism=c("Mm", "Hs"),
output=TRUE, plot=TRUE,resdir=".", dev="png", DEBUG = TRUE,...)
}

\arguments{
  
  \item{fnames}{A "character" string naming the input files.}
  \item{path}{a "character" string representing the data directory. By
    default this is set to the current working directory (".").}
  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only. By default, dev is set
    to "png".}
  \item{gprData}{A list of matrices containing information extracted
    from the input files.}
  \item{organism}{A "character" string defining the organism of the
    spotted probes. The default value is set to "Mm". See details for
    more information.}
  \item{resdir}{A "character" string representing the directory where
    the results will be saved. By default, this is set to the current
    working directory ("."). }
  \item{output}{Logical. If 'TRUE', the results are written to a file.}
  \item{plot}{Logical. If 'TRUE', qualitative diagnostics plots and
    quantitative quality plots are created and saved in "resdir".}
  \item{skip}{Number of lines to skip in the gpr files}
  \item{sep}{A "character" string defining the type of separation for
    the columns in the gpr files.}
  \item{quote}{A "character" string defining the type of quote in the
    gpr files.}
  \item{\dots}{additional arguments}
}

\details{
}

\value{
}

\author{Agnes Paquet}

\keyword{plot}