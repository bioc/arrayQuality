\name{PRvQCHyb}
\alias{PRvQCHyb}

\title{Print run Quality version Quality Control Hybridization}
\description{
    Qualitative diagnostic plots looking at print-run quality.  This
  component examine the QC hybridizations.
}

\usage{
PRvQCHyb(fnames, path=".", dev = "png", DEBUG=FALSE, prargs=NULL, samepr=TRUE, prname="xMm", save = TRUE, col,...)
}

\arguments{
  \item{fnames}{A "character" string naming the input files.}
  \item{path}{a character string representing the data directory. By
    default this is set to the current working directory (".").}
  \item{dev}{A "character" string naming the graphics device. This will
    take arguments "png", "jpeg" and "ps" only.}
  \item{DEBUG}{If 'TRUE', debug statements are printed.}
  \item{prargs}{A list with 4 components: Block, Row, Column and
    ncolumns. See Details for more information.}
  \item{samepr}{If 'TRUE', we assume everything in the directory are
    from the same print-run.}
  \item{prname}{A "character" string giving the name of the print-run.}
  \item{save}{If `TRUE', the figures will be saved to files.}
  \item{col}{color code for different control samples.}
  \item{\dots}{additional arguments.}
}

\details{
  The argument "prargs" is used to calculate the layout information
  about a print-run.  Components Block, Row and Columns denote the
  column names from the input data representing the print-tip location.
  The component "ncolumns" is an integer representing the number of
  print-tip columns in the data.  If the arguement is set to NULL,
  the the following default will be used.  list(Block="Block",
  Row="Row", Column="Column", ncolumns=4)
}

\value{
  Files of diagnostic plots.
}

\author{Jean Yee Hwa Yang}

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
   datadir <- system.file("data", package="arrayQuality")
   PRvQCHyb(fnames="9Mm137.gpr", path=datadir, prname="9Mm")
}
\keyword{hplot}% at least one, from doc/KEYWORDS
