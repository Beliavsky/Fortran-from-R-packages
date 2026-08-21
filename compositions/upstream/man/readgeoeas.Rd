\name{Read standard data files}
\alias{read.geoeas}
\alias{read.geoEAS}
\title{Reads a data file in a geoeas format}
\description{
  Reads a data file, which must be formatted either as a
  geoEAS file (described below).
}
\usage{
read.geoeas(file)
read.geoEAS(file)
          }
\arguments{
  \item{file}{a file name, with a specific format}
}
\value{
  A data set, with a "title" attribute.
}
\details{
  The data files must be in the adequate format: "read.geoEAS" and
  "read.geoeas" read geoEAS format.

  The geoEAS format has the following structure:
  \describe{
    \item{1}{a first row with a description of the data set}
    \item{2}{the number of variables (=nvars)}
    \item{3}{"nvars"  rows, each containing the name of a variable}
    \item{4}{the data set, in a matrix of "nvars" columns
      and as many rows as individuals}
    }
}

\note{ Labels and title should not contain tabs.
       This might produce an error when reading.}
\author{Raimon Tolosana-Delgado}
\seealso{
  \code{\link{read.table}}
}

\references{
  Missing references
}
\examples{
#
# Files can be found in the test-subdirectory of the package
#
\dontrun{
  read.geoeas("TRUE.DAT")
  read.geoEAS("TRUE.DAT")
}
}
\keyword{IO}
