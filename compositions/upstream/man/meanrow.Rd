\name{meanrow}
\alias{meanRow}
\alias{meanCol}
\alias{mean.row}
\alias{mean.col}
\title{The arithmetic mean of rows or columns}
\description{
  Computes the arithmetic mean.
}
\usage{
          meanRow(x,\dots, na.action=get(getOption("na.action")))
          meanCol(x,\dots, na.action=get(getOption("na.action")))
          }
\arguments{
  \item{x}{a numeric vector or matrix of data }
  \item{\dots}{arguments to \code{\link{mean}} }
  \item{na.action}{ The na.action to be used: one of
    \code{\link{na.omit}},\code{\link{na.fail}},\code{\link{na.pass}}}
}
\value{
  The arithmetic means of the rows
  (meanRow) or  columns (meanCol) of x.
}
\details{
  Computes the arithmetic means of the rows
  (meanRow) or  columns (meanCol) of x.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{mean.rplus}} }

\examples{
data(SimulatedAmounts)
meanCol(sa.tnormals)
meanRow(sa.tnormals)

}
\keyword{univar}
