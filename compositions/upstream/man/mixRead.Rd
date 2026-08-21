\name{mix.Read}
\alias{mix.Read}
\title{Reads a data file in a mixR format}
\description{
  Reads a data file, which is formatted in a simple
  compositional file including the first row with title, the second
  with data labels and afterwards the matrix with the data itself.
  In the first column of the matrix are cases labels. This is the format
  used in the mixR package.
}
\usage{
mix.Read(file,eps=1e-6)
          }
\arguments{
  \item{file}{a file name}
  \item{eps}{the epsilon to be used for checking null values.}
}
\value{
  A mixture object as a data frame with a title, row total, if constant, status
  (-2, -1, 0, 1, 2 or 3 -- see above) and class attributes and the data matrix.
}
\details{
  The data files must have the adequate structure:
  \describe{
    \item{1}{the first row with a title of the data set,}
    \item{2}{the second row with variables names}
    \item{3}{the data set in a matrix, rows as cases, variables in columns
      with the firs colum comprising cases labels.}
    }
  A mixture object 'm' consists of  \code{m$tit}  the title, \code{m$mat} the matrix with the data, 
  \code{m$sum} the value of the rows total, if constant and  \code{m$sta} the status of the mixture object
  with values:
  \tabular{rll}{
  \tab   -2   \tab     - matrix contains negative elements,\cr
  \tab   -1   \tab     - zero row sum exists,\cr
  \tab    0   \tab     - matrix contains zero elements,\cr
  \tab    1   \tab     - matrix contains positive elements, rows with different row sum(s), \cr
  \tab    2   \tab     - matrix with constant row sum   and\cr
  \tab    3   \tab     - closed mixture, the row sums are all equal to 1.
  }
}

\seealso{
  \code{\link{read.geoeas}}
  \code{\link{read.geoEAS}}
  \code{\link{read.table}}
}
\examples{
\dontrun{
  mix.Read("GLACIAL.DAT")
  mix.Read("ACTIVITY.DAT")
}
}
\keyword{IO}
