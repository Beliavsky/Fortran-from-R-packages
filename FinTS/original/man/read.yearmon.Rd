\name{read.yearmon}
\alias{read.yearmon}
\title{Reading Monthly zoo Series}
\description{
  Read a text file containing monthly data with a date column and return
  a zoo object with index = a yearmon series with the dates read as
  names.  
}
\usage{
read.yearmon(file, format = "", tz = "", FUN = NULL, regular = FALSE,
             index.column = 1, \dots)
}
\arguments{
  \item{file}{character giving the name of the file which the data
    are to be read from/written to. See \code{\link{read.table}} and
    \code{\link{write.table}} for more information.}
  \item{format}{date format argument passed to \code{\link{as.Date.character}}.}
  \item{tz}{time zone argument passed to \code{\link{as.POSIXct}}.}
  \item{FUN}{a function for computing the index from the first column
    of the data. See details.}
  \item{regular}{logical. Should the series be coerced to class \code{"zooreg"}
    (if the series is regular)?}
  \item{index.column}{integer. The column of the data frame in which the index/time
    is stored.}
  \item{\dots}{further arguments passed to \code{\link{read.table}} or
    \code{\link{write.table}}, respectively.}
}

\details{
  The arguments are as for \code{read.zoo}:  
  
  \code{TS <- read.zoo(...)}

  \code{zoo(coredata(TS), as.yearmon2(index(TS))) }
}

\value{
  an object of class \code{"zoo"} or \code{"zooreg"}.
}

\seealso{
  \code{\link[utils]{read.table}},
  \code{\link[zoo]{zoo}},
  \code{\link[zoo]{read.zoo}},
  \code{\link[zoo]{coredata}},
  \code{\link[zoo]{index}},
  \code{\link{as.yearmon2}}
}

\examples{
\dontrun{
## turn *numeric* first column into yearmon index
## where number is year + fraction of year represented by month
z <- read.zoo("foo.csv", sep = ",", FUN = as.yearmon2)
z2 <- read.yearmon("foo.csv", sep = ",")
}
}
\keyword{ts}
