\name{as.yearmon2}
\alias{as.yearmon2}

\title{
  Conditionally convert x to yearmon if the conversion is unique,
  retaining x as names.
}
\description{
  Convert \code{x} to class \code{"yearmon"}.  If duplicate months are
  found, return \code{x}.  Otherwise, return the conversion with names =
  \code{x}.
}
\usage{
as.yearmon2(x, \dots)
}
\arguments{
  \item{x}{
    object suitable for \code{as.yearmon}
  }
  \item{\dots}{
    additional argument(s) (e.g., a format) passed to
    \code{as.yearmon}.
  }
}
\details{
  Dates for some monthly data include the day of the month on which the
  data were published.  For many purposes, one would like to have the
  data as a \code{zoo} object with a \code{yearmon} index, while still
  retaining the full date for other purposes.

  If the \code{yearmon} form of the input is not unique,
  \code{as.yearmon2} returns the input unchanged with a warning.
  Otherwise, it returns the \code{yearmon} conversion with the input as
  names.

}
\value{
  Returns either its argument or its argument converted to class
  \code{yearmon} with names.
}

\seealso{
  \code{\link[zoo]{yearmon}}
}
\examples{
x1 <- as.Date(c("2000-01-01", "2000-01-01"))
as.yearmon2(x1)
#Warning message:
#In as.yearmon2(x1) :
#  1 duplicate months found in 'x';  returning 'x' unchanged

x2 <- as.Date(c("2000-01-01", "2000-02-01"))
as.yearmon2(x2)
# month of x2 with names x2
}
\keyword{ts}
