\name{optimget}
\alias{optimget}
\title{
  Queries an optimization option list
}
\description{
  This function allows to make queries on an existing optimization option list.
  This list must have been created and updated by the \code{optimset} function.
  The \code{optimget} allows to retrieve the value associated with a given key.
}
\usage{
  optimget(options = NULL, key = NULL, value = NULL)
}
\arguments{
  \item{options}{A list created or modifies by \code{optimset}.}
  \item{key}{A single character string, which should be the name of the field in
    \code{options} to query (case insensitive).}
  \item{value}{A default value.}
}
\details{
  \code{key} is matched against the field names of \code{options} using
  \code{grep} and a case-insensitive regular expression. If \code{key} is not
  found in \code{options}, the function returns NULL. If several matches are
  found, \code{optimget} is stopped.
}
\value{
  Return \code{options$key} if \code{key} is found in \code{options}. Return
  \code{value}, otherwise.
}
\author{
  Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

  Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{optimset}}
}
\examples{
  opt <- optimset(method='fminsearch')
  optimget(opt,'Display')
  optimget(opt,'abc','!@')
}
\keyword{ method }
