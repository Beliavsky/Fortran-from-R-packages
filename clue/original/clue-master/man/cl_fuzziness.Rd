\name{cl_fuzziness}
\alias{cl_fuzziness}
\title{Partition Fuzziness}
\description{
  Compute the fuzziness of partitions.
}
\usage{
cl_fuzziness(x, method = NULL, normalize = TRUE)
}
\arguments{
  \item{x}{a cluster ensemble of partitions, or an R object coercible to
    such.}
  \item{method}{a character string indicating the fuzziness measure to
    be employed, or \code{NULL} (default), or a function to be taken as
    a user-defined method.  Currently available built-in methods are
    \code{"PC"} (Partition Coefficient) and \code{"PE"} (Partition
    Entropy), with the default corresponding to the first one.  If
    \code{method} is a character string, its lower-cased version is
    matched against the lower-cased names of the available built-in
    methods using \code{\link{pmatch}}.}
  \item{normalize}{a logical indicating whether the fuzziness measure
    should be normalized in a way that hard partitions have value 0, and
    \dQuote{completely fuzzy} partitions (where for all objects, all
    classes get the same membership) have value 1.}
}
\details{
  If \eqn{m} contains the membership values of a partition, the
  (unnormalized) Partition Coefficient and Partition Entropy are given
  by \eqn{\sum_{n,i} m_{n,i}^2} and \eqn{\sum_{n,i} H(m_{n,i})},
  respectively, where \eqn{H(u) = u \log u - (1-u) \log(1-u)}{u log(u)
    - (1-u) log(1-u)}.

  Note that the normalization used here is different from the
  normalizations typically found in the literature.

  If a user-defined fuzziness method is to be employed, is must be a
  function taking a matrix of membership values and a logical to
  indicate whether normalization is to be performed as its arguments (in
  that order; argument names are not used).
}
\value{
  An object of class \code{"cl_fuzziness"} giving the fuzziness
  values.
}
\references{
  \bibshow{Bezdek:1981}
}
\seealso{
  Function \code{\link[e1071]{fclustIndex}} in package \pkg{e1071},
  which also computes several other \dQuote{fuzzy cluster indexes}
  (typically based on more information than just the membership
  values).
}
\examples{
if(require("e1071", quietly = TRUE)) {
    ## Use an on-line version of fuzzy c-means from package e1071 if
    ## available.
    data("Cassini")
    pens <- cl_boot(Cassini$x, B = 15, k = 3, algorithm = "cmeans",
                    parameters = list(method = "ufcl"))
    pens
    summary(cl_fuzziness(pens, "PC"))
    summary(cl_fuzziness(pens, "PE"))
}
}
\keyword{cluster}
