\name{gmp-ifiworkarounds}
\alias{outer}
\title{Base Functions in 'gmp'-ified Versions}
\description{
  Functions from \pkg{base} etc which need a \emph{copy} in the
  \pkg{gmp} namespace so they correctly dispatch.
}
\usage{
outer(X, Y, FUN = "*", ...)
}
\arguments{
  \item{X, Y, FUN, ...}{See \pkg{base} package help: \code{\link[base]{outer}}.}
}
\seealso{
 \code{\link[base]{outer}}.
}
\examples{
twop <- as.bigz(2)^(99:103)
(mtw <- outer(twop, 0:2))
stopifnot(
   identical(dim(mtw), as.integer(c(5,3)))
 ,
   mtw[,1] == 0
 ,
   identical(as.vector(mtw[,2]), twop)
)
}
\keyword{misc}
