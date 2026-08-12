\name{CKME}
\alias{CKME}
\title{Cassini Data Partitions Obtained by K-Means}
\description{
  A cluster ensemble of 50 \eqn{k}-means partitions of the Cassini data
  into three classes.
}
\usage{data("CKME")}
\format{
  A cluster ensemble of 50 (\eqn{k}-means) partitions.
}
\details{
  The ensemble was generated via
  \preformatted{
    require("clue")
    data("Cassini")
    set.seed(1234)
    CKME <- cl_boot(Cassini$x, 50, 3)
  }
}
\keyword{datasets}
