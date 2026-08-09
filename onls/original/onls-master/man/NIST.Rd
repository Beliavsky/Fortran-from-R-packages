\name{NIST}
\alias{NIST}
\encoding{latin1}

\title{Complicated nonlinear models from the NIST database}

\description{
Fits, by \code{\link{onls}}, a few complicated nonlinear models from the NIST database (\url{https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml}), \emph{i.e.} 'Chwirut2', 'Gauss1', 'ENSO' and 'Rat43'.
}

\usage{
NIST(verbose = FALSE) 
}

\arguments{
\item{verbose}{logical. If \code{TRUE}, steps from the \code{\link{onls}} output are printed to the console.}
}


\value{
Output of the models
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
\donttest{
NIST()
}
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

