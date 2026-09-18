\name{banknote}
\alias{banknote}
\docType{data}

\title{Swiss Banknotes Data}

\description{
  Data from Tables 1.1 and 1.2 (pp. 5-8) of Flury  and Riedwyl
  (1988). There are  six measurements made on 200 Swiss
  banknotes  (the old-Swiss 1000-franc). The banknotes belong to two
  classes of equal size:  \emph{genuine} and \emph{counterfeit}.
}

\usage{data(banknote)}

\format{A  \code{data.frame} of dimension \code{200x7} with the following variables:

  \describe{
    \item{Class}{a \code{factor} with classes: \code{genuine}, \code{counterfeit}}
    \item{Length}{Length of bill (mm)}
    \item{Left}{Width of left edge (mm)}
    \item{Right}{Width of right edge (mm)}
    \item{Bottom}{Bottom margin width (mm)}
    \item{Top}{Top margin width (mm)}
    \item{Diagonal}{Length of diagonal (mm)}
  }
}

\source{
  Flury, B. and Riedwyl, H. (1988).
  \emph{Multivariate Statistics: A practical approach.}
  London: Chapman & Hall.
}

\keyword{datasets}
