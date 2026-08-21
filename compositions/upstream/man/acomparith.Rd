\name{acomparith}
\alias{power.acomp}
\alias{*.acomp}
\alias{/.acomp}
\title{Power transform in the simplex}
\description{
  The Aitchison Simplex with its two operations perturbation as +  and
  power transform as * is a vector space. This vector space is
  represented by these operations. 
}
\usage{
power.acomp(x,s)
## Methods for class "acomp"
## x*y
## x/y
          }
\arguments{
  \item{x}{an acomp composition or dataset of compositions (or a number
    or a numeric vector)}
  \item{y}{a numeric vector of size 1 or nrow(x)}
  \item{s}{a numeric vector of size 1 or nrow(x)}
}
\value{
  An \code{"acomp"} vector or matrix.
}
\details{
  The power transform is the basic multiplication operation of the
  Aitchison simplex seen as a vector space. It is defined as:
  \deqn{(x*y)_i:=  clo( (x_i^{y_i})_i )_i }
  The division operation is just the multiplication with \eqn{1/y}.
}
\note{For \code{*} the arguments x and y can be exchanged. Note that 
this definition generalizes the power by a scalar, since \code{y} or 
\code{s} may be given as a scalar, or as a vector with as many components as
the composition in \code{\link{acomp}} \code{x}. The result is then a matrix 
where each row corresponds to the composition powered by one of the scalars 
in the vector.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilr}},\code{\link{clr}}, \code{\link{alr}},
}

\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

}
\examples{
acomp(1:5)* -1 + acomp(1:5)
data(SimulatedAmounts)
cdata <- acomp(sa.lognormals)
plot( tmp <- (cdata-mean(cdata))/msd(cdata) )
class(tmp)
mean(tmp)
msd(tmp)
var(tmp)
}
\keyword{multivariate}
