\name{nodeSpacer}
\alias{nodeSpacer}
\title{
Non-Uniform Finite Difference Node Spacer
}
\description{
\code{nodeSpacer} implements the spatial discretization scheme from Hout and Foulon (2010) with arbitrary left and right bounds. The function additionally includes logic for mesh shifting.
}
\usage{
nodeSpacer(K, leftBound, rightBound, nodes, density, kShift)
}
\arguments{
  \item{K}{
    scalar; option strike.
  }
  \item{leftBound}{
    scalar; near spatial boundary of the underlying domain.
  }
  \item{rightBound}{
    scalar; far spatial boundary of the underlying domain.
  }
  \item{nodes}{
    integer; number of nodes used in the spatial discretization.
  }
  \item{density}{
    scalar; impacts the concentration of nodes around the option strike. At 0, nodes are uniformly distributed between the \code{leftBound} and \code{rightBound}. Increasing the parameter increases the node concentration around the strike.
  }
  \item{kShift}{
    case; if 0, no mesh shifting, if 1, adjusts the node spacing such that the strike falls exactly between two nodes, if 2, adjusts the node spacing such that the strike falls exactly on a node. See Tavella and Randall (2000).
  }
}
\details{
Mesh shifting is accomplished by multiplying the vector of gridpoints by a scalar. For multi-asset options, this \code{nodeSpacer} is called iteratively to discretize each underlying's spatial domain.
}
\value{
\code{nodeSpacer} returns a vector of gridpoints used in spatial discretization in the finite difference method. The \code{nodes} input determines the length of the output vector.
}
\references{
Hout, K. J., Foulon, S., 2010. ADI finite difference schemes for option pricing in the Heston model with correlation. \emph{International Journal of Numerical Analysis and Modeling}, \bold{7} (2), 303--320.
\url{http://www.math.ualberta.ca/ijnam/Volume-7-2010/No-2-10/2010-02-06.pdf}

Pooley, D. M., Vetzal, K. R., Forsyth, P. A., 2002. Convergence remedies for non-smooth payoffs in option pricing.
\url{https://cs.uwaterloo.ca/~paforsyt/report.pdf}

Tavella, D., Randall, C., 2000. Pricing Financial Instruments: The Finite Difference Method. John Wiley & Sons, Inc., New York.
}
\author{
Michael Eichenberger and Carlo Rosa
}
\examples{
# sample mesh spacing
plot(nodeSpacer(100, 0, 500, 26, 5, 1), rep(0, times=26), main='Non-Uniform Mesh Spacing',
  xlab='Underlying Price (Strike = 100)', ylab='', yaxt='n', type='p', cex=0.5, pch=16)
}
