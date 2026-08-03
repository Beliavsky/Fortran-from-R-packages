\name{matrixFDM}
\alias{matrixFDM}
\title{
Finite Difference Matrix Generator
}
\description{
Generates a modified coefficient matrix (\bold{M}-matrix) used in the finite difference method from the option inputs. See Tavella and Randall (2000) for more on the standard formulation of the \bold{M}-matrix.
}
\usage{
matrixFDM(S, rf, q, vol, rho)
}
\arguments{
  \item{S}{
    list containing the vectors of spatial grid points associated with each underlying. Vector sizes of underlying spatial grid points need not be equal.
  }
  \item{rf}{
    scalar; applicable risk-free rate (domestic risk-free rate).
  }
  \item{q}{
    vector; holding costs of the option's underlyings (dividends, foreign risk-free rates, etc.).
  }
  \item{vol}{
    vector; volatilities of the option's underlyings.
  }
  \item{rho}{
    matrix; correlation matrix of the option's underlyings.
  }
}
\details{
\code{matrixFDM} first constructs the non-zero diagonals of the \bold{M}-matrix and stores them as columns.  The \code{bandSparse} function from the \pkg{Matrix} package then constructs a sparse banded matrix from the columns of the previously contructed matrix. Spatial domain boundaries are calculated first-order inwards with second difference terms dropped, maintaining block tridiagonality.
}
\value{
\code{matrixFDM} returns a \code{CsparseMatrix}-class matrix used for timestepping in the finite difference method.
}
\references{
Tavella, D., Randall, C., 2000. Pricing Financial Instruments: The Finite Difference Method. John Wiley & Sons, Inc., New York.
}
\author{
Michael Eichenberger and Carlo Rosa
}
\examples{
# finite difference matrix for uniformly-spaced two-asset option
S1 <- list(seq(0, 5, by=1), seq(0, 5, by=1))
rf <- 0.1
q <- c(0.05, 0.04)
vol <- c(0.20, 0.25)
rho <- matrix(c(1,-0.5,-0.5,1), 2, 2)
matrixFDM(S1, rf, q, vol, rho)
}
