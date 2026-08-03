\name{multiAssetOption}
\alias{multiAssetOption}
\title{
Finite Difference Method for Multi-Asset Option Valuation
}
\description{
\code{multiAssetOption} generalizes the standard finite difference method to handle mulitple underlying assets, non-uniform grid spacing, non-uniform timestepping, and American exercise. The implementation allows users to vary the option setup (number of underlying assets, call vs. put, European vs. American, etc.) as well as the features of the numerical method (grid spacing, timestepping scheme, etc.). Strike shifting the mesh and Rannacher smoothing are optionally included to remedy problems arising from potential spurious oscillations in the solution.
}
\usage{
multiAssetOption(X)
}
\arguments{
  \item{X}{
    list of inputs. List items given in the \bold{Details} section.
}
}
\details{
Items of the input list \code{X} are as follows:

\describe{
  \item{\code{X$opt$nAsset}}{
    integer; number of underlying assets.
}
  \item{\code{X$opt$payType}}{
    case; if 0, digital payoff, if 1, best-of payoff, if 2, worst-of payoff.
}
  \item{\code{X$opt$exerType}}{
    case; if 0, European exercise, if 1, American exercise.
}
  \item{\code{X$opt$pcFlag}}{
    case vector; if 0, call, if 1, put.
}
  \item{\code{X$opt$ttm}}{
    scalar; time to maturity.
}
  \item{\code{X$opt$strike}}{
    vector; option strikes.
}
  \item{\code{X$opt$rf}}{
    scalar; applicable risk-free rate (domestic risk-free rate).
}
  \item{\code{X$opt$q}}{
    vector; holding costs of the option's underlyings (dividends, foreign risk-free rates, etc.).
}
  \item{\code{X$opt$vol}}{
    vector; volatilities of the option's underlyings.
}
  \item{\code{X$opt$rho}}{
    matrix; correlation matrix of the option's underlyings.
}
  \item{\code{X$fd$m}}{
    vector; number of spatial steps for each underlying's domain discretization.
}
  \item{\code{X$fd$leftBound}}{
    vector; near spatial boundaries of each underlying's domain.
}
  \item{\code{X$fd$kMult}}{
    vector; right boundary strike multiples. If 0, far domain boundary calculated via formula given in Kangro and Nicolaides (2000). Otherwise, far domain boundary calculated as the strike multiplied by the strike multiple.
}
  \item{\code{X$fd$density}}{
    vector; impacts the concentration of nodes around the option strike. At 0, nodes are uniformly distributed between the near and far boundaries. Increasing the parameter increases the node concentration around the strike.
}
  \item{\code{X$fd$kShift}}{
    case vector; if 0, no mesh shifting, if 1, adjusts the node spacing such that the strike falls exactly between two nodes, if 2, adjusts the node spacing such that the strike falls exactly on a node. See Tavella and Randall (2000).
}
  \item{\code{X$fd$theta}}{
    scalar; implicitness parameter of the theta method. Chosen between 0 (fully explicit) and 1 (fully implicit).
}
  \item{\code{X$fd$maxSmooth}}{
    integer; number of Rannacher smoothing steps. See Rannacher (1984).
}
  \item{\code{X$fd$tol}}{
    scalar; error tolerance in penalty iteration for American exercise.
}
  \item{\code{X$fd$maxIter}}{
    integer; maximum number of iterations per penalty loop for American exercise.
}
  \item{\code{X$time$tsType}}{
    case; if 0, constant timestep size, if 1, adaptive timestep size. See Forsyth and Vetzal (2002).
}
  \item{\code{X$time$N}}{
    integer; number of total timesteps if not using adaptive timesteps.
}
  \item{\code{X$time$dtInit}}{
    scalar; inital timestep size for adaptive timesteps.
}
  \item{\code{X$time$dNorm}}{
    scalar; target relative change for adaptive timesteps.
}
  \item{\code{X$time$D}}{
    scalar; normalizing parameter for adaptive timesteps.
}
}

The classical order for the state vectors output from the function is illustrated by example. With two underlying assets, option values in each state vector are stored in the order: [11, 21, 31, ... , M1, 12, 22, ... , MN] with M being the total number of nodes used in the first asset spatial discretization and N being the total number of nodes in the second.
}
\value{
\code{multiAssetOption} returns a list:
  \item{value}{matrix of per-unit option values. Each column stores the state of the option value array (collection of option values for all nodes of the spatial mesh) as a vector following the classical order (see \bold{Details} section). The columns of the matrix are indexed over time, with the first column beginning at option maturity, and subsequent columns moving backward in time.}
  \item{S}{list containing the vectors of spatial grid points associated with each underlying. Vector sizes of underlying spatial grid points need not be equal.}
  \item{dimS}{dimension of option value array. This item can be used to reshape the column vectors in \code{value} into an appropriately dimensioned array using \code{array( ...  , dim=dimS)}.}
    \item{time}{vector of times associated with each column of the \code{value} item.}
  
For each column (time) of the \code{value} item, the option value at that time can be calculated as the option's notional amount multiplied by the unit option value interpolated over the \code{S} item at the current underlying prices.
}
\references{
Forsyth, P.A., Vetzal, K.R., 2002. Quadratic convergence for valuing American options using a penalty method. \emph{SIAM Journal on Scientific Computing}, \bold{23} (6), 2095--2122.

Kangro, R., Nicolaides, R., 2000. Far field boundary conditions for Black-Scholes equations. \emph{SIAM Journal on Numerical Analysis}, \bold{38} (4), 1357--1368.

Rannacher, R., 1984. Finite element solution of diffusion problems with irregular data. \emph{Numberische Mathematik}, \bold{43}, 309--327.

Tavella, D., Randall, C., 2000. Pricing Financial Instruments: The Finite Difference Method. John Wiley & Sons, Inc., New York.
}
\author{
Michael Eichenberger and Carlo Rosa
}
\examples{
# european dual-asset digital option example

# initialize inputs list
X <- list()

# option inputs
X$opt$nAsset <- 2
X$opt$payType <- 0
X$opt$exerType <- 0
X$opt$pcFlag <- c(1, 0)
X$opt$ttm <- 0.5
X$opt$strike <- c(110, 90)
X$opt$rf <- 0.10
X$opt$q <- c(0.05, 0.04)
X$opt$vol <- c(0.20, 0.25)
X$opt$rho <- matrix(c(1, -0.5, -0.5, 1), X$opt$nAsset, X$opt$nAsset)

# finite difference inputs
X$fd$m <- c(20, 20)
X$fd$leftBound <- c(0, 0)
X$fd$kMult <- c(0, 0)
X$fd$density <- c(5, 5)
X$fd$kShift <- c(1, 1)
X$fd$theta <- 0.5
X$fd$maxSmooth <- 2
X$fd$tol <- 1e-7
X$fd$maxIter <- 3

# timestep inputs
X$time$tsType <- 0
X$time$N <- min(X$fd$m) * 4
X$time$dtInit <- 0.1 / 4^log2(min(X$fd$m)/5)
X$time$dNorm <- 5 / 2^log2(min(X$fd$m)/5)
X$time$D <- 0.05

# function check
output <- multiAssetOption(X)
}