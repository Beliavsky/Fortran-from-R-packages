\name{Brownian}
\alias{rBrownian}
\alias{deBrowning}
\title{Brownian and Related Motions}
\description{
  Simulate paths of dependent Brownian motions, geometric
  Brownian motions and Brownian bridges based on given increment
  copula samples. And extract copula increments from paths of
  dependent Brownian motions and geometric Brownian motions.
}
\usage{
rBrownian(N, t, d = 1, U = matrix(runif(N * n * d), ncol = d),
          drift = 0, vola = 1, type = c("BM", "GBM", "BB"), init = 1)
deBrowning(x, t, drift = 0, vola = 1, type = c("BM", "GBM"))
}
\arguments{
  \item{N}{number \eqn{N} of paths to simulate (positive integer).}
  \item{x}{\eqn{n+1}-vector containing one path of the specified
    stochastic process or \eqn{(n+1, d)}-matrix containing one path of
    the specified \eqn{d} stochastic processes or \eqn{(N, n+1, d)}-array
    containing \eqn{N} paths of the specified \eqn{d} stochastic processes.}
  \item{t}{\eqn{n+1}-vector of the form
    \eqn{(t_0,\dots,t_n)}{(t_0,...,t_n)}
    with \eqn{0 = t_0 < \dots < t_n}{0 = t_0 < ... < t_n}
    containing the time points where the stochastic processes are considered.}
  \item{d}{number \eqn{d} of stochastic processes to simulate (positive
    integer).}
  \item{U}{\eqn{(N\cdot n, d)}{(N*n, d)}-\code{\link{matrix}} of copula
    realizations to be converted to the joint increments of
    the stochastic processes.}
  \item{drift}{\eqn{d}-vector or number
    (then recycled to a \eqn{d}-vector) of drifts
    (typically denoted by \eqn{\mu}{mu}). Note that risk-neutral
    drifts are \eqn{r - \sigma^2/2}{r - sigma^2/2}, where \eqn{r}
    is the risk-free interest rate and \eqn{\sigma}{sigma}
    the volatility.}
  \item{vola}{\eqn{d}-vector or number (then recycled to a \eqn{d}-vector)
    of volatilities (typically denoted by \eqn{\sigma}{sigma}).}
  \item{type}{\code{\link{character}} string indicating whether
    a Brownian motion (\code{"BM"}), geometric Brownian motion
    (\code{"GBM"}) or Brownian bridge (\code{"BB"}) is to be considered.}
  \item{init}{\eqn{d}-vector or number (then recycled to a
    \eqn{d}-vector) of initial values (typically stock prices at time 0)
    for \code{type = "GBM"}.}
}
\value{
  \code{rBrownian()} returns an \eqn{(N, n+1, d)}-array containing
  the \eqn{N} paths of the specified \eqn{d}
  stochastic processes simulated at the \eqn{n+1} time points
  (\eqn{t_0 = 0, t_1,\dots, t_n}{t_0 = 0, t_1, ..., t_n}).

  \code{deBrowning()} returns an \eqn{(N, n, d)}-array containing
  the \eqn{N} paths of the copula increments of the \eqn{d} stochastic
  processes over the \eqn{n+1} time points
  (\eqn{t_0 = 0, t_1,\dots, t_n}{t_0 = 0, t_1, ..., t_n}).
}
\author{Marius Hofert}
\examples{
## Setup
d <- 3 # dimension
library(copula)
tcop <- tCopula(iTau(tCopula(), tau = 0.5), dim = d, df = 4) # t_4 copula
vola <- seq(0.05, 0.20, length.out = d) # volatilities sigma
r <- 0.01 # risk-free interest rate
drift <- r - vola^2/2 # marginal drifts
init <- seq(10, 100, length.out = d) # initial stock prices
N <- 100 # number of replications
n <- 25 # number of time intervals
t <- 0:n/n # time points 0 = t_0 < ... < t_n

## Simulate N paths of a cross-sectionally dependent d-dimensional
## (geometric) Brownian motion ((G)BM) over n time steps
set.seed(271)
U <- rCopula(N * n, copula = tcop) # for dependent increments
X <- rBrownian(N, t = t, d = d, U = U, drift = drift, vola = vola) # BM
S <- rBrownian(N, t = t, d = d, U = U, drift = drift, vola = vola,
               type = "GBM", init = init) # GBM
stopifnot(dim(X) == c(N, n+1, d), dim(S) == c(N, n+1, d))

## DeBrowning
Z.X <- deBrowning(X, t = t, drift = drift, vola = vola) # BM
Z.S <- deBrowning(S, t = t, drift = drift, vola = vola, type = "GBM") # GBM
stopifnot(dim(Z.X) == c(N, n, d), dim(Z.S) == c(N, n, d))
## Note that for BMs, one loses one observation as X_{t_0} = 0 (or some other
## fixed value, so there is no random increment there that can be deBrowned.

\donttest{
    ## If we map the increments back to their copula sample, do we indeed
    ## see the copula samples again?
    U.Z.X <- pnorm(Z.X) # map to copula sample
    U.Z.S <- pnorm(Z.S) # map to copula sample
    stopifnot(all.equal(U.Z.X, U.Z.S)) # sanity check
    ## Visual check
    pairs(U.Z.X[,1,], gap = 0) # check at the first time point of the BM
    pairs(U.Z.X[,n,], gap = 0) # check at the last time point of the BM
    pairs(U.Z.S[,1,], gap = 0) # check at the first time point of the GBM
    pairs(U.Z.S[,n,], gap = 0) # check at the last time point of the GBM
    ## Numerical check
    ## First convert the (N * n, d)-matrix U to an (N, n, d)-array but in
    ## the right way (array(U, dim = c(N, n, d)) would use the U's in the
    ## wrong order)
    U. <- aperm(array(U, dim = c(n, N, d)), perm = c(2,1,3))
    ## Now compare
    stopifnot(all.equal(U.Z.X, U., check.attributes = FALSE))
    stopifnot(all.equal(U.Z.S, U., check.attributes = FALSE))
}

\donttest{
    ## Generate dependent GBM sample paths with quasi-random numbers
    library(qrng)
    set.seed(271)
    U.. <- cCopula(to_array(sobol(N, d = d * n, randomize = "digital.shift"), f = n),
                   copula = tcop, inverse = TRUE)
    S. <- rBrownian(N, t = t, d = d, U = U.., drift = drift, vola = vola,
                    type = "GBM", init = init)
    pairs(S [,2,], gap = 0) # pseudo-samples at t_1
    pairs(S.[,2,], gap = 0) # quasi-samples at t_1
    pairs(S [,n+1,], gap = 0) # pseudo-samples at t_n
    pairs(S.[,n+1,], gap = 0) # quasi-samples at t_n
}

\donttest{
    ## Generate paths from a Brownian bridge
    B <- rBrownian(N, t = t, type = "BB")
    plot(NA, xlim = 0:1, ylim = range(B),
         xlab = "Time t", ylab = expression("Brownian bridge path"~(B[t])))
    for(i in 1:N)
        lines(t, B[i,,], col = adjustcolor("black", alpha.f = 25/N))
}
}
\keyword{datagen}