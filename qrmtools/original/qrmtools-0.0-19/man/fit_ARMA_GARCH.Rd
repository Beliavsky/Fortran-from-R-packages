\name{fit_ARMA_GARCH}
\alias{fit_ARMA_GARCH}
\title{Fitting ARMA-GARCH Processes}
\description{
  Fail-safe componentwise fitting of univariate ARMA-GARCH processes.
}
\usage{
fit_ARMA_GARCH(x, ugarchspec.list = ugarchspec(), solver = "hybrid",
               verbose = FALSE, ...)
}
\arguments{
  \item{x}{\code{\link{matrix}}-like data structure, possibly an
    \code{xts} object.}
  \item{ugarchspec.list}{object of class \code{uGARCHspec} (as returned by
    \code{ugarchspec()}) or a list of such. In case of a list, its
    length has to be equal to the number of columns of
    \code{x}. \code{ugarchspec.list} provides the ARMA-GARCH
    specifications for each of the time series (columns of \code{x}).}
  \item{solver}{string indicating the solver used; see \code{?ugarchfit}.}
  \item{verbose}{\code{\link{logical}} indicating whether verbose
    output is given.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link[rugarch]{ugarchfit}()}.}
}
\value{
  If \code{x} consists of one column only (e.g. a vector),
  \code{ARMA_GARCH()} returns the fitted object; otherwise it returns
  a list of such.
}
\author{Marius Hofert}
\seealso{
  \code{\link{fit_GARCH_11}()} for fast(er) and numerically more
  robust fitting of GARCH(1,1) processes.
}
\examples{
library(rugarch)
library(copula)

## Read the data, build -log-returns
data(SMI.12) # Swiss Market Index data
stocks <- c("CSGN", "BAER", "UBSN", "SREN", "ZURN") # components we work with
x <- SMI.12[, stocks]
X <- -returns(x)
n <- nrow(X)
d <- ncol(X)

## Fit ARMA-GARCH models to the -log-returns
## Note: - Our choice here is purely for demonstration purposes.
##         The models are not necessarily adequate
##       - The sample size n is *too* small here for properly capturing GARCH effects.
##         Again, this is only for demonstration purposes here.
uspec <- c(rep(list(ugarchspec(distribution.model = "std")), d-2), # ARMA(1,1)-GARCH(1,1)
           list(ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,2)),
                           distribution.model = "std")),
           list(ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,1)),
                           mean.model = list(armaOrder = c(1,2), include.mean = TRUE),
                           distribution.model = "std")))
system.time(fitAG <- fit_ARMA_GARCH(X, ugarchspec.list = uspec))
str(fitAG, max.level = 1) # list with components fit, warning, error
## Now access the list to check

\dontrun{
## Pick out the standardized residuals, plot them and fit a t copula to them
## Note: ugarchsim() needs the residuals to be standardized; working with
##       standardize = FALSE still requires to simulate them from the
##       respective standardized marginal distribution functions.
Z <- sapply(fitAG$fit, residuals, standardize = TRUE)
U <- pobs(Z)
pairs(U, gap = 0)
system.time(fitC <- fitCopula(tCopula(dim = d, dispstr = "un"), data = U,
                              method = "mpl"))

## Simulate (standardized) Z
set.seed(271)
U. <- rCopula(n, fitC@copula) # simulate from the fitted copula
nu <- sapply(1:d, function(j) fitAG$fit[[j]]@fit$coef["shape"]) # extract (fitted) d.o.f. nu
Z. <- sapply(1:d, function(j) sqrt((nu[j]-2)/nu[j]) * qt(U.[,j], df = nu[j])) # Z

## Simulate from fitted model
X. <- sapply(1:d, function(j)
    fitted(ugarchsim(fitAG$fit[[j]], n.sim = n, m.sim = 1, startMethod = "sample",
                     rseed = 271, custom.dist = list(name = "sample",
                                                     distfit = Z.[,j, drop = FALSE]))))

## Plots original vs simulated -log-returns
opar <- par(no.readonly = TRUE)
layout(matrix(1:(2*d), ncol = d)) # layout
ran <- range(X, X.)
for(j in 1:d) {
    plot(X[,j],  type = "l", ylim = ran, ylab = paste(stocks[j], "-log-returns"))
    plot(X.[,j], type = "l", ylim = ran, ylab = "Simulated -log-returns")
}
    par(opar)
}
}
\keyword{ts}