\name{alloc}
\alias{alloc_ellip}
\alias{conditioning}
\alias{alloc_np}
\title{Computing allocations}
\description{
  Computing (capital) allocations.
}
\usage{
## For elliptical distributions under certain assumptions
alloc_ellip(total, loc, scale)

## Nonparametrically
conditioning(x, level, risk.measure = "VaR_np", ...)
alloc_np(x, level, risk.measure = "VaR_np", include.conditional = FALSE, ...)
}
\arguments{
  \item{total}{total to be allocated (typically the risk measure of the
    sum of the underlying loss random variables).}
  \item{loc}{location vector of the elliptical distribution of the loss
    random vector.}
  \item{scale}{scale (covariance) matrix of the elliptical distribution of
    the loss random vector.}
  \item{x}{\eqn{(n, d)}-matrix containing \eqn{n} iid
    \eqn{d}-dimensional losses.}
  \item{level}{either one or two confidence level(s) for
    \code{risk.measure}; in the former case the upper bound on the
    conditioning region is determined by confidence level 1.}
  \item{risk.measure}{\code{\link{character}} string or
    \code{\link{function}} specifying the risk measure to be computed
    on the row sums of \code{x} based on the given level(s) in order
    to determine the conditioning region.}
  \item{include.conditional}{\code{\link{logical}} indicating whether
    the computed sub-sample of \code{x} is to be returned, too.}
  \item{\dots}{additional arguments passed to \code{risk.measure}.}
}
\value{
  \eqn{d}-vector of allocated amounts (the allocation) according to the
  Euler principle under the assumption that the underlying loss random
  vector follows a \eqn{d}-dimensional elliptical distribution with
  location vector \code{loc} (\eqn{\bm{mu}}{mu} in the reference) and
  scale matrix \code{scale} (\eqn{\Sigma}{Sigma} in the reference, a
  covariance matrix) and that the risk measure is law-invariant,
  positive-homogeneous and translation invariant.
}
\details{
  The result of \code{alloc_ellip()} for \code{loc = 0} can be found in
  McNeil et al. (2015, Corollary 8.43). Otherwise, McNeil et al. (2015,
  Theorem 8.28 (1)) can be used to derive the result.
}
\author{Marius Hofert}
\references{
  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
### Ellipitical case ###########################################################

## Construct a covariance matrix
sig <- 1:3 # standard deviations
library(copula) # for p2P() here
P <- p2P(c(-0.5, 0.3, 0.5)) # (3, 3) correlation matrix
Sigma <- P * sig \%*\% t(sig) # corresponding covariance matrix
stopifnot(all.equal(cov2cor(Sigma), P)) # sanity check

## Compute the allocation of 1.2 for a joint loss L ~ E_3(0, Sigma, psi)
AC <- alloc_ellip(1.2, loc = 0, scale = Sigma) # allocated amounts
stopifnot(all.equal(sum(AC), 1.2)) # sanity check
## Be careful to check whether the aforementioned assumptions hold.


### Nonparametrically ##########################################################

## Generate data
set.seed(271)
X <- qt(rCopula(1e5, copula = gumbelCopula(2, dim = 5)), df = 3.5)

## Estimate an allocation via MC based on a sub-sample whose row sums have a
## nonparametric VaR with confidence level in ...
alloc_np(X, level = 0.9) # ... (0.9, 1]
CA  <- alloc_np(X, level = c(0.9, 0.95)) # ... in (0.9, 0.95]
CA. <- alloc_np(X, level = c(0.9, 0.95), risk.measure = VaR_np) # providing a function
stopifnot(identical(CA, CA.))
}
\keyword{models}