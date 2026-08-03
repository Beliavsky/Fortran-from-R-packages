\name{CLA}
\alias{CLA}
\title{Critical Line Algorithm for mean-variance optimal portfolio}
\description{
  The Critical Line Algorithm was first proposed by Markowitz(1987) to
  solve the mean-variance optimal portfolio problem.

  We solve the problem with \dQuote{box} constraints, i.e., allow to
  specify lower and upper bounds (via \code{lB} and \code{uB}) for each
  asset weight.

  Here we provide a pure \R{} implementation, quite fine tuned and
  debugged compared to earlier ones.
}
\usage{
CLA(mu, covar, lB, uB,
    check.cov = TRUE, check.f = TRUE,
    tol.lambda = 1e-07,
    give.MS = TRUE, keep.names = TRUE, trace = 0)
}
\arguments{
  \item{mu}{numeric vector of length \code{n} containing the expected
    return \eqn{E[R_i]} for \eqn{1=1,2,\dots,n}.}
  \item{covar}{the \eqn{n \times n}{n x n} covariance matrix of the
    returns, must be positive definite.}

  \item{lB, uB}{vectors of length \code{n} with lower and upper bounds
    for the asset weights.}

  \item{check.cov}{\code{\link{logical}} indicating if the \code{covar}
    matrix should be checked to be positive definite.}
  \item{check.f}{\code{\link{logical}} indicating if a warning should be
    produced when the algorithm cannot produce a new (smaller) lambda even
    though there are still \bold{f}ree weights to be chosen.}

  \item{tol.lambda}{the tolerance when checking for lambda changes or
    being zero.}
  \item{give.MS}{\code{\link{logical}} indicating if \code{\link{MS}()}
    should be computed (and returned) as well.}
  \item{keep.names}{\code{\link{logical}} indicating if the
    \code{weights_set} matrix should keep the (asset) \code{names(mu)}.}
  \item{trace}{an integer (or \code{\link{logical}}) indicating if and
    how much diagnostic or progress output should be produced.}
}
\details{
  The current implementation of the CLA is based (via Norring's)
  on Bailey et al.(2013).   We have found buglets in that implementation
  which lead them to introduce their \dQuote{purge} routines
  (\code{purgeNumErr}, \code{purgeExcess}),% in Y. Shi's masters' thesis
  which are no longer necessary.

  Even though this is a pure \R implementation, the algorithm is quite fast
  also when the number of assets \eqn{n} is large (1000s), though that
  depends quite a bit on the exact problem.
}
\note{
  The exact results of the algorithm, e.g., the assets with non-zero
  weights, may slightly depend on the (computer) platform, e.g., for the
  \acronym{S&P 500} example, differences between 64-bit or 32-bit, version
  of BLAS or Lapack libraries etc, do have an influence, see the \R script
  \file{tests/SP500-ex.R} %  ../tests/SP500-ex.R}
  in the package sources.
}
\value{
  an object of \code{\link{class}} \code{"CLA"} which is a
  \code{\link{list}} with components

  \item{weights_set}{a \eqn{n \times m}{n * m} matrix of asset weights,
    corresponding to the \eqn{m} steps that the CLA has completed or the
    \eqn{m} \dQuote{turning points} it has computed.}
  \item{free_indices}{a \code{\link{list}} of length \code{m}, the
    \eqn{k}-th component with the indices in \eqn{{1,\dots,n}} of
    those assets whose weights were not at the boundary after ... }% <<< FIXME
  \item{gammas}{numeric vector of length \eqn{m} of the values
    \eqn{\gamma_k}{gamma[k]} for CLA step \eqn{k}, \eqn{k=1,\dots,n}.}
  \item{lambdas}{numeric vector of length \eqn{m} of the Lagrange parameters
    \eqn{\lambda_k}{lambda[k]} for CLA step \eqn{k}, \eqn{k=1,\dots,n}.}
  \item{MS_weights}{the \eqn{\mu(W)} and \eqn{\sigma(W)} corresponding
    to the asset weights \code{weights_set}, i.e., simply the same as
    \code{\link{MS}(weights_set = weights_set, mu = mu, covar = covar)}.}
}
\references{
  Markowitz, H. (1952)
  Portfolio selection, \emph{The Journal of Finance} \bold{7}, 77--91;
  \doi{10.2307/2975974}.

  Markowitz, H. M. (1987, 1st ed.) and
  Markowitz, H. M. and Todd, P. G. (2000)
  \emph{Mean-Variance Analysis in Portfolio Choice and Capital Markets};
  chapters 7 and 13.

  %% Markowitz, H., Todd, P., Xu, G. and Yamane, Y. (1993)
  %% Computation of mean-semivariance efficient sets by the Critical Line Algorithm,
  %% \emph{Annals of Operations Research} \bold{45}(1): 307--317;
  %% \doi{10.1007/BF02282055}.

  Niedermayer, A. and Niedermayer, D. (2010)
  Applying Markowitz’s Critical Line Algorithm, in J. B. Guerard (ed.),
  Handbook of Portfolio Construction, Springer; chapter 12, 383--400;
  \doi{10.1007/978-0-387-77439-8_12}.

  Bailey, D. H. and López de Prado, M. (2013)
  An open-source implementation of the critical-line algorithm for portfolio
  optimization, \emph{Algorithms} \bold{6}(1), 169--196;
  \doi{10.3390/a6010169},

  Yanhao Shi (2017)
  Implementation and applications of critical line algorithm for
  portfolio optimization; unpublished Master's thesis, ETH Zurich.
  % ~/Betreute-Arbeiten/YanhaoShi/MasterThesis-YanhaoShi.pdf
}
%% \note{
%% }
\author{Alexander Norring did the very first version (unpublished master
  thesis). Current implementation: Yanhao Shi and Martin Maechler
}
\seealso{
  \code{\link{MS}};
  for plotting \code{CLA} results: \code{\link{plot.CLA}}.
}
\examples{
data(muS.sp500)
## Full data taking too much time for example
set.seed(47)
iS <- sample.int(length(muS.sp500$mu), 24)

CLsp.24 <- CLA(muS.sp500$mu[iS], muS.sp500$covar[iS, iS], lB=0, uB=1/10)
CLsp.24 # using the print() method for class "CLA"

plot(CLsp.24)

if(require(Matrix)) { ## visualize how weights change "along turning points"
  show(image(Matrix(CLsp.24$weights_set, sparse=TRUE),
             main = "CLA(muS.sp500 <random_sample(size=24)>) $ weights_set",
             xlab = "turning point", ylab = "asset number"))
}

## A 3x3 example (from real data) where CLA()'s original version failed
## and  'check.f = TRUE' produces a warning :
mc3 <- list(
    mu = c(0.0408, 0.102, -0.023),
    cv = matrix(c(0.00648, 0.00792, 0.00473,
                  0.00792, 0.0334,  0.0121,
                  0.00473, 0.0121, 0.0793), 3, 3,
           dimnames = list(NULL,
                           paste0(c("TLT", "VTI","GLD"), ".Adjusted"))))

rc3 <- with(mc3,  CLA(mu=mu, covar=cv, lB=0, uB=1, trace=TRUE))

}
\keyword{optimize}
\keyword{arith}

