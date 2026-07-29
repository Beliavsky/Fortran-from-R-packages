\name{dependencemeasures}
\alias{lambda_gStudent}
\alias{corgnvmix}
\title{Dependence Measures for grouped normal variance mixture copulas}
\description{
  Computation of rank correlation coefficients Spearman's rho and Kendall's tau
  for grouped normal variance mixture copulas as well as computation of the
  (lower and upper) tail dependence coefficient of a grouped \emph{t} copula.
}
\usage{
corgnvmix(scale, qmix, method = c("kendall", "spearman"), groupings = 1:2,
          ellip.kendall = FALSE, control = list(), verbose = TRUE, ...)

lambda_gStudent(df, scale, control = list(), verbose = TRUE)
}
\arguments{
  \item{scale}{\eqn{n}-\code{\link{vector}} giving the \eqn{\rho} parameters of
  the copula.}
  \item{qmix}{specification of the mixing variables; see \code{\link{pgnvmix}()}.}
  \item{method}{\code{\link{character}} indicating if Spearman's rho or
  Kendall's tau is to be computed.}
  \item{groupings}{\code{\link{vector}} specifying the grouping structure; either
  \code{rep(1, 2)} (ungrouped) or \code{1:2} (grouped case).}
  \item{ellip.kendall}{\code{\link{logical}} if the formula for Kendalll's tau
  for elliptical copulas shall be used; see details below.}
  \item{df}{either scalar or \eqn{2}-\code{\link{vector}} giving the degrees-of-
  freedoms for the t copula; if provided as scalar, the copula is an (ungrouped)
  t copula and \code{lambda_gStudent()} uses a closed formula.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    has not been reached.}
    \item{\dots}{additional arguments (for example, parameters) passed
    to the underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{lambda_gStudent()} and \code{corgnvmix()} return
  a \code{\link{numeric}} \eqn{n}-vector with the computed
  dependence measure with corresponding attributes
  \code{"abs. error"} and \code{"rel. error"}(error estimates of the RQMC estimator)
  and \code{"numiter"} (number of iterations).
}
\details{
  For grouped normal variance mixture copulas, including the grouped \emph{t},
  there is no closed formula for Kendall's tau and Spearman's rho. The function
  \code{corgnvmix()} approximates these dependence measures by numerically
  approximating an integral representation for these measures.

  If no grouping is present (i.e., when \code{groupings = rep(1, 2)}), the
  copula is an elliptical copula for which the formula \eqn{\tau = 2asin(\rho)/pi}
  holds. This formula holds only approximately in the grouped case; the quality
  of the approximation depends on how different the mixing variables for the
  two components are. When the mixing distributions are not too far apart and
  when the copula parameter is not close to 1, this approximation is ``very
  accurate``, as demonstrated in Daul et al (2003).

  In the ungrouped case, \code{lambda_gStudent()} computes the tail dependence
  coefficient \eqn{lambda} based on the known formula
  \code{2 * pt( -sqrt( (df + 1)*(1 - rho) / (1 + rho)), df = df + 1)} for
  the tail dependence coefficient of a \emph{t} copula.

  In the grouped case, RQMC methods are used to efficiently approximate the
  integral given in Eq. (26) of Luo and Shevchenko (2010).

}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2020),
  Grouped Normal Variance Mixtures.
  \emph{Risks} 8(4), 103.

  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.

  Luo, X. and Shevchenko, P. (2010).
  The \emph{t} copula with multiple parameters of degrees of freedom: bivariate
  characteristics and application to risk management.
  \emph{Quantitative Finance} 10(9),
  1039-1054.

  Daul, S., De Giorgi, E. G., Lindskog, F. and McNeil, A (2003).
  The grouped \emph{t} copula with an application to credit risk.
  \emph{Available at SSRN 1358956}.
}
\seealso{
  \code{\link{dgStudentcopula}()}, \code{\link{pgStudentcopula}()},
  \code{\link{rgStudentcopula}()}
}
\examples{
### Examples for corgnvmix() ###################################################

## Create a plot displaying Spearman's rho for a grouped t copula as a function
## of the copula parameter for various choices of the degrees-of-freedom
qmix <- "inverse.gamma"
df <- matrix( c(1, 2, 1, 5, 1, Inf), ncol = 2, byrow = TRUE)
l.df <- nrow(df)
scale <- seq(from = 0, to = 1, length.out = 99)
set.seed(1) # for reproducibility
kendalls <- sapply(seq_len(l.df), function(i)
   corgnvmix(scale, qmix = qmix, method = "kendall", df = df[i, ]))
## Include the elliptical approximation (exact when df1 = df2)
kendall_ell <- corgnvmix(scale, method = "kendall", ellip.kendall = TRUE)
## Plot
lgnd <- character(l.df + 1)
lgnd[1] <- "elliptical (equal df)"
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = expression(rho),
     ylab = "Kendall's tau")
lines(scale, kendall_ell, lty = 1)
for(i in 1:l.df){
   lines(scale, kendalls[, i], col = i + 1, lty = i + 1)
   lgnd[i+1] <- paste0("df1 = ", df[i, 1], ", df2 = ", df[i, 2])
}
legend("topleft", lgnd, col = 1:(l.df + 1), lty = 1:(l.df + 1), bty = 'n')


### Examples for lambda_gStudent() #############################################

## Create a plot displaying 'lambda' as a function of the copula parameter
## for various choices of the degrees-of-freedom
df <- c(3, 6, 9)
df_ <- list( rep(df[1], 2), rep(df[2], 2), rep(df[3], 2), # ungrouped
             c(df[1], df[2]), c(df[1], df[3]), c(df[2], df[3])) # grouped
l.df_ <- length(df_)
scale <- seq(from = -0.99, to = 0.99, length.out = 112) # scale parameters
set.seed(1) # for reproducibilty
lambdas <-
   sapply(seq_len(l.df_), function(i) lambda_gStudent(df_[[i]], scale = scale))
lgnd <- character(length(df_))
plot(NA, xlim = range(scale), ylim = range(lambdas), xlab = expression(rho),
     ylab = expression(lambda))
for(i in seq_len(l.df_)){
   lines(scale, lambdas[, i], col = i, lty = i)
   lgnd[i] <- if(df_[[i]][1] == df_[[i]][2]) paste0("df = ", df_[[i]][1]) else
      paste0("df1 = ", df_[[i]][1], ", df2 = ", df_[[i]][2])
}
legend("topleft", lgnd, col = seq_len(l.df_), lty = seq_len(l.df_),
       bty = 'n')
## If called with 'df' a 1-vector, closed formula for lambda is used => check
lambda.true <- sapply(1:3, function(i) lambda_gStudent(df_[[i]][1], scale = scale))
stopifnot(max(abs( lambda.true - lambdas[, 1:3])) < 4e-4)
}
\keyword{distribution}