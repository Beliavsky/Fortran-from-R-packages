\name{qqplot_maha}
\alias{qqplot_maha}
\title{QQ Plot for Multivariate Normal Variance Mixtures}
\description{
  Visual goodness-of-fit test for multivariate normal variance mixtures:
  Plotting squared Mahalanobis distances against their theoretical quantiles.
}
\usage{
qqplot_maha(x, qmix, loc, scale, fitnvmix_object,
            trafo.to.normal = FALSE, test = c("KS.AD", "KS", "AD", "none"),
            boot.pars = list(B = 500, level = 0.95),
            plot = TRUE, verbose = TRUE, control = list(),
            digits = max(3, getOption("digits") - 4), plot.pars = list(), ...)
}
\arguments{
  \item{x}{\eqn{(n, d)}-data \code{\link{matrix}}.}
  \item{qmix}{see \code{\link{pnvmix}()}.}
  \item{loc}{see \code{\link{pnvmix}()}.}
  \item{scale}{see \code{\link{pnvmix}()}.}
  \item{fitnvmix_object}{Optional. Object of class \code{"fitnvmix"} typically
  returned by \code{\link{fitnvmix}()}; if provided, \code{x}, \code{qmix},
  \code{loc} and \code{scale} are ignored.}
  \item{trafo.to.normal}{\code{\link{logical}}. If \code{TRUE}, the
  underlying Mahalanobis distances are mapped to normals by a probability-
  quantile-transform so that the resulting QQ plot is essentially a normal
  QQ plot. Defaults to \code{FALSE}.}
  \item{test}{\code{\link{character}} specifying if (and which) GoF test
  shall be performed. \code{"KS"} performs a Kolmogorov-Smirnoff (see
  \code{\link{ks.test}()}), \code{"AD"} an Anderson-Darling test (see
  \code{\link[ADGofTest]{ad.test}()} from the package \code{ADGofTest}
  and \code{"none"} performs no test. By default, \code{test = "KS.AD"}
  in which case both tests are performed.}
  \item{boot.pars}{\code{\link{list}} with elements \code{B}
  (Bootstrap sample size for computing CIs; if \code{B <= 1},
  no Bootstrap is performed) and
  \code{level} specifying the confidence level.}
  \item{plot}{\code{\link{logical}} specifying if the results should be
  plotted.}
  \item{verbose}{see \code{\link{pnvmix}()}.}
  \item{control}{see \code{\link{get_set_param}()}.}
  \item{digits}{integer specifying the number of digits of the test
  statistic and the p-value to be displayed.}
  \item{plot.pars}{\code{\link{list}} specifying plotting parameters such
  as logarithmic axes; see \code{\link{get_set_qqplot_param}()}.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{qqplot_maha()} (invisibly) returns an object of the class
  \code{"qqplot_maha"} for which the methods \code{plot()} and \code{print()}
  are defined. The return object contains, among others, the components
  \describe{
    \item{\code{maha2}}{Sorted, squared Mahalanobis distances of the data
    from \code{loc} wrt to \code{scale}.}
    \item{\code{theo_quant}}{The theoretical quantile function
    evaluated at \code{ppoints(length(maha2))}.}
    \item{\code{boot_CI}}{\eqn{(2,length(maha2))} matrix containing the
    Bootstrap CIs for the empirical quantiles.}
    \item{\code{asymptSE}}{\code{\link{vector}} of length \code{length(maha2)}
    with estimated, asympotic standard errors for the empirical quantiles.}
   }
}
\details{
  If \eqn{X} follows a multivariate normal variance mixture, the distribution of
  the Mahalanobis distance \eqn{D^2 = (X-\mu)^T \Sigma^{-1} (X-\mu)}
  is a gamma mixture whose distribution function can be approximated.

  The function \code{qqplot_maha()} first estimates the theoretical quantiles
  by calling \code{\link{qgammamix}()} and then plots those against
  the empirical squared Mahalanobis distances
  from the data in \code{x} (with \eqn{\mu=}\code{loc} and
  \eqn{\Sigma=}\code{scale}). Furthermore, the function computes
  asymptotic standard errors of the sample quantiles by using an asymptotic
  normality result for the order statistics which
  are used to plot the asymptotic CI; see Fox (2008, p. 35 -- 36).+

  If \code{boot.pars$B > 1} (which is the default), the function additionally
  performs Bootstrap to construct a CI. Note that by default, the plot contains
  both the asymptotic and the Bootstrap CI.

  Finally, depending on the parameter \code{test}, the function performs a
  univariate GoF test of the observed Mahalanobis distances as described above.
  The test result (i.e., the value of the statistic along with a p-value) are
  typically plotted on the second y-axis.

  The return object of class \code{"qqplot_maha"} contains all computed values
  (such as p-value, test-statistics, Bootstrap CIs and more). We highlight that
  storing this return object makes the QQ plot
  quickly reproducible, as in this case, the theoretical quantiles do not need
  to be recomputed.

  For changing plotting parameters (such as logarithmic axes or colors)
  via the argument \code{plot.pars}, see \code{\link{get_set_qqplot_param}()}.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{fitnvmix}()}, \code{\link{get_set_qqplot_param}()},
  \code{\link{rnvmix}()}, \code{\link{pnvmix}()}, \code{\link{dnvmix}()}
}
\examples{
## Sample from a heavy tailed multivariate t and construct QQ plot
set.seed(1)
d <- 2
n <- 1000
df <- 3.1
rho <- 0.5
loc <- rep(0, d)
scale <- matrix(c(1, rho, rho, 1), ncol = 2)
qmix <- "inverse.gamma"
## Sample data
x <- rnvmix(n, qmix = qmix, loc = loc, scale = scale, df = df)
# Construct QQ Plot with 'true' parameters and store result object
qq1 <- qqplot_maha(x, qmix = qmix, df = df, loc = loc, scale = scale)
## ... which is an object of class "qqplot_maha" with two methods
stopifnot(class(qq1) == "qqplot_maha", "plot.qqplot_maha" \%in\% methods(plot),
          "print.qqplot_maha" \%in\% methods(print))
plot(qq1) # reproduce the plot
plot(qq1, plotpars = list(log = "xy")) # we can also plot on log-log scale

## In fact, with the 'plotpars' argument, we can change a lot of things
plot(qq1, plotpars = list(col = rep("black", 4), lty = 4:6, pch = "*",
                          plot_test = FALSE, main = "Same with smaller y limits",
                          sub = "MySub", xlab = "MyXlab", ylim = c(0, 1.5e3)))

## What about estimated parameters?
myfit <- fitStudent(x)
## We can conveniently pass 'myfit', rather than specifying 'x', 'loc', ...
set.seed(1)
qq2.1 <- qqplot_maha(fitnvmix_object = myfit, test = "AD", trafo_to_normal = TRUE)
set.seed(1)
qq2.2 <- qqplot_maha(x, qmix = "inverse.gamma", loc = myfit$loc,
                     scale = myfit$scale, df = myfit$df,
                     test = "AD", trafo_to_normal = TRUE)
stopifnot(all.equal(qq2.1$boot_CI, qq2.2$boot_CI)) # check
qq2.2 # it mentions here that the Maha distances were transformed to normal


## Another example where 'qmix' is a function, so quantiles are internally
## estimated via 'qgammamix()'
n <- 15 # small sample size to have examples run fast
## Define the quantile function of an IG(nu/2, nu/2) distribution
qmix <- function(u, df) 1 / qgamma(1 - u, shape = df/2, rate = df/2)
## Sample data
x <- rnvmix(n, qmix = qmix, df = df, loc = loc, scale = scale)
## QQ Plot of empirical quantiles vs true quantiles, all values estimated
## via RQMC:
set.seed(1)
qq3.1 <- qqplot_maha(x, qmix = qmix, loc = loc, scale = scale, df = df)
## Same could be obtained by specifying 'qmix' as string in which case
## qqplot_maha() calls qf()
set.seed(1)
qq3.2 <- qqplot_maha(x, qmix = "inverse.gamma", loc = loc, scale = scale, df = df)
}
\keyword{distribution}