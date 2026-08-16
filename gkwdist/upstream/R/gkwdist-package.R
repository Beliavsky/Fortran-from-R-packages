#' Generalized Kumaraswamy Distribution Family
#'
#' @description
#' The \pkg{gkwdist} package provides a comprehensive, high-performance
#' implementation of the five-parameter Generalized Kumaraswamy (GKw)
#' distribution and its seven nested sub-families for modeling bounded
#' continuous data on the unit interval \eqn{(0,1)}. All distribution and
#' likelihood functions are implemented in \strong{C++} via
#' \pkg{RcppArmadillo} for maximum computational efficiency.
#'
#' @details
#' \subsection{Overview}{
#' The Generalized Kumaraswamy distribution, proposed by Carrasco, Ferrari and
#' Cordeiro (2010), extends the Kumaraswamy distribution (Jones, 2009) by
#' incorporating three additional shape parameters. This generalization provides
#' remarkable flexibility in accommodating various density shapes including
#' unimodality, bimodality, asymmetry, J-shapes, inverted J-shapes, U-shapes,
#' bathtub shapes, and heavy or light tails.
#' }
#'
#' \subsection{Mathematical Specification}{
#' The probability density function (PDF) of the GKw distribution with
#' parameters \eqn{\theta = (\alpha, \beta, \gamma, \delta, \lambda)}
#' (all positive) is:
#' \deqn{f(x; \theta) = \frac{\lambda\alpha\beta x^{\alpha-1}}{B(\gamma, \delta + 1)}
#' (1 - x^\alpha)^{\beta-1} [1 - (1 - x^\alpha)^\beta]^{\gamma\lambda-1}
#' \{1 - [1 - (1 - x^\alpha)^\beta]^\lambda\}^\delta}
#' for \eqn{0 < x < 1}, where \eqn{B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)}
#' is the beta function.
#'
#' The cumulative distribution function (CDF) is:
#' \deqn{F(x; \theta) = I_{[1-(1-x^\alpha)^\beta]^\lambda}(\gamma, \delta + 1)}
#' where \eqn{I_y(a,b)} denotes the regularized incomplete beta function ratio.
#' }
#'
#' \subsection{Distribution Family Hierarchy}{
#' The GKw distribution nests seven important sub-models:
#' \describe{
#'   \item{\strong{Beta-Kumaraswamy (BKw)}}{4 parameters:
#'     \eqn{(\alpha, \beta, \gamma, \delta)}. Set \eqn{\lambda = 1} in GKw.}
#'   \item{\strong{Kumaraswamy-Kumaraswamy (KKw)}}{4 parameters:
#'     \eqn{(\alpha, \beta, \delta, \lambda)}. Set \eqn{\gamma = 1} in GKw.}
#'   \item{\strong{Exponentiated Kumaraswamy (EKw)}}{3 parameters:
#'     \eqn{(\alpha, \beta, \lambda)}. Set \eqn{\gamma = 1, \delta = 0} in GKw.
#'     Has closed-form quantile function.}
#'   \item{\strong{McDonald (Mc)}}{3 parameters:
#'     \eqn{(\gamma, \delta, \lambda)}. Set \eqn{\alpha = \beta = 1} in GKw.
#'     Also known as Beta-Power distribution.}
#'   \item{\strong{Kumaraswamy (Kw)}}{2 parameters: \eqn{(\alpha, \beta)}.
#'     Set \eqn{\gamma = \delta = \lambda = 1} in GKw. Has closed-form CDF
#'     and quantile function.}
#'   \item{\strong{Beta}}{2 parameters: \eqn{(\gamma, \delta)}.
#'     Set \eqn{\alpha = \beta = \lambda = 1} in GKw. Classical beta
#'     distribution with \code{shape1 = gamma}, \code{shape2 = delta + 1}.}
#'   \item{\strong{Uniform}}{0 parameters. Set all shape parameters to 1.}
#' }
#' }
#'
#' @section Distribution Functions:
#' The package provides the standard R distribution function API for all
#' seven distributions. Each distribution has four core functions with
#' prefix \code{d}, \code{p}, \code{q}, \code{r}:
#'
#' \strong{Generalized Kumaraswamy (GKw):}
#' \itemize{
#'   \item \code{\link{dgkw}}: Density function
#'   \item \code{\link{pgkw}}: Distribution function (CDF)
#'   \item \code{\link{qgkw}}: Quantile function (inverse CDF)
#'   \item \code{\link{rgkw}}: Random generation
#' }
#'
#' \strong{Beta-Kumaraswamy (BKw):}
#' \code{\link{dbkw}}, \code{\link{pbkw}}, \code{\link{qbkw}}, \code{\link{rbkw}}
#'
#' \strong{Kumaraswamy-Kumaraswamy (KKw):}
#' \code{\link{dkkw}}, \code{\link{pkkw}}, \code{\link{qkkw}}, \code{\link{rkkw}}
#'
#' \strong{Exponentiated Kumaraswamy (EKw):}
#' \code{\link{dekw}}, \code{\link{pekw}}, \code{\link{qekw}}, \code{\link{rekw}}
#'
#' \strong{McDonald (Mc):}
#' \code{\link{dmc}}, \code{\link{pmc}}, \code{\link{qmc}}, \code{\link{rmc}}
#'
#' \strong{Kumaraswamy (Kw):}
#' \code{\link{dkw}}, \code{\link{pkw}}, \code{\link{qkw}}, \code{\link{rkw}}
#'
#' \strong{Beta:}
#' \code{\link{dbeta_}}, \code{\link{pbeta_}}, \code{\link{qbeta_}}, \code{\link{rbeta_}}
#'
#' All distribution functions are implemented in C++ for optimal performance.
#'
#' @section Likelihood Functions:
#' High-performance analytical log-likelihood, gradient, and Hessian functions
#' are provided for maximum likelihood estimation. These functions return
#' \strong{negative} values to facilitate direct use with optimization routines
#' like \code{\link[stats]{optim}}, which perform minimization by default.
#'
#' Function signature: \code{function(par, data)} where \code{par} is a numeric
#' vector of parameters and \code{data} is the observed sample.
#'
#' \describe{
#'   \item{\code{ll*}}{Negative log-likelihood:
#'     \eqn{-\ell(\theta) = -\sum_{i=1}^n \log f(x_i; \theta)}}
#'   \item{\code{gr*}}{Negative gradient (negative score vector):
#'     \eqn{-\nabla_\theta \ell(\theta)}}
#'   \item{\code{hs*}}{Negative Hessian matrix:
#'     \eqn{-\nabla^2_\theta \ell(\theta)}}
#' }
#'
#' Available for all distributions:
#' \itemize{
#'   \item \strong{GKw:} \code{\link{llgkw}}, \code{\link{grgkw}}, \code{\link{hsgkw}}
#'   \item \strong{BKw:} \code{\link{llbkw}}, \code{\link{grbkw}}, \code{\link{hsbkw}}
#'   \item \strong{KKw:} \code{\link{llkkw}}, \code{\link{grkkw}}, \code{\link{hskkw}}
#'   \item \strong{EKw:} \code{\link{llekw}}, \code{\link{grekw}}, \code{\link{hsekw}}
#'   \item \strong{Mc:} \code{\link{llmc}}, \code{\link{grmc}}, \code{\link{hsmc}}
#'   \item \strong{Kw:} \code{\link{llkw}}, \code{\link{grkw}}, \code{\link{hskw}}
#'   \item \strong{Beta:} \code{\link{llbeta}}, \code{\link{grbeta}}, \code{\link{hsbeta}}
#' }
#'
#' These analytical derivatives are exact (up to floating-point precision) and
#' provide substantial performance gains over numerical approximations.
#'
#' @section Utility Functions:
#' \itemize{
#'   \item \code{\link{gkwgetstartvalues}}: Compute starting values for
#'     optimization using method of moments or quantile matching
#' }
#'
#' @section Applications:
#' The GKw distribution family is particularly suitable for modeling:
#' \itemize{
#'   \item \strong{Proportions and rates}: Bounded continuous data in \eqn{(0,1)}
#'   \item \strong{Percentages}: Economic indices, market shares, conversion rates
#'   \item \strong{Fractions}: Parts of a whole, concentration measures
#'   \item \strong{Scores and indices}: Normalized measurements, standardized tests
#'   \item \strong{Reliability data}: Component lifetimes on bounded domains
#'   \item \strong{Hydrological phenomena}: Reservoir levels, soil moisture content
#'   \item \strong{Financial ratios}: Debt-to-equity, current ratio, profit margins
#'   \item \strong{Quality metrics}: Defect rates, efficiency scores, purity levels
#'   \item \strong{Biostatistics}: Survival probabilities, dose-response curves
#' }
#'
#' @section Advantages Over Standard Distributions:
#' Compared to the classical Beta distribution, the GKw family offers:
#' \itemize{
#'   \item Greater flexibility in density shapes (bimodal, U-shaped, bathtub)
#'   \item Better accommodation of asymmetry and tail behavior
#'   \item Closed-form CDF and quantile for some sub-models (Kw, EKw, KKw)
#'   \item Computational efficiency via C++ implementation
#'   \item Easy parameter interpretation through nested structure
#'   \item Superior performance for extreme parameter values
#' }
#'
#' @section Performance:
#' All functions are implemented in C++ using \pkg{RcppArmadillo}, providing:
#' \itemize{
#'   \item 10-100× speedup over pure R implementations
#'   \item Linear scaling with sample size
#'   \item Optimized memory allocation
#'   \item Numerical stability for extreme parameter values
#'   \item Efficient vectorized operations
#' }
#'
#' Typical benchmarks on modern hardware:
#' \itemize{
#'   \item Density evaluation: \eqn{>10^7} evaluations per second
#'   \item Log-likelihood: \eqn{n = 10^6} observations in \eqn{<100}ms
#'   \item Gradient computation: \eqn{<5}× slower than log-likelihood
#' }
#'
#' @section Model Selection Workflow:
#' Recommended strategy for choosing among distributions:
#' \enumerate{
#'   \item \strong{Exploratory Analysis}: Examine histograms, kernel density
#'     estimates, and summary statistics of your data.
#'   \item \strong{Start Simple}: Fit Beta and Kumaraswamy distributions
#'     (2 parameters). Use \code{\link[stats]{optim}} with \code{method = "BFGS"}
#'     and analytical gradients.
#'   \item \strong{Diagnostic Checking}: Assess fit quality using Q-Q plots,
#'     probability plots, and goodness-of-fit tests (e.g., Kolmogorov-Smirnov).
#'   \item \strong{Progressive Complexity}: If simple models inadequate, try
#'     3-parameter models (EKw or Mc), then 4-parameter models (BKw or KKw).
#'   \item \strong{Information Criteria}: Use AIC, BIC, or AICc to balance
#'     goodness-of-fit and model parsimony. Lower values indicate better models.
#'   \item \strong{Nested Testing}: Perform likelihood ratio tests when
#'     comparing nested models (e.g., Kw vs. EKw).
#'   \item \strong{Cross-Validation}: Validate final model on held-out data
#'     or via bootstrap procedures.
#'   \item \strong{Residual Analysis}: Examine probability integral transform
#'     residuals for uniformity and independence.
#' }
#'
#' @section Statistical Inference:
#' Maximum likelihood estimation is performed using numerical optimization:
#' \enumerate{
#'   \item Obtain starting values via \code{\link{gkwgetstartvalues}} or
#'     manual specification based on sample moments
#'   \item Minimize negative log-likelihood using \code{\link[stats]{optim}}
#'     with \code{method = "BFGS"} or \code{"L-BFGS-B"}
#'   \item Provide analytical gradient via \code{gr} argument for faster
#'     convergence and improved accuracy
#'   \item Set \code{hessian = TRUE} to obtain observed information matrix
#'   \item Compute standard errors as \code{sqrt(diag(solve(hessian)))}
#'   \item Construct confidence intervals using normal approximation or
#'     profile likelihood
#'   \item Perform hypothesis tests using Wald, score, or likelihood ratio
#'     statistics
#' }
#'
#' For large samples (\eqn{n > 10^4}), consider using \code{method = "L-BFGS-B"}
#' with box constraints to prevent parameter estimates from drifting to
#' boundary values or becoming numerically unstable.
#'
#' @author
#' J. E. Lopes \email{evandeilton@@gmail.com}
#' \href{https://orcid.org/0009-0007-5887-4084}{ORCID: 0009-0007-5887-4084}
#'
#' LEG - Laboratory of Statistics and Geoinformation
#'
#' PPGMNE - Graduate Program in Numerical Methods in Engineering
#'
#' Federal University of Paraná (UFPR), Brazil
#'
#' @references
#' Carrasco, J. M. F., Ferrari, S. L. P., and Cordeiro, G. M. (2010).
#' A new generalized Kumaraswamy distribution.
#' \emph{arXiv preprint arXiv:1004.0911}.
#' \doi{10.48550/arXiv.1004.0911}
#'
#' Jones, M. C. (2009).
#' Kumaraswamy's distribution: A beta-type distribution with some tractability
#' advantages.
#' \emph{Statistical Methodology}, \strong{6}(1), 70-81.
#' \doi{10.1016/j.stamet.2008.04.001}
#'
#' Kumaraswamy, P. (1980).
#' A generalized probability density function for double-bounded random processes.
#' \emph{Journal of Hydrology}, \strong{46}(1-2), 79-88.
#' \doi{10.1016/0022-1694(80)90036-0}
#'
#' Cordeiro, G. M., and de Castro, M. (2011).
#' A new family of generalized distributions.
#' \emph{Journal of Statistical Computation and Simulation}, \strong{81}(7), 883-898.
#' \doi{10.1080/00949650903530745}
#'
#' McDonald, J. B. (1984).
#' Some generalized functions for the size distribution of income.
#' \emph{Econometrica}, \strong{52}(3), 647-663.
#' \doi{10.2307/1913469}
#'
#' Eddelbuettel, D., and Francois, R. (2011).
#' Rcpp: Seamless R and C++ Integration.
#' \emph{Journal of Statistical Software}, \strong{40}(8), 1-18.
#' \doi{10.18637/jss.v040.i08}
#'
#' Eddelbuettel, D., and Sanderson, C. (2014).
#' RcppArmadillo: Accelerating R with high-performance C++ linear algebra.
#' \emph{Computational Statistics & Data Analysis}, \strong{71}, 1054-1063.
#' \doi{10.1016/j.csda.2013.02.005}
#'
#' @seealso
#' \code{\link[stats]{Beta}} for the standard beta distribution,
#' \code{\link[stats]{optim}} for numerical optimization,
#' \code{\link[stats]{dbeta}} for beta distribution functions
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Distribution Functions
#'
#' library(gkwdist)
#'
#' # Set parameters for GKw distribution
#' alpha <- 2.0
#' beta <- 3.0
#' gamma <- 1.5
#' delta <- 2.0
#' lambda <- 1.2
#'
#' # Create sequence of x values
#' x <- seq(0.01, 0.99, length.out = 200)
#'
#' # Compute density
#' dens <- dgkw(x, alpha, beta, gamma, delta, lambda)
#'
#' # Compute CDF
#' cdf <- pgkw(x, alpha, beta, gamma, delta, lambda)
#'
#' # Compute specific quantiles
#' probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' quants <- qgkw(probs, alpha, beta, gamma, delta, lambda)
#' print(round(quants, 4))
#'
#' # Generate random sample
#' set.seed(123)
#' sample <- rgkw(1000, alpha, beta, gamma, delta, lambda)
#'
#' # PDF
#' plot(x, dens,
#'   type = "l", lwd = 2, col = "darkblue",
#'   main = "GKw Probability Density Function",
#'   xlab = "x", ylab = "f(x)", las = 1
#' )
#' grid(col = "gray80", lty = 2)
#'
#' # CDF
#' plot(x, cdf,
#'   type = "l", lwd = 2, col = "darkred",
#'   main = "GKw Cumulative Distribution Function",
#'   xlab = "x", ylab = "F(x)", las = 1
#' )
#' grid(col = "gray80", lty = 2)
#'
#' # Histogram of random sample
#' hist(sample,
#'   breaks = 30, probability = TRUE,
#'   col = "lightblue", border = "white",
#'   main = "Random Sample from GKw",
#'   xlab = "x", ylab = "Density", las = 1
#' )
#' lines(x, dens, col = "darkblue", lwd = 2)
#' grid(col = "gray80", lty = 2)
#'
#' # Q-Q plot
#' theoretical <- qgkw(ppoints(length(sample)), alpha, beta, gamma, delta, lambda)
#' empirical <- sort(sample)
#' plot(theoretical, empirical,
#'   pch = 19, col = rgb(0, 0, 1, 0.3),
#'   main = "Q-Q Plot",
#'   xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", las = 1
#' )
#' abline(0, 1, col = "red", lwd = 2, lty = 2)
#' grid(col = "gray80", lty = 2)
#'
#' ## Example 2: Comparing Distribution Families
#'
#' # Create comparison plot
#' x <- seq(0.001, 0.999, length.out = 500)
#'
#' # GKw (5 parameters) - most flexible
#' plot(x, dgkw(x, 2, 3, 1.5, 2, 1.2),
#'   type = "l", lwd = 2, col = "black",
#'   main = "GKw Distribution (5 params)",
#'   xlab = "x", ylab = "Density", las = 1, ylim = c(0, 3)
#' )
#' grid(col = "gray80", lty = 2)
#'
#' # BKw (4 parameters)
#' plot(x, dbkw(x, 2, 3, 1.5, 2),
#'   type = "l", lwd = 2, col = "darkgreen",
#'   main = "BKw Distribution (4 params)",
#'   xlab = "x", ylab = "Density", las = 1, ylim = c(0, 3)
#' )
#' grid(col = "gray80", lty = 2)
#'
#' # EKw (3 parameters)
#' plot(x, dekw(x, 2, 3, 1.5),
#'   type = "l", lwd = 2, col = "darkred",
#'   main = "EKw Distribution (3 params)",
#'   xlab = "x", ylab = "Density", las = 1, ylim = c(0, 3)
#' )
#' grid(col = "gray80", lty = 2)
#'
#' # Kw (2 parameters) - base distribution
#' plot(x, dkw(x, 2, 3),
#'   type = "l", lwd = 2, col = "darkblue",
#'   main = "Kw Distribution (2 params)",
#'   xlab = "x", ylab = "Density", las = 1, ylim = c(0, 3)
#' )
#' grid(col = "gray80", lty = 2)
#'
#'
#' ## Example 3: Maximum Likelihood Estimation
#'
#' # Generate data from Kumaraswamy distribution
#' set.seed(2024)
#' n <- 2000
#' true_alpha <- 2.5
#' true_beta <- 3.5
#' data <- rkw(n, true_alpha, true_beta)
#'
#' # Obtain starting values
#' start_vals <- gkwgetstartvalues(data, family = "kw", n_starts = 3)
#'
#' # Maximum likelihood estimation with analytical gradient
#' fit <- optim(
#'   par = start_vals,
#'   fn = llkw, # Negative log-likelihood
#'   gr = grkw, # Analytical gradient
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE,
#'   control = list(maxit = 500)
#' )
#'
#' # Extract results
#' mle <- fit$par
#' se <- sqrt(diag(solve(fit$hessian)))
#'
#' # Construct confidence intervals (95%)
#' ci <- data.frame(
#'   Parameter = c("alpha", "beta"),
#'   True = c(true_alpha, true_beta),
#'   MLE = mle,
#'   SE = se,
#'   Lower = mle - 1.96 * se,
#'   Upper = mle + 1.96 * se
#' )
#'
#' print(ci, digits = 4)
#'
#' # Goodness-of-fit diagnostic
#' x_grid <- seq(0.001, 0.999, length.out = 200)
#' fitted_dens <- dkw(x_grid, mle[1], mle[2])
#' true_dens <- dkw(x_grid, true_alpha, true_beta)
#'
#' hist(data,
#'   breaks = 40, probability = TRUE,
#'   col = "lightgray", border = "white",
#'   main = "Kumaraswamy Distribution Fit",
#'   xlab = "Data", ylab = "Density", las = 1
#' )
#' lines(x_grid, fitted_dens, col = "red", lwd = 2, lty = 1)
#' lines(x_grid, true_dens, col = "blue", lwd = 2, lty = 2)
#' legend("topright",
#'   legend = c("Data", "Fitted", "True"),
#'   col = c("gray", "red", "blue"),
#'   lwd = c(8, 2, 2), lty = c(1, 1, 2),
#'   bty = "n"
#' )
#' grid(col = "gray80", lty = 2)
#'
#'
#' ## Example 4: Model Selection Using Information Criteria
#'
#' # Generate data from Exponentiated Kumaraswamy
#' set.seed(456)
#' n <- 1500
#' data <- rekw(n, alpha = 2, beta = 3, lambda = 1.5)
#'
#' # Define competing models
#' models <- list(
#'   Beta = list(
#'     ll = function(par) llbeta(par, data),
#'     gr = function(par) grbeta(par, data),
#'     start = gkwgetstartvalues(data, family = "beta", n_starts = 2),
#'     k = 2
#'   ),
#'   Kw = list(
#'     ll = function(par) llkw(par, data),
#'     gr = function(par) grkw(par, data),
#'     start = gkwgetstartvalues(data, family = "kw", n_starts = 2),
#'     k = 2
#'   ),
#'   EKw = list(
#'     ll = function(par) llekw(par, data),
#'     gr = function(par) grekw(par, data),
#'     start = gkwgetstartvalues(data, family = "ekw", n_starts = 2),
#'     k = 3
#'   ),
#'   Mc = list(
#'     ll = function(par) llmc(par, data),
#'     gr = function(par) grmc(par, data),
#'     start = gkwgetstartvalues(data, family = "mc", n_starts = 2),
#'     k = 3
#'   )
#' )
#'
#' # Fit all models
#' results <- lapply(names(models), function(name) {
#'   m <- models[[name]]
#'   fit <- optim(par = m$start, fn = m$ll, gr = m$gr, method = "BFGS")
#'
#'   loglik <- -fit$value
#'   aic <- -2 * loglik + 2 * m$k
#'   bic <- -2 * loglik + m$k * log(n)
#'
#'   data.frame(
#'     Model = name,
#'     k = m$k,
#'     LogLik = round(loglik, 2),
#'     AIC = round(aic, 2),
#'     BIC = round(bic, 2),
#'     stringsAsFactors = FALSE
#'   )
#' })
#'
#' # Combine and sort by AIC
#' comparison <- do.call(rbind, results)
#' comparison <- comparison[order(comparison$AIC), ]
#' rownames(comparison) <- NULL
#'
#' print(comparison)
#' cat("\nBest model by AIC:", comparison$Model[1], "\n")
#' cat("Best model by BIC:", comparison$Model[which.min(comparison$BIC)], "\n")
#' }
#'
#' @keywords internal
"_PACKAGE"
