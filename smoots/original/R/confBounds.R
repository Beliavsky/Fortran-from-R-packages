#' Asymptotically Unbiased Confidence Bounds
#'
#'@param obj an object returned by either \code{\link{msmooth}},
#'\code{\link{tsmooth}} or \code{\link{dsmooth}}.
#'@param alpha the confidence level; a single numeric value between \code{0}
#'and \code{1}; \code{0.95} is the default.
#'@param plot a logical value; for \code{plot = TRUE}, the default, a plot is
#'created.
#'@param p the order of polynomial used for the parametric polynomial
#'regression that is conducted as a benchmark for the trend function;
#'must satisfy \eqn{0 \leq}{0 \le} \code{p} \eqn{\leq 3}{\le 3}; set to
#'\code{1} by default; is irrelevant, if a derivative of the trend of order
#'greater than zero is being analyzed.
#'@param showPar set to \code{TRUE}, if the parametric fitted values are to be
#'shown against the unbiased estimates and the confidence bounds for
#'\code{plot = TRUE}; the default is \code{TRUE}.
#'@param rescale a single logical value; is set to \code{TRUE} by default;
#'if the output of a derivative estimation process is passed to \code{obj} and
#'if \code{rescale = TRUE}, the estimates and confidence bounds will be
#'rescaled according to \code{x} for the plot (see also the details on the
#'parameter \emph{...}); the numerical output stays unchanged.
#'@param ... further arguments that can be passed to the \code{plot} function;
#'if an argument \code{x} with time points is not given by the user,
#'\code{x = 1:length(obj$ye)} is used per default for the observation time
#'points.
#'
#'@details
#'This function is part of the \code{smoots} package and was implemented under
#'version 1.1.0. The underlying theory is based on the additive nonparametric
#'regression function
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series, \eqn{x_t} is the rescaled time
#'on the interval \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth trend function and
#'\eqn{\epsilon_t} are stationary errors with \eqn{E(\epsilon_t) = 0} and
#'short-range dependence.
#'
#'The purpose of this function is the estimation of reasonable confidence
#'intervals for the nonparametric trend function and its derivatives. The
#'optimal bandwidth minimizes the Asymptotic Mean Integrated Squared Error
#'(AMISE) criterion, however, local polynomial estimates are (usually) biased.
#'The bias is then (approximately)
#'\deqn{\frac{h^{k - v} m^{(k)}(x) \beta_{(\nu, k)}}{k!},}{h^[k-v]
#'m^(k)(x)\beta_(\nu,k) / k!,}
#'where \eqn{p} is the order of the local polynomials, \eqn{k = p + 1} is the
#'order of the asymptotically equivalent kernel, \eqn{\nu} is the order of the
#'of the trend function's derivative, \eqn{m^(v)} is the \eqn{\nu}-th order
#'derivative of the trend function and \eqn{\beta_{(\nu, k)} = \int_{-1}^{1}
#'u^k K_{(\nu, k)}(u) du}{\beta_(\nu, k) = int_[-1]^[1] u^k K_(\nu, k)(u) du}.
#'\eqn{K_{(\nu, k)}(u)}{K_(\nu, k)(u)} is the \eqn{k}-th order asymptotically
#'equivalent kernel function for estimating \eqn{m^{(\nu)}}{m^(\nu)}.
#'A renewed estimation with an adjusted bandwidth
#'\eqn{h_{ub} = o(n^{-1 / (2k + 1)})}{h_[ub] = o[n^[-1 / (2k + 1)]]}, i.e., a
#'bandwidth with a smaller order than the optimal bandwidth, is conducted.
#'\eqn{h = h_{A}^{(2k + 1) / (2k)}}{h = (h_[A])^[(2k + 1) / (2k)]}, where
#'\eqn{h_{A}}{h_[A]} is the optimal bandwidth, is implemented.
#'
#'Following this idea, we have that
#'\deqn{\sqrt{nh}[m^{(\nu)}(x) - \hat{m}^{(\nu)}(x)]}{(nh)^[1/2] (m^(v)(x) -
#'hat(m)^(v)(x))}
#'converges to
#'\deqn{N(0,2\pi c_f R(x))}
#'in distribution, where \eqn{2\pi c_f} is the sum of autocovariances.
#'Consequently, the trend (or derivative) estimates are asymptotically unbiased
#'and normally distributed.
#'
#'To make use of this function, an object of class \code{smoots} can be given
#'as input that was created by either \code{\link{msmooth}},
#'\code{\link{tsmooth}} or \code{\link{dsmooth}}. Based on the optimal
#'bandwidth saved within \code{obj}, an adjustment to the bandwidth is made so
#'that the estimates following the adjusted bandwidth are (relatively)
#'unbiased.
#'
#'Based on the input argument \code{alpha}, the level of confidence between
#'\code{0} and \code{1}, the respective confidence bounds are calculated for
#'each observation point.
#'
#'From the input argument \code{obj}, the order of derivative is automatically
#'obtained. By means of the argument \code{p}, an order of polynomial is
#'selected for a parametric regression of the trend function. This is only
#'meaningful, if the trend (and not its derivatives) is analyzed. Otherwise,
#'the argument is automatically dropped by the function. Furthermore, if
#'\code{plot = TRUE}, a plot of the unbiased trend (or derivative) estimates
#'alongside the confidence bounds is created. If also \code{showPar = TRUE},
#'the estimated parametric trend (or parametric constant value for the
#'derivatives) is added to the confidence bound plot for comparison.
#'
#'NOTE:
#'
#'The values that are returned by the function are obtained with respect to
#'the rescaled time points on the interval \eqn{[0, 1]}. While the plot can be
#'adjusted and rescaled by means of a given vector with the actual time points,
#'the numeric output is not rescaled. For this purpose we refer the user to
#'the \code{\link{rescale}} function of the \code{smoots} package.
#'
#'This function implements C++ code by means of the
#'\code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages for
#'better performance.
#'
#'@md
#'
#'@export
#'
#'@importFrom graphics lines polygon
#'@importFrom stats as.formula lm qnorm
#'
#'@references
#' Beran, J. and Feng, Y. (2002). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. Unpublished.
#'
#' Feng, Y., Gries, T., Fritz, M., Letmathe, S. and Schulz, D. (2020).
#' Diagnosing the trend and bootstrapping the forecasting intervals using a
#' semiparametric ARMA. Discussion Paper. Paderborn University. Unpublished.
#'
#'@return
#'A plot is created in the plot window and a list with different components
#'is returned.
#'
#'\describe{
#'\item{alpha}{a numeric vector of length 1; the level of confidence; input
#'argument.}
#'\item{b.ub}{a numeric vector with one element that represents the adjusted
#'bandwidth for the unbiased trend estimation.}
#'\item{p.estim}{a numeric vector with the estimates following the parametric
#'regression defined by \code{p} that is conducted as a benchmark for the trend
#'function; for the trend's derivatives or for \code{p = 0}, a constant value
#'is the benchmark; the values are obtained with respect to the rescaled time
#'points on the interval \eqn{[0, 1]}.}
#'\item{n}{the number of observations.}
#'\item{np.estim}{a data frame with the three (numeric) columns \strong{ye.ub},
#'\strong{lower} and \strong{upper}; in \strong{ye.ub} the unbiased trend
#'estimates, in \strong{lower} the lower confidence bound and in \strong{upper}
#'the upper confidence bound can be found; the values are obtained with respect
#'to the rescaled time points on the interval \eqn{[0, 1]}.}
#'\item{v}{the order of the trend's derivative considered for the test.}
#'}
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'log_gdp <- log(smoots::gdpUS$GDP)
#'est <- msmooth(log_gdp)
#'confBounds(est)
#'

confBounds <- function(obj, alpha = 0.95, p = c(0, 1, 2, 3), plot = TRUE,
                showPar = TRUE, rescale = TRUE, ...) {

  # The plot is defined by the input object of class "smoots"
  if (!inherits(obj, "smoots")) {
    stop("Input object 'obj' not recognized. It needs to be of class",
      " 'smoots'.")
  }
  if (!(attr(obj, "function") %in% c("msmooth", "tsmooth", "dsmooth"))) {
    stop("Input object 'obj' of class 'smoots' must have been returned by
    either 'msmooth()', 'tsmooth()' or 'dsmooth()'.")
  } else {
    n.plot <- obj[["v"]] + 1
    if (attr(obj, "function") %in% c("msmooth", "tsmooth")) {
      bb <- obj$bb
    } else {
      bb <- 1
    }
  }
  if (length(alpha) != 1 || is.na(alpha) || !is.numeric(alpha) ||
      !(alpha > 0 && alpha < 1)) {
    stop("The argument 'alpha' must be 0 < alpha < 1.")
  }
  if (n.plot != 1) {
    p <- 0
  }
  if (!(length(p) %in% c(1, 4)) || !all(!is.na(p)) || !is.numeric(p)) {
    stop("The argument 'p' must be a single integer value with ",
         "0 <= p <= 3.")
  }
  p <- floor(p)
  if (all(p == (0:3))) p <- 1
  if (length(p) != 1 || !(p %in% (0:3))) {
    stop("The argument 'p' must be a single integer value with ",
         "0 <= p <= 3.")
  }
  if (length(plot) != 1 || is.na(plot) || !(plot %in% c(TRUE, FALSE))) {
    stop("The argument 'plot' must be a single logical value (TRUE or FALSE).")
  }
  if (length(showPar) != 1 || is.na(plot) || !(showPar %in% c(TRUE, FALSE))) {
    stop("The argument 'showPar' must be a single logical value (TRUE or ",
      "FALSE).")
  }

  if (length(rescale) != 1 || is.na(rescale) ||
    !(rescale %in% c(TRUE, FALSE))) {
    stop("The argument 'rescale' must be a single logical value (TRUE or ",
      "FALSE).")
  }

  # Set default value for p. The default is p = 1 (Linear Regression)
  if (n.plot == 1) {
    pp <- obj$p
  } else {
    pp <- obj$pp
    if (obj$pp == 1) {
      alg <- "A"
    } else {
      alg <- "B"
    }
  }

  # Call optimal bandwidth and order of polynomial from smoots object.
  b.opt <- obj$b0
  p.obj <- obj$p
  # Define the order of kernel.
  k <- p.obj + 1
  # Get the original observation series.
  y <- obj$orig
  # Calculate the (asymptotically) unbiased trend estimates.
  mu <- obj$mu

  # Calculate the adjusted bandwidth for the unbiased trend estimation.
  b.ub <- b.opt ^ ((2 * k + 1) / (2 * k))
  if (obj[["bvc"]] == "Y") {
    bv_func <- lookup$bvc_lookup[[as.character(pp), mu + 1]]
  } else {
    bv_func <- function(b) b
  }
  cf0Method <- c("NP" = function(X) cf0Cpp(X)$cf0.LW,
    "AR" = function(X) cf0.AR.est(X)$cf0.AR,
    "MA" = function(X) cf0.MA.est(X)$cf0.MA,
    "ARMA" = function(X) cf0.ARMA.est(X)$cf0.ARMA)
  cf0Calc <- cf0Method[[obj[["Mcf"]]]]

  ub.est <- gsmoothCalc2Cpp(y, n.plot - 1, p.obj, mu, b.ub, bb)
  ye.ub <- c(ub.est$ye)

  n <- length(y)
  t <- 1:n


  # Conduct a parametric regression for p > 0.
  if (p %in% 1:3) {
    par.formula <- as.formula(paste0("y ~ ", paste0("I(t ^", 1:p,
                              " )", collapse = " + ")))
    par.est <- lm(par.formula)
    ye.par <- par.est$fitted.values
  } else if (p == 0) {
    if (n.plot == 1) {
      ye.par <- mean(y)
    } else if (n.plot == 2) {
      ye.par <- lm(y ~ I(t / n))$coef[[2]]
    } else if (n.plot == 3) {
      ye.par <- 0
    }
  } else {
    ye.par <- NA
  }

  if (n.plot == 1) {
    bv <- bv_func(b.ub)
    ye.ub.cf <- c(gsmoothCalcCpp(y, n.plot - 1, pp, mu, bv, bb))
    cf <- cf0Calc(y - ye.ub.cf)
  } else {
    pilot.est <- msmoothCalc(y, p = pp, mu = mu, bStart = obj$bStart.p,
      alg = alg, method = "lpr")
    b.pilot.ub <- pilot.est$b0 ^ ((2 * (pilot.est$p + 1) + 1) /
      (2 * (pilot.est$p + 1)))
    bv.pilot <- bv_func(b.pilot.ub)
    pilot.ub.est.cf <- c(gsmoothCalcCpp(y, 0, pp, mu, bv.pilot, 1))
    cf <- cf0Calc(y - pilot.ub.est.cf)
  }

  # Get the weighting system from the input smoots object.
  ws <- ub.est$ws

  R.x.sub <- rowSums(ws^2)
  n.ws <- length(R.x.sub)
  boundary <- (n.ws - 1) / 2
  R.x <- c(R.x.sub[1:boundary], rep(R.x.sub[boundary + 1], n - 2 * boundary),
    R.x.sub[(boundary + 2):n.ws])
  norm.value <- qnorm(alpha + ((1 - alpha) / 2), mean = 0, sd = 1)

  est.norm <- norm.value * sqrt(cf * R.x)
  conf1 <- ye.ub - est.norm
  conf2 <- ye.ub + est.norm
  v <- n.plot - 1

  outp <- list(np.estim = data.frame(ye.ub = ye.ub, lower = conf1,
    upper = conf2), p.estim = ye.par, b.ub = b.ub, alpha = alpha,
    v = v, n = obj$n)

  if (plot == TRUE) {
    dots <- list(...)
    if (is.null(dots[["col"]])) dots[["col"]] <- "red"
    if (is.null(dots[["main"]])) {
      if (v == 0) {
        title.der <- "zeroth"
      } else if (v == 1) {
        title.der <- "first"
      } else if (v == 2) {
        title.der <- "second"
      }
      color.n <- color.name(dots[["col"]])
      if (p != 0 && showPar == TRUE) {
        plot.title <- paste0("The parametric (blue) and the ",
          "nonparametric unbiased ", title.der, " derivative (",
          color.n, ")\nof the trend together ", "with ", alpha * 100,
          "%-confidence intervals (grey)")
      } else {
        plot.title <- paste0("The nonparametric unbiased ",
          title.der, " derivative (", color.n, ")\nof the trend together ",
          "with ", alpha * 100, "%-confidence intervals (grey)")
      }
      dots[["main"]] <- plot.title
    }

    if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "Values"
    if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "Time"
    x_coord <- dots[["x"]]
    if (is.null(x_coord)) x_coord <- t
    dots[["x"]] <- 1
    if (rescale == TRUE && v > 0) {
      ye.ub <- smoots::rescale(ye.ub, x = x_coord, v = v)
      conf1 <- smoots::rescale(conf1, x = x_coord, v = v)
      conf2 <- smoots::rescale(conf2, x = x_coord, v = v)
      ye.par <- smoots::rescale(ye.par, x = x_coord, v = v)
    }
    if (showPar == TRUE && p == 0) {
      ye.par <- rep(ye.par, n)
    }
    if (is.null(dots[["lty"]])) dots[["lty"]] <- 1
    if (is.null(dots[["lwd"]])) dots[["lwd"]] <- 1
    type.sub <- dots[["type"]]
    if (is.null(type.sub)) type.sub <- "l"
    dots[["type"]] <- "n"
    if (is.null(dots[["xlim"]])) dots[["xlim"]] <- c(x_coord[1], x_coord[n])
    x.low <- sum((x_coord) < dots$xlim[1]) + 1
    x.up <- sum((x_coord) <= dots$xlim[2])
    if (is.null(dots[["ylim"]])) {
      if (showPar == TRUE) {
          dots[["ylim"]] <- c(min(conf1[x.low:x.up], ye.par[x.low:x.up]),
                         max(conf2[x.low:x.up], ye.par[x.low:x.up]))
      } else {
        dots[["ylim"]] <- c(min(conf1[x.low:x.up]),
          max(conf2[x.low:x.up]))
      }
    }
    do.call(what = graphics::plot, args = dots)
    polygon(c(x_coord, rev(x_coord)), c(conf1, rev(conf2)),
            col = col.alpha("gray", alpha = 0.8), border = NA)
    dots[["x"]] <- x_coord
    dots[["y"]] <- ye.ub
    dots[["type"]] <- type.sub
    do.call(what = lines, args = dots)
    if (showPar == TRUE & !is.na(ye.par[1])) {
      lwd.sub <- dots[["lwd"]]
      lines(x_coord, ye.par, lty = 1, col = "deepskyblue4", lwd = lwd.sub)
    }
  }
  class(outp) <- "smoots"
  attr(outp, "function") <- "confBounds"
  outp
}

