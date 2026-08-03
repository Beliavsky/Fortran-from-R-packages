#' @title Density Quantile Function
#' @description Estimate the density quantile function
#' @param y Vector of dependent data
#' @param x Matrix of covariates
#' @param u Quantile residuals
#' @param alpha Probability level
#' @param sparsity iid or ind
#' @param bandwidth_estimator Bofinger, Chamberlain or Hall-Sheather
#' @references
#' For the iid and nid method, see Koenker (1994), and Hendricks and Koenker (1992).
#' For the bandwidth types, see Bofinger (1975), Chamberlain (1994), and Hall and Sheather(1988).
#' @keywords internal
#' @export
density_quantile_function <- function(y, x, u, alpha, sparsity, bandwidth_estimator) {
  n <- nrow(x)
  k <- ncol(x)
  eps <- .Machine$double.eps^(2/3)

  # Get the bandwidth
  if (bandwidth_estimator == "Bofinger") {
    bandwidth <- n^(-1/5) * ((9/2 * stats::dnorm(stats::qnorm(alpha))^4)/(2 * stats::qnorm(alpha)^2 + 1)^2)^(1/5)
  } else if (bandwidth_estimator == "Chamberlain") {
    tau <- 0.05
    bandwidth <- stats::qnorm(1 - alpha/2) * sqrt(tau * (1 - tau)/n)
  } else if (bandwidth_estimator == "Hall-Sheather") {
    tau <- 0.05
    bandwidth <- n^(-1/3) * stats::qnorm(1 - tau/2)^(2/3) * ((3/2 * stats::dnorm(stats::qnorm(alpha))^2)/(2 * stats::qnorm(alpha)^2 + 1))^(1/3)
  } else {
    stop("Not a valid bandwith method!")
  }

  # Compute the density
  if (sparsity == "iid") {
    # Koenker (1994)
    h <- max(k + 1, ceiling(n * bandwidth))
    ir <- (k + 1):(h + k + 1)
    ord.resid <- sort(u[order(abs(u))][ir])
    xt <- ir/(n - k)
    density <- 1/suppressWarnings(quantreg::rq(ord.resid ~ xt)$coef[2])
    density <- rep(as.numeric(density), n)
  } else if (sparsity == "nid") {
    # Hendricks and Koenker (1992)
    bu <- quantreg::rq(y ~ x - 1, tau = alpha + bandwidth)$coef
    bl <- quantreg::rq(y ~ x - 1, tau = alpha - bandwidth)$coef
    density <- pmax(0, 2 * bandwidth/(x %*% (bu - bl) - eps))
  } else {
    stop("Not a valid density quantile function estimator!")
  }

  density
}

#' @title Conditional truncated variance
#' @description Estimate the variance of y given x and given y <= 0.
#' If approach is:
#' \enumerate{
#'   \item ind -  Variance of all y where y <= 0.
#'   \item scl_N or scl_sp - Assumes a location-scale model:
#'   y = x'b + (x'g)e. First, it estimates b and g using PMLE.
#'   Then, it computes the conditional truncated variance by integrating the truncated density of y.
#' }
#' @param y Vector of dependent data
#' @param x Matrix of covariates including the intercept
#' @param approach ind, scl_N or scl_sp
#' @keywords internal
#' @export
conditional_truncated_variance <- function(y, x, approach) {
  if (sum(y <= 0) <= 2) {
    stop("Not enough negative quantile residuals!")
  }

  if (approach == "ind") {
    cv <- rep(stats::var(y[y <= 0]), length(y))
  } else {
    cv <- tryCatch({
      # Get conditional mean and sigma
      mu_sigma <- conditional_mean_sigma(y, x)
      mu <- mu_sigma$mu
      sigma <- mu_sigma$sigma

      # Truncated conditional variance
      if (approach == "scl_N") {
        beta <- -mu / sigma
        beta[beta < -30] <- -30
        cv <- sigma^2 * (1 - beta * stats::dnorm(beta)/stats::pnorm(beta) -
                           (stats::dnorm(beta)/stats::pnorm(beta))^2)
      } else if (approach == "scl_sp") {
        # Kernel density estimate of the standardized residuals
        df <- stats::density((y - mu) / sigma, bw = "SJ")
        lower_int_boundary <- min(df$x)
        pdf <- stats::approxfun(df, yleft = 0, yright = 0)

        # Truncation points
        beta <- -mu / sigma

        # Integrate the truncated pdf by the trapezoidal rule
        div <- 1000
        h <- (max(beta) - lower_int_boundary) / (div - 1)
        b_approx <- seq(lower_int_boundary, max(beta) + h, h)
        midpoint <- b_approx[-div] + h/2
        y0 <- pdf(b_approx)
        y1 <- b_approx * y0
        y2 <- b_approx^2 * y0
        cb <- cumsum(y0[-1] + y0[-div]) / 2 * h
        m1 <- cumsum(y1[-1] + y1[-div]) / 2 * h
        m2 <- cumsum(y2[-1] + y2[-div]) / 2 * h
        cv_approx <- m2 / cb - (m1 / cb)^2
        cv_approx[cv_approx < 0] <- NA

        # Approximate the conditional truncated variance
        cv <- sigma^2 * stats::approx(x = midpoint, y = cv_approx, xout = beta)$y
      }
      if (any(is.na(cv)) | any(!is.finite(cv)) | any(cv < 0)) stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit the ", approach, " estimator, switching to the ind approach!"))
      rep(stats::var(y[y <= 0]), length(y))
    })
  }

  cv
}

#' @title Cumulative Density Function at Quantile
#' @description Returns the cumulative density function evaluated at quantile predictions.
#' For a correctly specified model this should yield a value close to the quantile level.
#' @param y Vector of dependent data
#' @param x Matrix of covariates including the intercept
#' @param q Vector of quantile predictions
#' @keywords internal
#' @export
cdf_at_quantile <- function(y, x, q) {
  # Get conditional mean and sigma
  mu_sigma <- conditional_mean_sigma(y, x)
  mu <- mu_sigma$mu
  sigma <- mu_sigma$sigma

  # Empirical CDF of standardized data
  cdf <- function(x) stats::ecdf((y - mu) / sigma)(x)

  # CDF of standardized quantile predictions
  z <- (q - mu) / sigma
  cdf(z)
}

#' @title Conditional Mean and Sigma
#' @description Estimate the conditional mean and sigma of the dependent data conditional on covariates x
#' @param y Vector of dependent data
#' @param x Matrix of covariates including the intercept
#' @keywords internal
#' @export
conditional_mean_sigma <- function(y, x) {
  # Starting values and ensure positive fitted standard deviations
  fit1 <- stats::lm(y ~ x - 1)
  fit2 <- stats::lm(abs(fit1$residuals) ~ x - 1)
  fit2$coefficients[1] <- fit2$coefficients[1] - min(0.001, min(fit2$fitted.values))
  b0 <- c(fit1$coefficients, fit2$coefficients)

  # Estimate the model under normality
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    ifelse(all(sigma > 0), -sum(stats::dnorm(x=y, mean=mu, sd=sigma, log=TRUE)), NA)
  }
  fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="BFGS"), silent=TRUE)
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="Nelder-Mead",
                            control=list(maxit=10000)), silent=TRUE)
  }
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    warning("Conditional variance estimation: cannot optimize model, returning starting values instead")
    b <- b0
  } else {
    b <- fit$par  
  }
  
  # Estimated means and standard deviations
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])

  list(mu = mu, sigma = sigma)
}
