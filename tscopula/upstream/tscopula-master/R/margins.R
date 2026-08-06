#' Marginal model for time series
#'
#' Class of objects for marginal models for stationary time series. The
#' object is given a name and there must exist functions pname, qname,
#' dname and rname. As well as the parameters of the distribution, dname must have the
#' logical argument log specifying whether log density should be computed.
#'
#' @slot name name of the marginal model class.
#' @slot pars a numeric vector containing the named parameters of the distribution
#' which are passed as arguments to pname, qname, dname and rname.
#'
#' @export
#'
#' @examples
#' new("margin", name = "gauss", pars = c(mu = 0, sigma = 1))
setClass("margin", slots = list(
  name = "character",
  pars = "numeric"
))

#' @describeIn margin Coef method for margin class
#'
#' @param object an object of the class.
#'
#' @export
#'
#'
setMethod("coef", "margin", function(object) {
  if (is(object, "marginfit")) {
    object <- object@margin
  }
  object@pars
})


#' Constructor function for margin
#'
#' @param name character string giving name of distribution
#' @param pars parameters of the distribution
#'
#' @return An object of class \linkS4class{margin}.
#' @export
#'
#' @examples
#' margin("sst")
margin <- function(name, pars = NULL) {
  if (name == "norm")
    name <- "gauss"
  pfunc <- eval(parse(text = paste("p", name, sep = "")))
  defaults <- unlist(formals(pfunc)[-1])
  if (!is.null(pars)) {
    parnames <- names(pars)
    parset <- parnames[parnames %in% names(defaults)]
    if (length(parset) == 0) {
      stop("Unrecognized parameter names")
    }
    defaults[parset] <- pars[parset]
  }
  new("margin", name = name, pars = defaults)
}

#' Compute CDF of marginal model
#'
#' Compute the cumulative distribution function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param q vector of values at which CDF should be computed.
#'
#' @return A vector of values for the CDF.
#' @export
#'
#' @examples
#' margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
#' pmarg(margmod, c(-2, 0, 2))
pmarg <- function(x, q) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("p", x@name, sep = "")))
  do.call(func, append(x@pars, list(q = q)))
}

#' Compute quantiles of marginal model
#'
#' Compute the quantile function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param p vector of probabilities for which quantiles should be computed.
#'
#' @return A vector of values for the quantile function.
#' @export
#'
#' @examples
#' margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
#' qmarg(margmod, c(0.05, 0.5, 0.95))
qmarg <- function(x, p) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("q", x@name, sep = "")))
  do.call(func, append(x@pars, list(p = p)))
}

#' Compute density of marginal model
#'
#' Compute the density function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param y vector of values for which density should be computed.
#' @param log logical variable specifying whether log density should be
#' returned.
#'
#' @return A vector of values for the density.
#' @export
#'
#' @examples
#' margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
#' dmarg(margmod, c(-2, 0, 2), log = TRUE)
dmarg <- function(x, y, log = FALSE) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("d", x@name, sep = "")))
  do.call(func, append(x@pars, list(x = y, log = log)))
}

#' Gaussian distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param log flag for log density.
#' @name gauss
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#'
#' @rdname gauss
#' @export
dgauss <- function(x, mu = 0, sigma = 1, log = FALSE){
  if (sigma < 0)
    return(NA)
  dnorm(x, mean = mu, sd = sigma, log = log)
}
#' @rdname gauss
#' @export
pgauss <- function(q, mu = 0, sigma = 1){
  if (sigma < 0)
    return(NA)
  pnorm(q, mean = mu, sd = sigma)
}
#' @rdname gauss
#' @export
qgauss <- function(p, mu = 0, sigma = 1){
  if (sigma < 0)
    return(NA)
  qnorm(p, mean = mu, sd = sigma)
}
#' @rdname gauss
#' @export
rgauss <- function(n, mu = 0, sigma = 1){
  rnorm(n, mean = mu, sd = sigma)
}

#' Centred Gaussian distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param sigma scale parameter.
#' @param log flag for log density.
#' @name gauss0
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#'
#' @rdname gauss0
#' @export
dgauss0 <- function(x, sigma = 1, log = FALSE){
  dgauss(x, mu = 0, sigma = sigma, log = log)
}
#' @rdname gauss0
#' @export
pgauss0 <- function(q, sigma = 1){
  pgauss(q, mu = 0, sigma = sigma)
}
#' @rdname gauss0
#' @export
qgauss0 <- function(p, sigma = 1){
  qgauss(p, mu = 0, sigma = sigma)
}
#' @rdname gauss0
#' @export
rgauss0 <- function(n, sigma = 1){
  rgauss(n, mu = 0, sigma = sigma)
}

#' Laplace distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter.
#' @param scale scale parameter.
#' @param log flag for log density.
#' @name laplace
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#'
#' @rdname laplace
#' @export
dlaplace <- function(x, mu = 0, scale = 1, log = FALSE){
  if (scale <= 0){
    return(NA)
  }
  y <- (x - mu)/scale
  tmp <- -abs(y) - log(2*scale)
  arg <- rep(NA, length(x))
  if (!log)
    tmp <- exp(tmp)
  tmp
}
#' @rdname laplace
#' @export
plaplace <- function(q, mu = 0, scale = 1){
  if (scale <= 0)
    return(NA)
  y <- (q - mu)/scale
  tmp <- exp(-abs(y))/2
  tmp[y >= 0] <- 1 - tmp[y >= 0]
  tmp
}
#' @rdname laplace
#' @export
qlaplace <- function(p, mu = 0, scale = 1){
  if (scale <= 0)
    return(NA)
  tmp <- rep(NA, length(p))
  tmp[p < 0.5] <- log(2*p[p < 0.5])
  tmp[p >= 0.5] <- -log(2 - 2*p[p >= 0.5])
  mu + scale*tmp
}
#' @rdname laplace
#' @export
rlaplace <- function(n, mu = 0, scale = 1){
  qlaplace(runif(n), mu, scale)
}

#' Centred Laplace distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param scale scale parameter.
#' @param log flag for log density.
#' @name laplace0
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#'
#' @rdname laplace0
#' @export
dlaplace0 <- function(x, scale = 1, log = FALSE){
 dlaplace(x, mu = 0, scale = scale, log = log)
}
#' @rdname laplace0
#' @export
plaplace0 <- function(q, scale = 1){
 plaplace(q, mu = 0 , scale = scale)
}
#' @rdname laplace0
#' @export
qlaplace0 <- function(p, scale = 1){
 qlaplace(p, mu = 0, scale = scale)
}
#' @rdname laplace0
#' @export
rlaplace0 <- function(n, scale = 1){
  qlaplace(runif(n), mu = 0, scale)
}


#' Skew Laplace distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter.
#' @param scale scale parameter.
#' @param gamma skewness parameter.
#' @param log flag for log density.
#' @name slaplace
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#'
#' @rdname slaplace
#' @export
dslaplace <- function(x, mu = 0.05, scale = 1, gamma = 1, log = FALSE){
  dsdoubleweibull(x, mu = mu, shape = 1, scale = scale, gamma = gamma, log = log)
}
#' @rdname slaplace
#' @export
pslaplace <- function(q, mu = 0.05, scale = 1, gamma = 1){
  psdoubleweibull(q, mu = mu, shape = 1, scale = scale, gamma = gamma)
}
#' @rdname slaplace
#' @export
qslaplace <- function(p, mu = 0.05, scale = 1, gamma = 1){
  qsdoubleweibull(p, mu = mu, shape = 1, scale = scale, gamma = gamma)
}
#' @rdname slaplace
#' @export
rslaplace <- function(n, mu = 0.05, scale = 1, gamma = 1){
  qslaplace(runif(n), mu, scale, gamma)
}

#' Double Weibull distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter.
#' @param shape shape parameter.
#' @param scale scale parameter.
#' @param log flag for log density.
#' @name doubleweibull
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#' @rdname doubleweibull
#' @export
ddoubleweibull <- function(x, mu = 0.05, shape = 1, scale = 1, log = FALSE){
  dsdoubleweibull(x, mu = mu, shape = shape, scale = scale, gamma = 1, log = log)
}
#' @rdname doubleweibull
#' @export
pdoubleweibull <- function(q, mu = 0.05, shape = 1, scale = 1){
  psdoubleweibull(q, mu = mu, shape = shape, scale = scale, gamma = 1)
}
#' @rdname doubleweibull
#' @export
qdoubleweibull <- function(p, mu = 0.05, shape = 1, scale = 1){
  qsdoubleweibull(p, mu = mu, shape = shape, scale = scale, gamma = 1)
}
#' @rdname doubleweibull
#' @export
rdoubleweibull <- function(n, mu = 0.05, shape = 1, scale = 1){
  qdoubleweibull(runif(n), mu = mu, shape = shape, scale = scale)
}

#' Skew double Weibull distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu location parameter.
#' @param shape shape parameter.
#' @param scale scale parameter.
#' @param gamma skewness parameter.
#' @param log flag for log density.
#' @name sdoubleweibull
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#' @rdname sdoubleweibull
#' @export
dsdoubleweibull <- function(x, mu = 0.05, shape = 1, scale = 1, gamma = 1, log = FALSE)
{
  if ((scale <= 0) | (shape <= 0)) {
    return(NA)
  }
  arg <- rep(NA, length(x))
  y <- (x - mu)/scale
  arg[y < 0] <- gamma * abs(y[y < 0])
  arg[y >= 0] <- y[y >=0]/gamma
  tmp <- (shape - 1)*log(arg) - arg^shape + log(shape) -log(scale) - log(gamma + 1/gamma)
  if (!log)
    tmp <- exp(tmp)
  tmp
}
#' @rdname sdoubleweibull
#' @export
psdoubleweibull <- function(q, mu = 0.05, shape = 1, scale = 1, gamma = 1)
{
  arg <- rep(NA, length(q))
  y <- (q - mu)/scale
  arg[y < 0] <- gamma * abs(y[y < 0])
  arg[y >= 0] <- y[y >=0]/gamma
  cumy <- exp(-arg^shape)/(1+gamma^2)
  cumy[y > 0] <- (1-(gamma^2)*cumy[y > 0])
  cumy
}
#' @rdname sdoubleweibull
#' @export
qsdoubleweibull <- function(p, mu = 0.05, shape = 1, scale = 1, gamma = 1)
{
  tmp <- rep(NA, length(p))
  wt <- (1+gamma^2)
  lower <- (p <= 1/wt)
  tmp[lower] <- mu - scale * ((-log(wt*p[lower]))^(1/shape))/gamma
  tmp[!lower] <- mu + scale * ((-log(wt*(1-p[!lower])/(wt-1)))^(1/shape))*gamma
  tmp
}
#' @rdname sdoubleweibull
#' @export
rsdoubleweibull <- function(n, mu = 0.05, shape = 1, scale = 1, gamma = 1)
{
  qsdoubleweibull(runif(n), mu, shape, scale, gamma)
}

#' Student t distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param df degrees of freedom.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param log flag for log density.
#' @name st
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#' @rdname st
#' @export
pst <- function(q, df = 10, mu = 0, sigma = 1) {
  pt((q - mu) / sigma, df)
}
#' @rdname st
#' @export
qst <- function(p, df, mu, sigma) {
  qt(p, df) * sigma + mu
}
#' @rdname st
#' @export
dst <- function(x, df, mu, sigma, log = FALSE) {
  if ((sigma < 0) | (df < 0)) {
    return(NA)
  }
  dens <- dt((x - mu) / sigma, df = df, log = log)
  if (log) {
    return(dens - log(sigma))
  } else {
    return(dens / sigma)
  }
}
#' @rdname st
#' @export
rst <- function(n, df, mu, sigma) {
  rt(n, df) * sigma + mu
}

#' Centred Student t distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param df degrees of freedom.
#' @param sigma scale parameter.
#' @param log flag for log density.
#' @name st0
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#' @rdname st0
#' @export
pst0 <- function(q, df = 10, sigma = 1) {
  pst(q, df = df, mu = 0, sigma = sigma)
}
#' @rdname st0
#' @export
qst0 <- function(p, df, sigma) {
  qst(p, df = df, mu = 0, sigma = sigma)
}
#' @rdname st0
#' @export
dst0 <- function(x, df, sigma, log = FALSE) {
  dst(x, df = df, mu = 0, sigma = sigma, log = log)
}
#' @rdname st0
#' @export
rst0 <- function(n, df, sigma) {
  rt(n, df = df) * sigma
}

#' Skew Student t distribution
#'
#' @param x vector of values.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param df degrees of freedom.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param gamma skewness parameter.
#' @param log flag for log density.
#' @name sst
#' @return A vector of density, distribution function, quantile or random values.
NULL
#> NULL
#' @rdname sst
#' @export
psst <- function(q, df = 10, gamma = 1, mu = 0, sigma = 1) {
  result <- rep(NA, length(q))
  x <- (q - mu) / sigma
  result[x < 0] <- 2 / (gamma^2 + 1) * pt(gamma * x[x < 0], df)
  result[x >= 0] <- 1 / (gamma^2 + 1) + 2 / (1 + (1 / gamma^2)) * (pt(x[x >= 0] / gamma, df) -
    1 / 2)
  result
}
#' @rdname sst
#' @export
qsst <- function(p, df, gamma, mu, sigma) {
  result <- rep(NA, length(p))
  probzero <- 1 / (gamma^2 + 1)
  result[p < probzero] <- 1 / gamma * qt(((gamma^2 + 1) * p[p < probzero]) / 2, df)
  result[p >= probzero] <- gamma * qt((1 + 1 / gamma^2) / 2 * (p[p >= probzero] - probzero) +
    1 / 2, df)
  result * sigma + mu
}
#' @rdname sst
#' @export
dsst <- function(x, df, gamma, mu, sigma, log = FALSE) {
  if ((sigma < 0) | (df < 0)) {
    return(NA)
  }
  result <- rep(NA, length(x))
  x <- (x - mu) / sigma
  result[x < 0] <- dt(gamma * x[x < 0], df, log = log)
  result[x >= 0] <- dt(x[x >= 0] / gamma, df, log = log)
  if (log) {
    return(result + log(2 / (gamma + 1 / gamma)) - log(sigma))
  } else {
    return(result * (2 / (gamma + 1 / gamma)) / sigma)
  }
}
#' @rdname sst
#' @export
rsst <- function(n, df, gamma, mu, sigma) {
  p <- runif(n)
  result <- rep(NA, n)
  probzero <- 1 / (gamma^2 + 1)
  result[p < probzero] <- 1 / gamma * qt(((gamma^2 + 1) * p[p < probzero]) / 2, df)
  result[p >= probzero] <- gamma * qt((1 + 1 / gamma^2) / 2 * (p[p >= probzero] - probzero) +
    1 / 2, df)
  result * sigma + mu
}

#' @describeIn margin Simulation method for margin class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
#' sim(margmod, n = 500)
setMethod("sim", c(object = "margin"), function(object, n = 1000) {
  func <- eval(parse(text = paste("r", object@name, sep = "")))
  do.call(func, append(object@pars, list(n = n)))
})

#' Fitted marginal model for time series
#'
#' @slot margin an object of class \linkS4class{margin}.
#' @slot data numeric vector or time series of data.
#' @slot fit a list containing details of the maximum likelihood fit.
#'
#' @export
#'
setClass("marginfit",
  contains = "margin",
  slots = list(
    margin = "margin",
    data = "ANY",
    fit = "list"
  )
)

#' Fit method for margin class
#'
#' @param x an object of class \linkS4class{margin}.
#' @param y a vector or time series of data.
#' @param tsoptions list of optional arguments:
#' hessian is logical variable specifying whether Hessian matrix should be returned;
#' start is vector od named starting values
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{marginfit}.
#' @export
#'
#' @examples
#' margmod <- margin("gauss", pars = c(mu = 0, sigma = 1))
#' data <- sim(margmod, n = 500)
#' fit(margmod, data)
setMethod("fit", c(x = "margin", y = "ANY"), function(x, y,
                                                      tsoptions = list(),
                                                      control = list()) {
  defaults <- list(hessian = FALSE, method = "Nelder-Mead")
  cdefaults <- list(maxit = 1000, warn.1d.NelderMead = FALSE)
  tsoptions <- setoptions(tsoptions, defaults)
  control <- setoptions(control, cdefaults)
  dens <- eval(parse(text = paste("d", x@name, sep = "")))
  objective <- function(theta, dens, y) {
    dx <- do.call(dens, append(theta, list(x = y, log = TRUE)))
    -sum(dx)
  }
  fit <- optim(x@pars,
    fn = objective,
    dens = dens,
    y = y,
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  x@pars <- fit$par
  new("marginfit", margin = x, data = y, fit = fit)
})

#' @describeIn margin Show method for margin class
#'
#' @param object an object of the class.
#'
#' @export
#'
#'
setMethod("show", "margin", function(object) {
  if (!(is(object, "marginfit"))){
    cat("name: ", object@name, "\n", sep = "")
    cat("parameters: \n")
    print(coef(object))
  }
  if (is(object, "marginfit")) {
    cat("name: ", object@margin@name, "\n", sep = "")
    ests <- object@fit$par
    if (is.element("hessian", names(object@fit))) {
      ses <- safe_ses(object@fit$hessian)
      ests <- rbind(ests, ses)
      dimnames(ests)[[1]] <- c("par", "se")
    }
    cat("estimates:\n")
    print(ests)
    cat("convergence status: ", object@fit$convergence, ", log-likelihood: ", -object@fit$value,
      "\n",
      sep = ""
    )
  }
})

#' @describeIn marginfit logLik method for marginfit class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("logLik", "marginfit", function(object) {
  ll <- -object@fit$value
  attr(ll, "nobs") <- length(object@data)
  attr(ll, "df") <- length(object@fit$par)
  class(ll) <- "logLik"
  ll
})

#' Plot method for marginfit class
#'
#' @param x an object of class \linkS4class{marginfit}.
#' @param bw logical variable specifying whether black-white options should be chosen.
#'
#' @return No return value, generates plot.
#' @export
#'
setMethod("plot", c(x = "marginfit", y = "missing"),
          function(x, bw = FALSE) {
            colchoice <- ifelse(bw, "gray50", "red")
              qus <- qmarg(x@margin, ppoints(x@data))
              plot(qus,
                   sort(as.numeric(x@data)),
                   xlab = "Theoretical",
                   ylab = "data"
              )
              abline(0, 1, col = colchoice)
            })

#' Construct empirical margin
#'
#' @return An object of class \linkS4class{margin} signifying an empirical distribution function.
#' @export
#'
edf <- function(){
  new("margin", name = "edf")
}
