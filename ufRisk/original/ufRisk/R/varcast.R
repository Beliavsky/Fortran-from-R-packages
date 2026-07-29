#' Calculation of one-step ahead forecasts of Value at Risk and Expected Shortfall (parametric and
#' semiparametric)
#'
#'One-step ahead forecasts of Value at Risk and Expected Shortfall for a selection of short-memory
#'and long-memory parametric as well as semiparametric GARCH-type models are
#'computed.
#'
#' @importFrom stats 'dt' 'qt' 'sd' 'arima'
#' @importFrom utils 'tail'
#'
#' @export
#'
#' @param x a vector containing the price series.
#' @param a.v confidence level for calculating VaR; is set to \code{0.99} by default.
#' @param a.e confidence level for calculating ES; is set to \code{0.975} by default.
#' @param model model to be estimated. Options are 'sGARCH', 'eGARCH', 'apARCH',
#' 'lGARCH', 'fiGARCH' and 'filGARCH'; is set to \code{'sGARCH'} by default.
#' @param garchOrder orders to be estimated; c(1, 1), i.e. p = q = 1, is the
#' default.
#' @param distr distribution to use for the innovations of the respective GARCH model;
#' is set to \code{'std'} by default
#' @param n.out size of out-sample; is set to \code{250} by default.
#' @param smooth a character object; defines the data-driven smoothing approach
#' for the estimation of the nonparametric scale function; for
#' \code{smooth = 'lpr'}, the scale function is obtained from the logarithm
#' of the squared centralized returns by means of the \code{msmooth()} function or
#' \code{tsmoothlm()} function if \code{model} is set to \code{'sGARCH'}, \code{'eGARCH'},
#' \code{'apARCH'} and \code{lGARCH'} or \code{'fiGARCH} and \code{'filGARCH'},
#' respectively; is set to \code{smooth = 'none'} by default.
#' @param ... depending on the choice of \code{model}, further arguments can be
#' passed to either \code{smoots::msmooth()} or to \code{tsmoothlm()};
#' if no other arguments are given, the default settings
#' are used for both functions with the exception of \code{p = 3}.
#'
#' @details
#' Let \eqn{Y_t} be a (demeaned) return series. The semiparametric extension of
#' the GARCH(p,q) model (Bollerslev, 1986) is called a Semi-GARCH model
#' (Feng, 2004) and is defined by
#'
#' \deqn{Y_t = s(x_t)\sigma_t \eta_t,}
#'  with \eqn{\eta_t \sim IID(0,1)} and
#'
#' \deqn{\sigma^2_t = \alpha_0 + \sum_{i=1}^p \alpha_i Y^2_{t-i}
#' + \sum_{j=1}^q \beta_j \sigma^2_{t-j},}
#'
#' where \eqn{\sigma_t > 0} are the conditional standard deviations, \eqn{s(x_t) > 0} is
#' a nonparametric scale function with \eqn{x_t} being the rescaled observation
#' time points on the interval [0,1] and \eqn{\alpha_i}
#' and \eqn{\beta_j} are non-negative real valued coefficients, except for
#' \eqn{\alpha_0}, which must satisfy \eqn{\alpha_0 > 0}. Furthermore, it is assumed that
#' Var\eqn{(\sigma_t \eta_t) = 1}. In this functions, different
#' short-memory and long-memory GARCH-type models are selectable for the parametric part of the
#' model. More specifically, the standard GARCH (Bollerslev, 1986), the
#' Log-GARCH (Pantula, 1986; Geweke, 1986; Milhoj, 1988), the eGARCH
#' (Nelson, 1991), the APARCH (Ding et al., 1993), the FIGARCH (Baillie et al., 1996)
#' and the FI-Log-GARCH (Feng et al., 2020) model are implemented. For more
#' information on the representations of the last three models mentioned, we
#' refer the reader to the corresponding references listed in the references
#' section.
#'
#' While the innovations \eqn{\eta_t} must be i.i.d. (independent and identically
#' distributed) with zero-mean and unit-variance and while any
#' distribution that satisfies these conditions is suitable, the standardized
#' t-distribution is selected for the estimation of the models and computation
#' of the Value at Risk (VaR) as well as the Expected Shortfall (ES) within
#' this function.
#'
#' For a given level \eqn{\alpha \in (0, 1)},
#'
#' \deqn{VaR(\alpha) = inf \{z \in R: F_L(z) \geq \alpha\}}
#'
#' defines the VaR at level alpha. In this definition, \eqn{L} is the
#' loss variable (making a loss is denoted as a positive value, whereas gains
#' are negative values) and \eqn{F_L} is its cumulative distribution function.
#' Explained differently, \eqn{VaR(\alpha)} is the \eqn{\alpha}-quantile of the loss
#' distribution.
#'
#' The ES for a level \eqn{\alpha}, however, is given by
#'
#' \deqn{ES(\alpha) = (1 / (1 - \alpha)) \int_{\alpha}^1 VaR(u)du,}
#'
#' i.e. it is the expected loss in case \eqn{VaR(\alpha)} is exceeded. More
#' information on these risk measures can be found on pp. 64-72 in McNeil et
#' al. (2015).
#'
#' To apply the function, a numeric vector that contains the price series that
#' is to be analyzed ordered from past to present must be passed to the
#' argument \code{x}. Furthermore, the user can set different levels of alpha
#' for the VaR and the ES via the arguments \code{a.v} and \code{a.e},
#' respectively. A parametric short-memory or long-memory GARCH-type model can be
#' selected by means of \code{model}, which only accepts a single-element character vector
#' as input. At the time of the release of package version 1.0.0, a standard
#' GARCH ('sGARCH'), a Log-GARCH ('lGARCH'), an eGARCH ('eGARCH'), an APARCH
#' ('apARCH'), a FIGARCH ('fiGARCH') and a FI-Log-GARCH ('filGARCH') model can be selected,
#' each with conditional t-distribution. By default, a standard GARCH model is applied.
#' The orders of the GARCH-type models can be defined with \code{garchOrder},
#' which is a numeric vector with two elements. Its first element is the
#' ARCH order p, whereas the GARCH order q can be adjusted via the second
#' element. If no adjustments are made, the orders p = q = 1 are selected. The
#' number of out-sample observations is set via the argument \code{n.out}. If n
#' is the total number of observations of the
#' whole price series, the model is estimated for the first n - n.out
#' observations (in-sample), while the VaR and the ES are obtained for the last
#' n.out observations (out-sample) based on the estimated model for the
#' in-sample. Moreover, the data-driven estimation method of the underlying
#' scale function can be adjusted via the argument \code{smooth}. If
#' \code{smooth = 'lpr'} is selected, the scale function is obtained by
#' applying an iterative plug-in algorithm logarithm of the squared
#' centralized returns. Depending on the setting of \code{model} an algorithm
#' proposed by Feng, Gries and Fritz (2020) or by Letmathe, Feng and Uhde
#' (2021) is employed. In the former case, the function \code{msmooth()} of the
#' \code{smoots} package is applied and for the latter the \code{tsmoothlm()}
#' function of the \code{esemifar} package is used. An ellipsis \code{...} is
#' implemented to allow for additional arguments for \code{msmooth()} and
#' \code{tsmoothlm()}.
#'
#' NOTE:
#'
#' This function makes use of the \code{arima()} function of the stats package,
#' of the \code{fracdiff()} function of the \code{fracdiff} package, of the
#' \code{ugarchspec()} and \code{ugarchfit()} functions of the \code{rugarch}
#' package, of the \code{msmooth()} function of the \code{smoots} package
#' and of the \code{esemifar()} function of the \code{esemifar} for estimation.
#' Moreover, Log-GARCH and FI-Log-GARCH models in the parametric part of the
#' complete models are estimated via their ARMA and FARIMA representations,
#' respectively, and must therefore satisfy \eqn{p \geq q}.
#'
#' @return
#' This function returns a list with the following elements.
#' \describe{
#' \item{model}{selected model for estimation}
#' \item{mean}{the estimated mean of the in-sample returns}
#' \item{model.fit}{estimated model parameters for the parametric part of the
#' in-sample}
#' \item{np.est}{the estimation results for the nonparametric part of the
#' in-sample model}
#' \item{ret.in}{in-sample return series}
#' \item{ret.out}{out-sample return series}
#' \item{sig.in}{estimated in-sample total volatility}
#' \item{sig.fc}{out-sample forecasts of the total volatility}
#' \item{scale}{the estimated nonparametric scale function values for the
#' in-sample}
#' \item{scale.fc}{the scale function forecast for the out-sample}
#' \item{VaR.e}{out-sample forecasts of the (1-\code{a.e})100\% VaR}
#' \item{VaR.v}{out-sample forecasts of the (1-\code{a.v})100\% VaR}
#' \item{ES}{out-sample forecasts of the (1-\code{a.e})100\% ES}
#' \item{dfree}{estimated degrees of freedom for the standardized returns}
#' \item{a.v}{coverage level for the 99 \% VaR}
#' \item{a.e}{coverage level for 97.5 \% VaR}
#' \item{garchOrder}{the orders p and q of the implemented GARCH-type model}
#' }
#'
#' @author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University) \cr
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}
#'
#' @references
#' Baillie, R. T., Bollerslev, T., & Mikkelsen, H. O. (1996). Fractionally
#' integrated generalized autoregressive conditional heteroskedasticity.
#' In: Journal of Econometrics, 74.1, pp. 3-30.
#'
#' Bollerslev, T. (1986) Generalized autoregressive conditional
#' heteroskedasticity. In: Journal of Econometrics 31.3, pp. 307-327.
#'
#' Ding, Z., Granger, C.W., and Engle, R.F. (1993). A long memory property
#' of stock market returns and a new model. In: Journal of Empirical Finance
#' 1.1, pp. 83-106.
#'
#' Feng, Y. (2004). Simultaneously modeling conditional heteroskedasticity and
#' scale change. In: Econometric Theory 20.3, pp. 563-596.
#'
#' Feng, Y., Beran, J., Letmathe, S., & Ghosh, S. (2020). Fractionally
#' integrated Log-GARCH with application to value at risk and expected
#' shortfall (No. 137). Paderborn University, CIE Center for
#' International Economics.
#'
#' Pantula, S.G. (1986). Modeling the persistence of conditional variances:
#' a comment. In: Econometric Reviews 5, pp. 79-97.
#'
#' Geweke, J. (1986). Comment on: Modelling the persistence of conditional
#' variances. In: Econometric Reviews 5, pp. 57-61.
#'
#' Letmathe, S., Feng, Y., & Uhde, A. (2021). Semiparametric GARCH models with
#' long memory applied to Value at Risk and Expected Shortfall (No. 141).
#' Paderborn University, CIE Center for International Economics.
#'
#' McNeil, A.J., Frey, R., and Embrechts, P. (2015). Quantitative risk
#' management: concepts, techniques and tools - revised edition. Princeton
#' University Press.
#'
#' Milhoj, A. (1988). A Multiplicative parametrization of ARCH models.
#' Universitetets Statistiske Institut.
#'
#' Nelson, D. B. (1991). Conditional heteroskedasticity in asset returns: A
#' new approach. In: Econometrica: Journal of the Econometric Society, 347-370.
#'
#'
#' @examples
#'
#'# Example for Walmart Inc. (WMT)
#' prices = WMT$price.close
#'
#' # forecasting VaR and ES
#' results = varcast(prices, model = 'sGARCH', n.out = 250)
#' ret.out = results$ret.out
#' n.out = length(ret.out)
#' VaR97.5 = results$VaR.e
#' VaR99 = results$VaR.v
#' ES = results$ES
#'
#'# plotting VaR at 99% coverage
#' matplot(1:n.out, cbind(-ret.out, VaR99),
#'   type = 'hl',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = '99% VaR (red) for the WMT return series')
#'
#'# plotting VaR at 97.5% coverage and corresponding ES
#' matplot(1:n.out, cbind(-ret.out, ES, VaR97.5),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = '97.5% VaR (green) and ES (red) for the WMT return series')
#'

varcast <- function(x, a.v = 0.99, a.e = 0.975,
                    model = c("sGARCH", "lGARCH", "eGARCH", "apARCH", "fiGARCH",
                              "filGARCH"),
                    garchOrder = c(1, 1),
                    distr = c("norm", "std"),
                    n.out = 250,
                    smooth = c("none", "lpr"), ...) {

  if (length(x) <= 1 || !all(!is.na(x)) || !is.numeric(x)) {
    stop("A numeric vector of length > 1 and without NAs must be passed to",
         " 'x'.")
  }
  if (length(a.v) != 1 || is.na(a.v) || !is.numeric(a.v) || a.v <= 0 || a.v >= 1) {
    stop("A single numeric value that satisfies >0 and <1 must be passed to",
         " 'a.v'")
  }
  if (length(a.e) != 1 || is.na(a.e) || !is.numeric(a.e) || a.e <= 0 || a.e >= 1) {
    stop("A single numeric value that satisfies >0 and <1 must be passed to",
         " 'a.e'")
  }
  if (!(length(model) %in% c(1, 6)) || !all(!is.na(model)) || !is.character(model)) {
    stop("A single character value must be passed to 'model'.")
  }
  if (all(model == c("sGARCH", "lGARCH", "eGARCH", "apARCH", "fiGARCH", "filGARCH"))) {
    model <- "sGARCH"
  }
  if (length(garchOrder) != 2 || !all(!is.na(garchOrder)) || !is.numeric(garchOrder) ||
      (garchOrder[[1]] == 0 && garchOrder[[2]] == 0)) {
    stop("A vector of length 2 must be passed to 'garchOrder' giving the",
         " first and second order of the GARCH-type model.")
  }
  garchOrder <- floor(garchOrder)
  if (length(n.out) != 1 || is.na(n.out) || !is.numeric(n.out) || n.out <= 0) {
    stop("A single positive integer must be passed to 'n.out'.")
  }
  n.out <- floor(n.out)

  if (!(length(distr) %in% c(1, 2)) || !all(!is.na(distr)) ||
      !is.character(distr)) {
    stop("Input argument 'distr' must be a single character ",
         "value. Valid choices are 'norm' and 'std'.")
  }

  if (all(distr == c("norm", "std"))) distr <- "std"
  if (length(distr) != 1 || !(distr) %in% c("norm", "std")) {
    stop("Input argument 'distr' must be a single character ",
         "value. Valid choices are 'norm' and 'std'.")
  }

  if (!(length(smooth) %in% c(1, 2)) || !all(!is.na(smooth)) ||
      !is.character(smooth)) {
    stop("Input argument 'smooth' must be a single character ",
         "value. Valid choices are 'none' and 'lpr'.")
  }

  if (all(smooth == c("none", "lpr"))) smooth <- "none"
  if (length(smooth) != 1 || !(smooth) %in% c("none", "lpr")) {
    stop("Input argument 'smooth' must be a single character ",
         "value. Valid choices are 'none' and 'lpr'.")
  }

  if (smooth == "lpr" && model %in% c("sGARCH", "lGARCH", "eGARCH", "apARCH")) {
    smooth.method = "smoots"
  }

  if (smooth == "lpr" && model %in% c("fiGARCH", "filGARCH")) {
    smooth.method = "esemifar"
  }

  dots <- list(...)
  ret <- diff(log(x))
  n.ret <- length(ret)
  n.in <- n.ret - n.out
  ret.in.t <- ret[1:n.in]
  mean.ret.in <- mean(ret.in.t)
  ret.in <- ret.in.t - mean.ret.in
  ret.out.t <- ret[(n.in + 1):n.ret]
  ret.out <- ret.out.t - mean.ret.in
  p.true <- p <- garchOrder[1]
  q.true <- q <- garchOrder[2]
  if (p.true == 0) p <- 1
  if (q.true == 0) q <- 1
  l <- max(p, q)

  if (is.null(dots[["p"]])) {
    dots[["p"]] <- 3
  }
  dots[["y"]] <- log(ret.in^2)

  switch(smooth,
         none = {
           np.est <- NA
           zeta.in <- ret.in
           zeta.out <- ret.out
           res.in <- log(ret.in^2)
           res.out <- log(ret.out^2)
           sxt <- 1
           sfc <- 1
           mule <- mean(res.in)
           Csig <- 1
         },
         lpr = {
           switch(smooth.method,
                  smoots = {
                    if (is.null(dots[["alg"]])) {
                      dots[["alg"]] <- "A"
                    }
                    np.est <- suppressMessages(do.call(what = smoots::msmooth,
                                                       args = dots))
                    res.in <- np.est[["res"]]
                    res.out <- log(ret.out^2) - np.est[["ye"]][n.in]
                    mule <- -log(mean(exp(res.in)))
                    sxt <- exp(0.5 * (np.est[["ye"]] - mule))
                    sfc <- sxt[n.in]
                    zeta.in <- ret.in / sxt
                    zeta.out <- ret.out / sfc
                  },
                  esemifar = {
                    np.est <- suppressMessages(do.call(what = esemifar::tsmoothlm,
                                                       args = dots))
                    res.in <- np.est[["res"]]
                    res.out <- log(ret.out^2) - np.est[["ye"]][n.in]
                    Csig <- sd(ret.in / exp(0.5 * np.est[["ye"]]))
                    sxt <- exp(0.5 * np.est[["ye"]]) * Csig
                    sfc <- sxt[n.in]
                    zeta.in <- ret.in / sxt
                    zeta.out <- ret.out / sfc
                  })
         }
  )

  zeta.n <- tail(zeta.in, n = l)
  zeta.fc <- c(zeta.n, zeta.out[1:n.out])

  if (model %in% c("sGARCH", "eGARCH", "apARCH", "fiGARCH")) {
    spec <- rugarch::ugarchspec(variance.model =
                                  list(model = model,
                                       garchOrder = c(p.true, q.true)),
                                mean.model = list(armaOrder = c(0, 0),
                                                  include.mean = FALSE),
                                distribution.model = distr)
    model_fit <- rugarch::ugarchfit(spec = spec, data = c(zeta.in, zeta.out),
                                    out.sample = n.out,
                                    fit.control = list(trunclag = min(n.in, 1000)))

    if (distr == "std") dfree <- model_fit@fit$coef[["shape"]]
    if (distr == "norm") dfree <- NA

    c.sig <- as.numeric(rugarch::sigma(model_fit))
    sig.in <- c.sig * sxt
    sig.fc <- rugarch::ugarchforecast(model_fit, n.ahead = 1, n.roll = n.out - 1)
    sig.fc <- as.numeric(rugarch::sigma(sig.fc)) * sfc
  }


  switch(model,
         lGARCH = {
           if (p.true < q.true) {
             stop("p >= q must be satisfied for the estimation of a Log-GARCH
                  model ", "via its ARMA representation.")
           }
           model_fit <- arima(res.in,
                              order = c(p.true, 0, q.true),
                              include.mean = FALSE)
           y.cent <- c(res.in, res.out) - mule * isTRUE(smooth == "none")
           mulz <- -log(mean(exp(model_fit[["residuals"]])))
           ar <- model_fit$model$phi
           ma <- -model_fit$model$theta
           k <- 50
           d <- 0
           coef.all <- arfilt(ar, ma, d, k)

           pre0 <- c()
           add0 <- 0
           if (n.in < k) {
             add0 <- k - n.in
             pre0 <- rep(0, add0)
           }
           y.fc.out <- (1:n.out) * 0
           y.cent <- c(pre0, y.cent)

           for (i in (add0 + n.in + 1):(add0 + n.ret)) {
             y.fc.out[i - (add0 + n.in)] <- coef.all %*%
               y.cent[(i - 1):(i - k)]
           }

           sig.fc <- exp(0.5 * (y.fc.out + mule - mulz))
           sig.fc.in <- exp(0.5 * (y.cent[1:n.in] - model_fit[["residuals"]]
                                   + mule - mulz))
           sig.fc <- sig.fc * sfc
           sig.in <- sig.fc.in * sxt
           ret.sd <- zeta.in / sig.in
           if (distr == "std") {
             dfree <- as.numeric(rugarch::fitdist("std", ret.sd)$pars[[3]])
           }
           if (distr == "norm") {
             dfree <- NA
           }
         },
         filGARCH = {
           if (p.true < q.true) {
             stop("p >= q must be satisfied for the estimation of a Log-GARCH
                  model ", "via its ARMA representation.")
           }
           model_fit <-  suppressWarnings(fracdiff::fracdiff(res.in,
                                                             nar = p.true,
                                                             nma = q.true))
           ar <- model_fit[["ar"]]
           ma <- model_fit[["ma"]]
           d <- model_fit[["d"]]
           k <- 50
           coef.all <- arfilt(ar, ma, d, k)

           pre0 <- c()
           add0 <- 0
           y.cent <- c(res.in, res.out)
           y.cent <- y.cent - mean(res.in)
           if (n.in < k) {
             add0 <- k - n.in
             pre0 <- rep(0, add0)
           }
           y.fc.out <- (1:n.out) * 0
           y.cent <- c(pre0, y.cent)

           for (i in (add0 + n.in + 1):(add0 + n.ret)) {
             y.fc.out[i - (add0 + n.in)] <- coef.all %*% y.cent[(i - 1):(i - k)]
           }
           sig.fc <- exp(0.5 * (y.fc.out + mean(res.in)))
           y.fc.in <- model_fit[["fitted"]]
           sig.fc.in <- exp(0.5 * y.fc.in)
           ret.sd <- zeta.in / (sig.fc.in * sxt)
           if (distr == "std") {
             dfree <- as.numeric(rugarch::fitdist("std", ret.sd)$pars[3])
           }
           if (distr == "norm") {
             dfree <- NA
           }
           if (smooth == "lpr") {
             Ch <- sd(zeta.in / sig.fc.in)
           }
           if (smooth == "none") {
             Ch <- sd(ret.sd)
           }
           sig.fc <- sig.fc * sfc * Ch
           sig.in <- sig.fc.in * sxt * Ch
         }
  )

  switch(distr,
         norm = {
           sdev <- 1
           ES0t <- 1
           arglst <- list(p = c(a.e, a.v), distr = distr)
         },
         std = {
           sdev <- sqrt(dfree / (dfree - 2))
           ES0t <- (dfree + quant(a.e, df = dfree, distr = distr)^2) / (dfree - 1)
           arglst <- list(p = c(a.e, a.v), df = dfree, distr = distr)
         }
  )

  quants <- do.call(what = quant, args = arglst)
  q.e <- quants[1]
  q.a <- quants[2]
  VaR.fc <- -mean.ret.in + sig.fc * q.a / sdev
  VaR.e.fc <- -mean.ret.in + sig.fc * q.e / sdev

  arglst[[1]] <- q.e; names(arglst)[1] <- "x"
  d.es <- do.call(what = dens, args = arglst)
  ES0 <- d.es / (1 - a.e) * ES0t / sdev
  Esfc <- -mean.ret.in + sig.fc * ES0

  if (smooth == "lpr") {
    model <- paste0("Semi-", model)
  }
  results <- list(model = model, mean = mean.ret.in,
                  model.fit = model_fit, ret.in = ret.in.t,
                  ret.out = ret.out.t, sig.in = sig.in,
                  sig.fc = sig.fc, scale = sxt, scale.fc = sfc, VaR.e = VaR.e.fc,
                  VaR.v = VaR.fc, ES = Esfc, dfree = dfree, a.v = 1 - a.v,
                  a.e = 1 - a.e, garchOrder = garchOrder, np.est = np.est)

  class(results) <- "ufRisk"
  attr(results, "function") <- "varcast"
  results
}
