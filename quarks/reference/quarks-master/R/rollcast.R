#' Rolling one-step ahead forecasts of Value at Risk and Expected Shortfall
#'
#' Computes rolling one-step ahead forecasts of Value at Risk and Expected
#' Shortfall (Conditional Value at Risk) by means of plain historical
#' simulation age- and volatility-weighted historical simulation as well as
#' filtered historical simulation.
#'
#'
#' @param x a numeric vector of asset returns
#' @param p confidence level for VaR calculation; default is \code{0.975}
#' @param model model for estimating conditional volatility; options are \code{'EWMA'}
#' and \code{'GARCH'}; if \code{model = 'GARCH'}, additional arguments can be adjusted
#' via \code{...}; default is \code{'EWMA'}
#' @param method method to be used for calculation; default is \code{'plain'}
#' @param lambda decay factor for the calculation of weights; default is \code{0.98}
#' for \code{method = 'age'} and \code{0.94} for \code{method = 'vwhs'} or
#' \code{method = 'fhs'}
#' @param nout number of out-of-sample observations; most recent observations are used;
#' default is \code{NULL}
#' @param nwin window size for rolling one-step forecasting; most recent observations
#' before out-of-sample are used; default is \code{NULL}
#' @param nboot size of bootstrap sample; must be a single non-NA integer value
#' with \code{nboot > 0}; default is \code{NULL}
#' @param smoothscale a character object; defines the smoothing approach
#' for the unconditional variance from the logarithm of the squared centralized
#' returns; for \code{smoothscale = 'lpr'}, the unconditional
#' variance is smoothed via the \code{smoots::gsmooth()} function from the
#' \code{smoots} package; the bandwidth has to be chosen manually; otherwise the
#' default is used; if \code{smoothscale = 'auto'}, the function \code{smoots::msmooth()}
#' is employed and the bandwidth is chosen automatically (data-driven); see the
#' documentation of the \code{smoots} package for more information; is set to
#' \code{smoothscale = 'none'} by default
#' @param smoothopts additional arguments of \code{smoots::gsmooth()} and
#' \code{smoots::msmooth()}; see the documentation of the \code{smoots}
#' package for more information; is set to customized default settings
#' @param ... additional arguments of the \code{ugarchspec} function from the
#' \code{rugarch}-package; only applied if \code{model = 'GARCH'}; default
#' settings for the arguments \code{variance.model} and \code{mean.model} are:
#' \describe{
#' \item{\code{variance.model} = \code{list(model = 'sGARCH', garchOrder = c(1, 1))}}{}
#' \item{\code{mean.model} = \code{list(armaOrder = c(0, 0))}}{}
#' }
#'
#' @export
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{VaR}{Numerical vector containing out-of-sample forecasts of Value at
#' Risk}
#' \item{ES}{Numerical vector containing out-of-sample forecasts of Expected
#' Shortfall (Conditional Value at Risk)}
#' \item{xout}{Numerical vector containing out-of-sample returns}
#' \item{p}{Confidence level for VaR calculation}
#' \item{model}{Model for estimating conditional volatility}
#' \item{method}{Method to be used for calculation}
#' \item{nout}{Number of out-of-sample observations}
#' \item{nwin}{Window size for rolling one-step forecasting}
#' \item{nboot}{Size of bootstrap sample}
#' }
#' @examples
#'
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' n <- length(returns)
#' nout <- 250 # number of obs. for out-of-sample forecasting
#' nwin <- 1000 # window size for rolling forecasts
#'
#'
#' ### Example 1 - plain historical simulation
#' results1 <- rollcast(x = returns, p = 0.975, method = 'plain', nout = nout,
#'                      nwin = nwin)
#' matplot(1:nout, cbind(-results1$xout, results1$VaR, results1$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Plain HS - 97.5% VaR and ES for the DAX30 return series')
#'
#' ### Example 2 - age weighted historical simulation
#' results2 <- rollcast(x = returns, p = 0.975, method = 'age', nout = nout,
#'                      nwin = nwin)
#' matplot(1:nout, cbind(-results2$xout, results2$VaR, results2$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Age weighted HS - 97.5% VaR and ES for the DAX30 return series')
#'
#' ### Example 3 - volatility weighted historical simulation - EWMA
#' results3 <- rollcast(x = returns, p = 0.975, model = 'EWMA',
#'                      method = 'vwhs', nout = nout, nwin = nwin)
#' matplot(1:nout, cbind(-results3$xout, results3$VaR, results3$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Vol. weighted HS (EWMA) - 97.5% VaR and ES for the DAX30 return
#'   series')
#' \donttest{
#' ### Example 4 - volatility weighted historical simulation - GARCH
#' results4 <- rollcast(x = returns, p = 0.975, model = 'GARCH',
#'                      method = 'vwhs', nout = nout, nwin = nwin)
#' matplot(1:nout, cbind(-results4$xout, results4$VaR, results4$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Vol. weighted HS (GARCH) - 97.5% VaR and ES for the DAX30 return
#'   series')
#' }
#' ### Example 5 - filtered historical simulation - EWMA
#' results5 <- rollcast(x = returns, p = 0.975, model = 'EWMA',
#'                      method = 'fhs', nout = nout, nwin = nwin, nboot = 10000)
#' matplot(1:nout, cbind(-results5$xout, results5$VaR, results5$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Filtered HS (EWMA) - 97.5% VaR and ES for the DAX30 return
#'   series')
#' \donttest{
#' ### Example 6 - filtered historical simulation - GARCH
#' results6 <- rollcast(x = returns, p = 0.975, model = 'GARCH',
#'                      method = 'fhs', nout = nout, nwin = nwin, nboot = 10000)
#' matplot(1:nout, cbind(-results6$xout, results6$VaR, results6$ES),
#'   type = 'hll',
#'   xlab = 'number of out-of-sample obs.', ylab = 'losses, VaR and ES',
#'   main = 'Filtered HS (GARCH) - 97.5% VaR and ES for the DAX30 return
#'   series')
#' }

rollcast <- function(x, p = 0.975,
                     model = c("EWMA", "GARCH"),
                     method = c("plain", "age", "vwhs", "fhs"),
                     lambda = c(0.94, 0.98),
                     nout = NULL,
                     nwin = NULL,
                     nboot = NULL,
                     smoothscale = c("none", "lpr", "auto"),
                     smoothopts = list(),
                     ...) {

    if (length(x) <= 1 || any(is.na(x)) || !is.numeric(x)) {
        stop("A numeric vector of length > 1 and without NAs must be passed to",
             " 'x'.")
    }
    if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0)) {
        stop("The argument 'p' must be a single non-NA double value with ",
             "0 < p < l.")
    }
    if (!(length(method) %in% c(1, 4)) || any(is.na(method)) ||
        !is.character(method) || !all(method %in% c("plain", "age", "vwhs",
                                                    "fhs")))
        {
        stop("A single character value must be passed to 'method'.",
             "Valid choices are 'plain', 'age', 'vwhs' or 'fhs'.")
    }
    if (!(length(lambda) %in% c(1, 2)) || any(is.na(lambda)) ||
        !is.numeric(lambda) || all(lambda < 0) || all(lambda >= 1)) {
        stop("The argument 'lambda' must be a single non-NA double value with ",
             "0 < lambda < 1.")
    }
    if (length(nout) != 1 || is.na(nout) || !is.numeric(nout) || nout < 0 ||
        is.null(nout)) {
        stop("The argument 'nout' must be a single non-NA integer value.")
    }
    if (length(nwin) != 1 || is.na(nwin) || !is.numeric(nwin) || nwin <= 1 ||
        is.null(nwin)) {
        stop("The argument 'nwin' must be a single non-NA integer value
             with nwin > 1.")
    }
    if (nwin > (length(x) - nout) || nout > (length(x) - nwin) ||
       (nwin + nout) > length(x)) {
        stop("Window size and (or) out-of-sample size are too large.")
    }
    if (all(method == c("plain", "age", "vwhs", "fhs"))) {
      method <- "plain"
    }
    if (method == 'fhs' && (length(nboot) != 1 || is.na(nboot) ||
        !is.numeric(nboot) || nboot <= 0)) {
        stop("The argument 'nboot' must be a single non-NA integer value with ",
             "nboot > 0.")
    }
    if (!(length(model) %in% c(1, 2)) || any(is.na(model)) ||
        !is.character(model) || !all(model %in% c("EWMA", "GARCH"))) {
        stop("A single character value must be passed to 'model'. ",
              "Valid choices are 'EWMA' or 'GARCH'.")
    }
    if (!(length(smoothscale) %in% c(1, 3)) || any(is.na(smoothscale)) ||
        !is.character(smoothscale) || !all(smoothscale %in% c("none", "lpr", "auto"))) {
        stop("A single character value must be passed to 'smoothscale'. ",
             "Valid choices are 'none', 'lpr' or 'auto'.")
    }

    if (method != "fhs")
        nboot = "N/A"
    if (all(lambda == c(0.94, 0.98)) && method == "age")
        lambda <- 0.98
    if (all(lambda == c(0.94, 0.98)) && method == "vwhs" || method == "fhs")
        lambda <- 0.94
    if (all(model == c("EWMA", "GARCH")))
        model <- "EWMA"
    if (all(smoothscale == c("none", "lpr", "auto")))
        smoothscale <- "none"
    if (smoothscale == "auto") {
        if (is.null(smoothopts[["alg"]])) smoothopts[["alg"]] <- "A"
        if (is.null(smoothopts[["p"]])) smoothopts[["p"]] <- 3
        message("\n", "Please beware that smoothing of the unconditional ",
                "scale function is still in an experimental stage!", "\n")
    }
    if (smoothscale == "lpr") {
        if (is.null(smoothopts[["b"]])) smoothopts[["b"]] <- 0.15
        message("\n", "Please beware that smoothing of the unconditional ",
                "scale function is still in an experimental stage!", "\n")
    }
    # defining progress bar
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsed || Estimated time remaining: :eta]",
                           total = nout,
                           complete = "=",
                           incomplete = "-",
                           current = ">",
                           show_after = 0.2,
                           clear = TRUE)

    ### Main code
    n <- length(x)
    nin <- n - nout
    xin <- x[1:nin]
    if (nout == 0){
        xout <- NA
    }
    else{
        xout <- x[(nin + 1):n]
    }

    xstart <- xcast <- xin[(nin - nwin + 1):nin]
    fcasts <- matrix(NA, max(nout, 1), 2)
    sfc <- 1

    switch(method,
           plain = {
               fcastfun <- hs
               funargs <- list(x = xcast, p = p, method = method)
               },
           age = {
               fcastfun <- hs
               funargs <- list(x = xcast, p = p, method = method, lambda = lambda)
                },
           vwhs = {
               fcastfun <- vwhs
               funargs <- list(x = xcast, p = p, lambda = lambda,
                               model = model, ...)
                },
           fhs = {
               fcastfun <- fhs
               funargs <- list(x = xcast, p = p,
                               lambda = lambda, nboot = nboot,
                               model = model, ...)
                }
           )

    if (smoothscale %in% c("lpr", "auto")) {
        npest <- smooth.help(x = xstart, smoothscale = smoothscale,
                             smoothopts = smoothopts)
        xcast <- npest[["xstd"]]
        sfc <- npest[["sfc"]]
    }

    fcasts[1, ] <- as.double(unlist(do.call(what = fcastfun, args = funargs),
                          use.names = FALSE)[1:2]) * sfc
    if (nout > 1) {
        for (i in 2:nout) {
            pb$tick() # progress bar
            if (i <= nwin) {
                xcast <- c(xstart[i:nwin], xout[1:(i - 1)])
                funargs[["x"]] = xcast
                if (smoothscale %in% c("lpr", "auto")) {
                    npest <- smooth.help(x = xcast, smoothscale = smoothscale,
                                         smoothopts = smoothopts)
                    xcast <- npest[["xstd"]]
                    sfc <- npest[["sfc"]]
                }
                fcasts[i, ] <- as.double(unlist(do.call(what = fcastfun, args = funargs),
                                      use.names = FALSE)[1:2]) * sfc
            }
            if (i > nwin) {
                xcast <- xout[(i - nwin):(i - 1)]
                funargs[["x"]] = xcast
                if (smoothscale %in% c("lpr", "auto")) {
                    npest <- smooth.help(x = xcast, smoothscale = smoothscale,
                                         smoothopts = smoothopts)
                    xcast <- npest[["xstd"]]
                    sfc <- npest[["sfc"]]
                }
                fcasts[i, ] <- as.double(unlist(do.call(what = fcastfun, args = funargs),
                                      use.names = FALSE)[1:2]) * sfc
            }
        }
    }

    pb$terminate(); invisible() # terminate progress bar


    VaR <- fcasts[, 1]
    ES <- fcasts[, 2]

    if (model == "GARCH" && method %in% c("vwhs", "fhs"))
        model <- list(...)$variance.model$model
    if (is.null(model))
        model <- "sGARCH"
    if (method %in% c("plain", "age"))
        model <- "EWMA"

    results <- list(VaR = VaR, ES = ES, xout = xout, p = p, model = model,
                    method = method, nout = nout, nwin = nwin, nboot = nboot)

    message("\n", "Calculations completed.", "\n")

    class(results) <- "quarks"
    attr(results, "function") <- "rollcast"
    results
}
