### ARMA-GARCH #################################################################

##' @title Fitting ARMA-GARCH processes
##' @param x matrix-like data structure
##' @param ugarchspec.list object of class uGARCHspec (as returned by
##'        ugarchspec()) or a list of such
##' @param solver string indicating the solver used; see ?ugarchfit
##' @param verbose logical indicating whether verbose output is given
##' @param ... additional arguments passed to the underlying ugarchfit() or HoltWinters()
##' @return list of length equal to the number of columns of x containing
##'         the fitted objects
##' @author Marius Hofert
##' @note The default ugarchspec.list fits an ARMA(1,1)-GARCH(1,1) with N(0,1)
##'       standardized residuals
fit_ARMA_GARCH <- function(x, ugarchspec.list = ugarchspec(), solver = "hybrid",
                           verbose = FALSE, ...)
{
    ## Checking and expanding ugarchspec.list to a list of length d
    if(!is.matrix(x)) x <- cbind(x) # is.matrix() is also true for 'xts' objects
    d <- ncol(x)
    isGARCHspec <- function(spec) is(spec, class2 = "uGARCHspec")
    if(isGARCHspec(ugarchspec.list)) ugarchspec.list <- rep(list(ugarchspec.list), d) # repeat
    if(is.list(ugarchspec.list))
        stopifnot(length(ugarchspec.list) == d && all(sapply(ugarchspec.list, isGARCHspec)))
    stopifnot(is.logical(verbose))

    ## Iterate over all time series
    fit  <- vector("list", length = d)
    warn <- vector("list", length = d)
    err  <- vector("list", length = d)
    ## if(verbose) {
    ##     pb <- txtProgressBar(max = d, style = if(isatty(stdout())) 3 else 1)
    ##     on.exit(close(pb)) # on exit, close progress bar
    ## }
    for(j in seq_len(d)) {
        res <- catch(ugarchfit(ugarchspec.list[[j]], data = x[,j], solver = solver,
                               ...)) # fitting
        if(!is.null(res$value)) fit[[j]] <- res$value
        if(!is.null(res$warning)) warn[[j]] <- res$warning
        if(!is.null(res$error)) err[[j]]  <- res$error
        ## Progress
        if(verbose) {
            ## setTxtProgressBar(pb, j) # update progress bar
            if(j %% ceiling(d/20) == 0)
                cat(sprintf("%2d%% ", ceiling(j/d * 100)))
            if(j == d) cat("\n")
        }
    }

    ## Return
    list(fit = fit, warning = warn, error = err)
}


### GARCH(1,1) #################################################################

## Note:
## 1) The process is given by:
##    X_t = \sigma_t * Z_t,
##    \sigma_t^2 = \alpha_0 + \alpha_1 * X_{t-1}^2 + \beta_1 * \sigma_{t-1}^2
## 2) Zumbach (2000, "The pitfalls in fitting GARCH (1,1) processes"), see
##    Advances in Quantitative Asset Management, 179--200, suggested to use
##    a method-of-moment estimator for the variance parameter and to maximize
##    a reparameterized log-likelihood for better numerical robustness in
##    fitting the GARCH(1,1) process.

## Idea reparameterization:
## - Let \mu_{cor} = \alpha_1 + \beta_1, \mu_{ema} = \beta_1/\mu_{cor}
##   and \sigma^2 = \alpha_0/(1-\mu_{cor}).
## - Basic reparameterization (see (3), partly):
##   \sigma_t^2 = \alpha_0 + \alpha_1 * X_{t-1}^2 + \beta_1 * \sigma_{t-1}^2; see McNeil et al. (2015, p. 124)
##              = \alpha_0/(1-\mu_{cor}) (1-\mu_{cor})
##                + \mu_{cor} ((\beta_1/\mu_{cor})\sigma_{t-1}^2 + (1-\mu_{ema}) X_{t-1}^2)
##              = \sigma^2 * (1-\mu_{cor}) + \mu_{cor} * (1-\mu_{ema}) * X_{t-1}^2
##                + \mu_{cor} * \mu_{ema} * \sigma_{t-1}^2
##              = <nonrecursive part> + \mu_{cor} * \mu_{ema} * \sigma_{t-1}^2
## - Time interval: delta = 1 (= 1 day = daily data; for yearly use delta = 250)
##   => Annualized: \sigma^2_{ann} = (250/delta) * \sigma^2_{daily}
## - Further reparameterization:
##   + \mu_{cor} = \exp(-\delta/\tau_{cor})
##   + \mu_{ema} = \exp(-\delta/\tau_{ema})
##   + z_{cor} = log(\tau_{cor})
##   + z_{ema} = log(\tau_{ema})
## - Parameters finally used: \sigma_{ann} (> 0), z_{cor}, z_{ema} (both unconstrained)

##' @title Reparameterized Log-likelihood of a GARCH(1,1) as a Function of
##'        Two Parameters
##' @param param numeric(2) providing the two parameters z_{cor}, z_{ema}; see (13)
##' @param sig2 numeric(1), > 0, third parameter providing the annualized variance
##' @param x data (e.g., risk-factor changes)
##' @param x2.tm1 squared data x_{t-1}^2 shifted back one unit in time
##' @param delta unit of time (defaults to 1 = daily data; for yearly data, use 250)
##' @param distr standardized innovation distribution (distribution of Z; E(Z) = 0, Var(Z) = 1);
##'        available are N(0,1) and standardized t (t_nu(0, 1))
##' @return log-likelihood (computed based on the given param and sig2)
##' @author Marius Hofert
##' @note The likelihood in Zumbach (2000) only addresses distr = "norm" and has
##'       a wrong factor of 1/n in it.
logLik_GARCH_11 <- function(param, # (z_{cor}, z_{ema}) parameters (two out of three)
                            sig2 = mean(x^2), # third parameter
                            x, x2.tm1 = c(mean(x^2), head(x, n = -1)^2), # data; X_0^2 (sample mean of squared data), ..., X_{t-1}^2
                            delta = 1, distr = c("norm", "st"))
{
    ## Compute standardized residuals based on given parameters
    mu <- exp(-delta / exp(param)) # exp(-delta/tau_.) for tau_. = exp(z_.); see Formulas (12) and (13)
    nonrec <- sig2 * (1-mu[1]) + mu[1] * (1-mu[2]) * x2.tm1 # nonrecursive part (see above)
    sig2.t <- as.numeric(filter(nonrec, filter = prod(mu), "recursive", init = sig2)) # filter out \sigma_1^2, ..., \sigma_t^2 and multiply with mu[1] * mu[2] = \mu_{cor} * \mu_{ema}
    sig.t <- sqrt(sig2.t)
    Z.t <- x / sig.t
    ## Log-likelihood in the remaining parameters based on X_t
    ## (but expressed in terms of the (standardized) innovation distribution)
    switch(match.arg(distr),
    "norm" = { # numerical robustness shown in Zumbach (2000)
        stopifnot(length(param) == 2)
        ## Note: X_t ~ N(0, sig2.t), so the log-likelihood contribution is
        ##       log((1/\sqrt{2\pi}\sigma_t) \exp(-(X_t/\sigma_t)^2/2)), or
        ##       log((1/\sqrt{2\pi}) \exp(-Z_t^2/2)) - log(\sigma_t)
        ##       => can work with dnorm(, mean = 0, sd = 1) for the first part
        sum(dnorm(Z.t, log = TRUE) - log(sig.t))
    },
    "st" = { # numerical robustness unclear (assumed to be better than without reparameterization, as for N(0,1))
        stopifnot(length(param) == 3) # last one is the degrees of freedom
        df <- param[3]
        if(df > 2) {
            ## Note:
            ## - X ~ t_nu (Var(X) = nu/(nu-2)) => Y = \sqrt{(nu-2)/nu} X ~ t_nu(0, (nu-2)/nu) with Var(Y) = 1
            ##   Since F_Y(x) = F_X(x/sqrt{(nu-2)/nu}), f_Y(x) = f_X(x/sqrt{(nu-2)/nu}) / sqrt{(nu-2)/nu}.
            ## - With Z_t ~ t_nu(0, (nu-2)/nu), X_t ~ t_nu(0, \sigma_t^2 (nu-2)/nu)
            ##   and \tilde{\sigma}_t^2 = \sigma_t^2 (nu-2)/nu, a likelihood contribution
            ##   is f_{t_nu(0, \tilde{\sigma}_t^2)}(X_t) = f_{t_nu}(X_t/\tilde{\sigma}_t)/\tilde{\sigma}_t
            sig.t.tilde <- sig.t * sqrt((df-2)/df)
            sum(dt(x/sig.t.tilde, df = df, log = TRUE) - log(sig.t.tilde)) # X.t/sig.t.tilde = Z_t/\sqrt{(\nu-2)/\nu}
        } else {
            -Inf
        }
    },
    stop("Wrong 'distr'"))
}

##' @title Fitting GARCH(1,1) Processes with a Reparameterized Likelihood
##' @param x data (e.g., risk-factor changes)
##' @param init numeric(2) giving the initial values for z_{cor}, z_{ema}
##' @param sig2 numeric(1), > 0, third parameter providing the annualized variance
##' @param delta unit of time (defaults to 1 = daily data; for yearly data, use 250)
##' @param distr standardized innovation distribution (distribution of Z;
##'        E(Z) = 0, Var(Z) = 1); available are N(0,1) and standardized t (t_nu(0, (nu-2)/nu))
##' @param control control list passed to optim()
##' @param ... additional arguments passed to optim()
##' @return a list with components
##'         coef: the estimated alpha_0, alpha_1 and beta_1 and (for distr == "st")
##'               the estimated df
##'         logLik: the log-likelihood (in the two parameters) at the estimated parameters
##'                 (z_{cor}, z_{ema})
##'         counts, convergence, message: see ?optim()
##'         sig.t: the conditional volatility \sigma_t
##'         Z.t: the standardized residuals Z_t
##' @author Marius Hofert
##' @note 1) No 'estimate.cov' here as it wouldn't be meaningful anyways (only two
##'          of the three GARCH parameters estimated per likelihood anyways)
##'       2) Code with ideas from Marcel Braeutigam (brautigam@essec.edu)
fit_GARCH_11 <- function(x, init = NULL, # z_{cor}, z_{ema}, (d.o.f. nu)
                         sig2 = mean(x^2), # \sigma_{ann}
                         delta = 1, distr = c("norm", "st"), control = list(), ...)
{
    ## Checks
    distr <- match.arg(distr)
    stopifnot(is.numeric(x), is.null(init) || length(init) == if(distr == "norm") 2 else 3,
              sig2 > 0, delta > 0, is.list(control))

    ## Compute initial values (hardcoded, not ideal)
    if(is.null(init)) init <- if(distr == "norm") c(1, 1) else c(1, 1, 4) # 4 d.o.f.

    ## Fitting
    mx2 <- mean(x^2)
    x2.tm1 <- c(mx2, head(x, n = -1)^2) # X_0^2 (sample mean of squared data), ..., X_{t-1}^2
    control <- c(as.list(control), fnscale = -1) # maximization (overwrites possible additionally passed 'fnscale')
    fit <- optim(init, fn = logLik_GARCH_11,
                 sig2 = sig2, x = x, x2.tm1 = x2.tm1, delta = delta,
                 distr = distr, control = control, ...)

    ## Reparametrize into alpha_0, alpha_1 and beta_1; compare (3) with (2)
    mu <- exp(-delta / exp(fit$par[1:2])) # fit$par[1:2] = c(z_{cor}, z_{ema})
    alpha0 <- mx2 * (1 - mu[1]) # \alpha_0 = \sigma^2 * (1-\mu_{cor})
    alpha1 <- mu[1] * (1 - mu[2]) # \alpha_1 = \mu_{cor} * (1-\mu_{ema})
    beta1  <- mu[1] * mu[2] # \beta_1 = \mu_{cor} * \mu_{ema}

    ## Extract \sigma_t^2 and the standardized residuals Z; see logLik_GARCH_11()
    nonrec <- sig2 * (1-mu[1]) + mu[1] * (1-mu[2]) * x2.tm1 # nonrecursive part (see above)
    sig2.t <- as.numeric(filter(nonrec, filter = prod(mu), "recursive", init = sig2)) # filter out \sigma_1^2, ..., \sigma_t^2 and multiply with mu[1] * mu[2] = \mu_{cor} * \mu_{ema}
    Z.t <- x / sqrt(sig2.t)

    ## Return
    list(coef = c(alpha0 = alpha0, alpha1 = alpha1, beta1 = beta1,
                  df = if(distr == "st") fit$par[3]), # estimated degrees of freedom (only for standardized t innovations)
         logLik = fit$value, # maximized log-likelihood; see ?optim
         counts = fit$counts, # see ?optim; number of calls to fn and gr
         convergence = fit$convergence, # see ?optim; '0' indicates successful completion
         message = fit$message, # see ?optim
         sig.t = sqrt(sig2.t), # conditional volatility \sigma_t
         Z.t = Z.t) # standardized residuals Z_t
}

##' @title GARCH(1,1) Tail Index
##' @param innovations realizations of the innovations Z to estimate the mean
##'        (obtained, e.g., via rnorm(n) or rt(n, df = nu) * sqrt((nu-2)/nu))
##' @param alpha1 GARCH(1,1) coefficient, >= 0 with alpha1 + beta1 < 1
##' @param beta1 GARCH(1,1) coefficient, >= 0 with alpha1 + beta1 < 1
##' @param interval initial interval for root finding
##' @param ... additional arguments passed to the underlying uniroot()
##' @return approximate tail index alpha (or NA if not found)
##' @author Marius Hofert
##' @note - E((alpha_1 * Z^2 + beta_1)^(alpha/2)) = 1 according to
##'         McNeil et al. (2015, p. 576); see also p. 118, Definition 4.20,
##'         p. 119, Proposition 4.21
##'       - For checking:
##'         a <- seq(0, 10, length.out = 101)
##'         y <- sapply(a, function(x) f(x))
##'         plot(a, y, type = "l")
##'       - Careful: The tail index formula implemented assumes include.mean = FALSE
##'         (so the GARCH(1,1) has no additional constant mean).
tail_index_GARCH_11 <- function(innovations, alpha1, beta1,
                                interval = c(1e-6, 1e2), ...)
{
    if(!(alpha1 >= 0 && beta1 >= 0 && alpha1 + beta1 < 1)) NA
    f <- function(a) mean((alpha1 * innovations^2 + beta1)^(a/2)) - 1
    f. <- f(interval[1])
    f.. <- f(interval[2])
    if(f. * f.. >= 0) NA else uniroot(f, interval = interval, ...)$root
}
