## Note that the implementation of the weighted ML estimation is
## influenced and inspired by the function fitdistr from package MASS
## and function fitdist from package fitdistrplus.

FLXMCdist1 <- function(formula=.~., dist, ...) {
    foo <- paste("FLXMC", dist, sep = "")
    if (!exists(foo))
        stop("This distribution has not been implemented yet.")
    get(foo)(formula, ...)
}

prepoc.y.pos.1 <- function(x) {
    if (ncol(x) > 1)
        stop("for the inverse gaussian family y must be univariate")
    if (any(x < 0))
        stop("values must be >= 0")
    x
}

FLXMClnorm <- function(formula=.~., ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "lnorm", name="model-based log-normal clustering")
    z@preproc.y <- prepoc.y.pos.1
    
    z@defineComponent <- function(para) {
      predict <-  function(x, ...)
          matrix(para$meanlog, nrow = nrow(x), ncol = 1, byrow = TRUE)

      logLik <- function(x, y)
          dlnorm(y, meanlog = predict(x, ...), sdlog = para$sdlog, log = TRUE)

      new("FLXcomponent", parameters = list(meanlog = para$meanlog, sdlog = para$sdlog),
          predict = predict, logLik = logLik, df = para$df)
    }
    
    z@fit <- function(x, y, w, ...) {
        logy <- log(y); meanw <- mean(w)
        meanlog <-  mean(w * logy) / meanw
        sdlog <- sqrt(mean(w * (logy - meanlog)^2) / meanw)
        z@defineComponent(list(meanlog = meanlog, sdlog = sdlog, df = 2))
    }
    z
}

FLXMCinvGauss <- function(formula=.~., ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based inverse Gaussian clustering",
             dist = "invGauss")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$nu, nrow = nrow(x), ncol = length(para$nu),
                   byrow = TRUE)

        logLik <- function(x, y, ...)
            SuppDists::dinvGauss(y, nu = predict(x, ...), lambda = para$lambda, log = TRUE)

        new("FLXcomponent", parameters = list(nu = para$nu, lambda = para$lambda),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, ...){
        nu <- mean(w * y) / mean(w)
        lambda <- mean(w) / mean(w * (1 / y - 1 / nu))
        z@defineComponent(list(nu = nu, lambda = lambda, df = 2))
    }
    z
}

FLXMCgamma <- function(formula=.~., method = "Nelder-Mead", warn = -1, ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based gamma clustering",
             dist = "gamma")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$shape, nrow = nrow(x), ncol = length(para$shape),
                   byrow = TRUE)

        logLik <- function(x, y, ...)
            dgamma(y, shape = predict(x, ...), rate = para$rate, log = TRUE)
        
        new("FLXcomponent", parameters = list(shape = para$shape, rate = para$rate),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, component){
        if (!length(component)) {
            sw <- sum(w)
            mean <- sum(y * w) / sw
            var <- (sum(y^2 * w) / sw - mean^2) * sw / (sw - 1)
            start <- c(mean^2/var, mean/var)
        } else start <- unname(unlist(component))
        control <- list(parscale = c(1, start[2]))
        f <- function(parms) -sum(dgamma(y, shape = parms[1], rate = parms[2], log = TRUE) * w)
        oop <- options(warn = warn)
        on.exit(oop)
        parms <- optim(start, f, method = method, control = control)$par
        z@defineComponent(list(shape = parms[1], rate = parms[2], df = 2)) 
      }
    z
}
    
FLXMCexp <- function(formula=.~., ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based exponential clustering",
             dist = "exp")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$rate, nrow = nrow(x), ncol = length(para$rate),
                   byrow = TRUE)

        logLik <- function(x, y, ...)
            dexp(y, rate = predict(x, ...), log = TRUE)
        
        new("FLXcomponent", parameters = list(rate = para$rate),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, component)
        z@defineComponent(list(rate = mean(w) / mean(w * y), df = 1))
    z
}
        
FLXMCweibull <- function(formula=.~., method = "Nelder-Mead", warn = -1, ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based Weibull clustering",
             dist = "weibull")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$shape, nrow = nrow(x), ncol = length(para$shape),
                   byrow = TRUE)

        logLik <- function(x, y, ...) 
            dweibull(y, shape = predict(x, ...), scale = para$scale, log = TRUE)
        
        new("FLXcomponent", parameters = list(shape = para$shape, scale = para$scale),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, component){
        if (!length(component)) {
            ly <- log(y)
            sw <- sum(w)
            mean <- sum(ly * w) / sw
            var <- (sum(ly^2 * w) / sw - mean^2) * sw / (sw - 1)
            shape <- 1.2/sqrt(var)
            scale <- exp(mean + 0.572/shape)
            start <- c(shape, scale)
        } else start <- unname(unlist(component))
        f <- function(parms) -sum(dweibull(y, shape = parms[1], scale  = parms[2], log = TRUE) * w)
        oop <- options(warn = warn)
        on.exit(oop)
        parms <- optim(start, f, method = method)$par
        z@defineComponent(list(shape = parms[1], scale = parms[2], df = 2))
      }
    z
}

FLXMCburr <- function(formula=.~., start = NULL, method = "Nelder-Mead", warn = -1, ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based Burr clustering",
             dist = "burr")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$shape1, nrow = nrow(x), ncol = length(para$shape1),
                   byrow = TRUE)

        logLik <- function(x, y, ...)
            actuar::dburr(y, shape1 = predict(x, ...), shape2 = para$shape2, scale = para$scale, log = TRUE)
        
        new("FLXcomponent", parameters = list(shape1 = para$shape1, shape2 = para$shape2, scale = para$scale),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, component){
        if (!length(component)) {
            if (is.null(start)) start <- c(1, 1)
        } else start <- unname(unlist(component[2:3]))
        f <- function(parms) {
            shape1 <- sum(w) / sum(w * log(1 + (y/parms[2])^parms[1]))
            -sum(actuar::dburr(y, shape1 = shape1, shape2 = parms[1], scale = parms[2], log = TRUE) * w)
        }
        oop <- options(warn = warn)
        on.exit(oop)
        parms <- optim(start, f, method = method)$par
        z@defineComponent(list(shape1 = sum(w) / sum(w * log(1 + (y/parms[2])^parms[1])), shape2 = parms[1], scale = parms[2], df = 3))
      }
    z
}
    
FLXMCinvburr <- function(formula=.~., start = NULL, warn = -1, ...)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name = "model-based Inverse Burr clustering",
             dist = "invburr")
    z@preproc.y <- prepoc.y.pos.1

    z@defineComponent <- function(para) {
        predict <- function(x, ...) 
            matrix(para$shape1, nrow = nrow(x), ncol = length(para$shape1),
                   byrow = TRUE)

        logLik <- function(x, y, ...)
            actuar::dinvburr(y, shape1 = predict(x, ...), shape2 = para$shape2, scale = para$scale, log = TRUE)
        
        new("FLXcomponent", parameters = list(shape1 = para$shape1, shape2 = para$shape2, scale = para$scale),
            predict = predict, logLik = logLik, df = para$df)
    }

    z@fit <- function(x, y, w, component){
        if (!length(component)) {
            if (is.null(start)) start <- c(1, 1)
        } else start <- unname(unlist(component[2:3]))
        f <- function(parms) {
            shape1 <- sum(w) / sum(w * log(1 + (parms[2]/y)^parms[1]))
            -sum(actuar::dinvburr(y, shape1 = shape1, shape2 = parms[1], scale  = parms[2], log = TRUE) * w)
        }
        oop <- options(warn = warn)
        on.exit(oop)
        parms <- optim(start, f, method = "Nelder-Mead")$par
        z@defineComponent(list(shape1 = sum(w) / sum(w * log(1 + (parms[2]/y)^parms[1])), shape2 = parms[1], scale = parms[2], df = 3))
      }
    z
}
