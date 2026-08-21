#' @importFrom grDevices gray
#' @importFrom graphics abline axis box image lines par plot points text title
#' @importFrom stats arima ARMAtoMA anova ar ar.ols as.ts ccf coef dchisq deltat df.kernel dnorm
#' @importFrom stats end filter frequency is.ts lm lm.fit ls.print lsfit makeARIMA na.fail na.omit
#' @importFrom stats na.pass optim pchisq pf pnorm predict qchisq qnorm quantile residuals rnorm
#' @importFrom stats rstandard spec.ar spec.pgram start time ts tsdiag tsp "tsp<-" var window
#' @importFrom mgcv gam
#' @importFrom leaps regsubsets
#' @importFrom locfit locfit
#' @export 
acf <-
function (x, lag.max = NULL, type = c("correlation", "covariance", 
    "partial")[1], plot = TRUE, na.action = na.fail, demean = TRUE, 
    drop.lag.0 = TRUE,  ...) 
{
    acf.out <- stats::acf(x=x, lag.max=lag.max, type=type, plot=F, na.action=na.action,
        demean=demean,...)
    acf.out$series <- deparse(substitute(x))
    if (drop.lag.0) {
        if (type == "correlation") {
            acf.out$acf = acf.out$acf[-1, , , drop = FALSE]
            acf.out$lag = acf.out$lag[-1, , , drop = FALSE]
        }
    }
        if (plot) {
        plot1.acf(acf.out, ...)
        return(invisible(acf.out))
    }
    else return(acf.out)
}
