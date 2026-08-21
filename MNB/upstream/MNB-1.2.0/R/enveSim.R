#' Simulation envelope
#'
#' @description Simulated envelopes in normal probability plots
#' @param star Initial values for the parameters to be optimized over.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @param n.r Indicator which residual type graphics. 1 - weighted, 2 - Standardized weighted, 3 - Pearson, 4 - Standardized Pearson, 5 - standardized deviance component residuals and 6 - randomized quantile residuals.
#' @param nsim Number of Monte Carlo replicates.
#' @param plot TRUE or FALSE. Indicates if a graph should be plotted.
#' @details Atkinson (1985), suggests the use of simulated envelopes in normal probability
#' plots to facilitate the goodness of fit.
#' @return L, residuals and simulation envelopes in normal probability plots
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Atkinson A.C. (1985). Plots, Transformations and Regression: An Introduction to Graphical Methods of Diagnostic
#' Regression Analysis. Oxford University Press, New York.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#'
#' }
#' @examples
#'
#' \donttest{
#'
#' data(seizures)
#' head(seizures)
#'
#' star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)
#'
#' envelope.MNB(formula=Y ~ trt + period + trt:period +
#' offset(weeks),star=star,nsim=21,n.r=6,
#' dataSet=seizures,plot=FALSE)
#'
#' data(alzheimer)
#' head(alzheimer)
#'
#' star <- list(phi=10,beta1=2, beta2=0.2)
#' envelope.MNB(formula=Y ~ trat, star=star, nsim=21, n.r=6,
#' dataSet = alzheimer,plot=FALSE)
#'
#' }
#' @export

envelope.MNB <- function(star, formula, dataSet, n.r, nsim, plot=TRUE) {
    Y <- stats::model.response(data = stats::model.frame(formula,
        dataSet))
    X <- stats::model.matrix(formula, dataSet)
    off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

    dataSet.ind <- split(dataSet, f = dataSet$ind)

    n <- length(dataSet.ind)
    p <- dim(X)[2]
    mi <- dim(dataSet.ind[[1]])[1]
    N <- n * mi

    op <- fit.MNB(star = star, formula = formula,
        dataSet = dataSet,tab=FALSE)

    r <- re.MNB(star = star, formula = formula,
        dataSet = dataSet)

    r.quantil <- qMNB(par=op$par,formula=formula,dataSet=dataSet)

    if(n.r==6) {tp <- r.quantil}else{tp <- r[[n.r]]}

    e <- matrix(NA, length(tp), nsim)

    Y <- numeric(N)

    sigma <- (op$par[1])^(-0.5)

    for (k in 1:nsim) {
        ui <- log(flexsurv::rgengamma(n = n, mu = 0, sigma = sigma,
            Q = sigma))
        uij <- rep(ui, each = mi)
        eta <- X %*% (op$par[2:(p + 1)])
        if(methods::is(off,"NULL")){
        zij <- exp(eta + uij)}else{zij <- exp(eta + uij+off)}

        Y <- stats::rpois(N, zij)

        newDataSet <- data.frame(Y = Y, dataSet[,
            2:ncol(dataSet)])

        opEnv <- fit.MNB(star = star, formula = formula,
            dataSet = newDataSet,tab=FALSE)
        r.boot <- re.MNB(star = star, formula = formula,
            dataSet = newDataSet)

        r.qBoot <- qMNB(par=opEnv$par,formula=formula,dataSet=dataSet)

        if(n.r==6){ tp.boot <- r.qBoot}else{tp.boot <- r.boot[[n.r]]}
        e[, k] <- sort(tp.boot)
    }

    e1 <- apply(e, 1, min)
    e2 <- apply(e, 1, max)
    med <- apply(e, 1, mean)
    faixa <- range(tp, e1, e2,na.rm=TRUE,finite=TRUE)

    result <- list()
    result$mE <- cbind(e1, med, e2)
    result$residual <- tp

    if (plot==TRUE) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")

    # Plot - envelope
        stats::qqnorm(sort(tp), xlab = "Normal quantiles", ylab = "residual",
        pch = 15, ylim = faixa, main = "", cex.axis = 1.2,
        cex.lab = 1.2, cex = 0.6, bg = 5)
        graphics::par(new = T)
        stats::qqnorm(e1, axes = F, xlab = "", ylab = "", type = "l",
        ylim = faixa, lty = 1, main = "")
        graphics::par(new = T)
        stats::qqnorm(e2, axes = F, xlab = "", ylab = "", type = "l",
        ylim = faixa, lty = 1, main = "")
        graphics::par(new = T)
        stats::qqnorm(med, axes = F, xlab = "", ylab = "", type = "l",
        ylim = faixa, lty = 2, main = "")

    # Plot - residual
        graphics::plot(tp, ylab = "residual", xlab = "Index", ylim = faixa, cex.axis = 1.2, cex.lab = 1.2, pch = 15,
        cex = 0.6, bg = 5)
        graphics::abline(h = c(-3, 0, 3), lwd = 2, lty = 2)}

    return(result)
}
# ................................................................................................
