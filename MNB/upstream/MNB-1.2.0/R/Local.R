
#------------------------------------------------
# Case weights perturbation - likelihood
#------------------------------------------------

cases <- function(star, formula, dataSet) {
    Y <- stats::model.response(data = stats::model.frame(formula,dataSet))
    X <- stats::model.matrix(formula, dataSet)
    off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

    dataSet.ind <- split(dataSet, f = dataSet$ind)
    n <- length(dataSet.ind)
    p <- dim(X)[2]
    mi <- dim(dataSet.ind[[1]])[1]
    N <- n * mi

    phi <- star[1]
    beta <- star[2:(p + 1)]
    w <- star[(p + 2):(n + p + 1)]

    if(methods::is(off,"NULL")){
    eta <- X %*% beta} else{eta <- X %*% beta + off}
    YI <- apply(matrix(Y, mi, n), 2, sum)
    mu.i <- apply(matrix(exp(eta), mi, n), 2, sum)

    logver <- (sum(unlist(w) * lgamma(phi + YI)) -
        sum(unlist(w) * lgamma(phi)) - sum(unlist(w) *
        lfactorial(Y)) + sum(unlist(w) * phi * (log(phi))) +
        t(unlist(w) * Y) %*% eta - sum(unlist(w) *
        (phi + YI) * log(phi + mu.i)))
    logver
}

#------------------------------------------------
# Explanatory variable perturbation
#------------------------------------------------

cases.obs <- function(star, formula, dataSet) {
    Y <- stats::model.response(data = stats::model.frame(formula,
        dataSet))
    X <- stats::model.matrix(formula, dataSet)
    dataSet.ind <- split(dataSet, f = dataSet$ind)
    off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

    n <- length(dataSet.ind)
    p <- dim(X)[2]
    mi <- dim(dataSet.ind[[1]])[1]
    N <- n * mi

    phi <- star[1]
    beta <- star[2:(p + 1)]
    w <- star[(p + 2):(N + p + 1)]

    if(methods::is(off,"NULL")){
      eta <- X %*% beta} else{eta <- X %*% beta + off}

    YI <- rep(apply(matrix(Y, mi, n), 2, sum), each = mi)
    mu.i <- rep(apply(matrix(exp(eta), mi, n), 2,
        sum), each = mi)

    logver <- numeric()

    for (i in 1:N) {
        func1 <- (1/mi) * lgamma(phi + YI[i]) - (1/mi) *
            lgamma(phi) - lfactorial(Y[i]) + (1/mi) *
            phi * log(phi) - (1/mi) * phi * log(phi +
            mu.i[i]) + Y[i] * eta[i] - Y[i] * log(phi +
            mu.i[i])
        logver[i] <- w[i] * func1
    }

    sum(logver)
}

cova.pertu<- function (star,formula,dataSet,cova){
  Y <-  stats::model.response(data=stats::model.frame(formula, dataSet))
  X <-  stats::model.matrix(formula, dataSet)
  off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

  dataSet.ind <- split(dataSet, f= dataSet$ind)
  n <-  length(dataSet.ind)
  p <-  dim(X)[2]
  mi <- dim(dataSet.ind[[1]])[1]
  N <-  n*mi

  phi <- star[1]
  beta <- star[2:(p+1)]
  w <- star[(p+2):(n+p+1)]

  X.star<-X[,cova]+w*stats::sd(unique(X[,cova]))

  if(p==cova){X.new <- (cbind(X[,1:(cova-1)],X.star))}
  else{X.new <- (cbind(X[,1:(cova-1)],X.star,X[,(cova+1):p]))}

  if(methods::is(off,"NULL")){
    eta <- X %*% beta} else{eta <- X %*% beta + off}

  YI <- apply(matrix(Y,mi,n),2,sum)
  mu.i <- apply(matrix(exp(eta),mi,n), 2,sum)

  logver <- (
    sum (lgamma(phi+YI))-sum(lgamma(phi))-
      sum(lfactorial(Y))+sum(phi*(log(phi)))+
      t(Y)%*%eta- sum((phi+YI)*log(phi+mu.i))
  )
  logver
}

#--------------------------------------------------
# Dispersion parameters perturbation
#--------------------------------------------------

dispersion <- function(star, formula, dataSet) {
    Y <- stats::model.response(data = stats::model.frame(formula,
        dataSet))
    X <- stats::model.matrix(formula, dataSet)
    off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

    dataSet.ind <- split(dataSet, f = dataSet$ind)
    n <- length(dataSet.ind)
    p <- dim(X)[2]
    mi <- dim(dataSet.ind[[1]])[1]
    N <- n * mi

    phi <- star[1]
    beta <- star[2:(p + 1)]
    w <- star[(p + 2):(n + p + 1)]

    if(methods::is(off,"NULL")){
      eta <- X %*% beta} else{eta <- X %*% beta + off}

    YI <- apply(matrix(Y, mi, n), 2, sum)
    mu.i <- apply(matrix(exp(eta), mi, n), 2, sum)

    logver <- (sum(lgamma(unlist(w) * phi + YI)) -
        sum(lgamma(unlist(w) * phi)) - sum(lfactorial(Y)) +
        sum(unlist(w) * phi * (log(unlist(w) * phi))) +
        t(Y) %*% eta - sum((unlist(w) * phi + YI) *
        log(unlist(w) * phi + mu.i)))
    logver
}


#' Local influence
#'
#' @description It performes influence analysis by a local influence approach by Cook (1986). It is considering
#' three perturbation schemes: Case weights, explanatory variable and dispersion parameter perturbation. Another
#' procedure which considering is the total local curvature corresponding to the ith element approach by
#' Lesaffre and Verbeke (1998).
#' @param star Initial values for the parameters to be optimized over.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @param schemes Perturbation scheme. Possible values: "cases" for Case weights perturbation on ith subject
#' or cluster, "cases.obs" for Case weights perturbation on jth measurement
#' taken on the ith subject or cluster, "cova.pertu" for explanatory variable perturbation, "dispersion" for
#' dispersion parameter perturbation
#' @param cova Indicator which column from dataset (continuous covariate) must be perturbation.
#' @param plot TRUE or FALSE. Indicates if a graph should be plotted.
#' @details The function returns a list (L) with the eigenvector associated with the maximum curvature, the total
#' local influence and the index plot.
#' @return L and graphics
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Cook, R. D. (1986). Assessment of local influence (with discussion). Journal of
#' the Royal Statistical Society B, 48, 133-169.
#' \item Lesaffre E. and Verbeke G. (1998). Local influence in linear mixed models. Biometrics, 54, 570-582.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
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
#' local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
#' schemes="weight",plot=FALSE)
#'
#' local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
#' schemes="weight.obs",plot=FALSE)
#'
#' local.MNB(formula=Y ~ trt + period + trt:period + offset(log(weeks)),star=star,dataSet=seizures,
#' schemes="dispersion",plot=FALSE)
#'
#' }
#' @export
#' @import numDeriv

local.MNB <- function(star, formula, dataSet, schemes,cova,plot=TRUE) {
    dataSet.ind <- split(dataSet, f = dataSet$ind)
    nID <- length(dataSet.ind)
    mi <- dim(dataSet.ind[[1]])[1]
    n <- length(dataSet.ind)
    N <- n * mi

    if (schemes == "weight") {
        es <- fit.MNB(star = star, formula = formula,
            dataSet = dataSet,tab=FALSE)
        parPert <- c(es$par, rep(1, nID))
        mhessW <- numDeriv::hessian(cases, x = parPert, formula = formula,
            dataSet = dataSet)
        delta <- mhessW[1:length(es$par), (length(es$par) +
            1):(nID + length(es$par))]
        FF <- t(delta) %*% solve(es$hessian) %*%
            delta
        Ci <- diag(FF)
        Cdmax <- eigen(FF)$values[nID]
        vect <- eigen(FF)$vector[, nID]

    }

    if (schemes == "weight.obs") {
        es <- fit.MNB(star = star, formula = formula,
            dataSet = dataSet,tab=FALSE)
        parPert <- c(es$par, rep(1, N))
        mhessW <- numDeriv::hessian(cases.obs, x = parPert,
            formula = formula, dataSet = dataSet)
        delta <- mhessW[1:length(es$par), (length(es$par) +
            1):(N + length(es$par))]
        FF <- t(delta) %*% solve(es$hessian) %*%
            delta
        Ci <- diag(FF)
        Cdmax <- eigen(FF)$values[N]
        vect <- eigen(FF)$vector[, N]

    }

    if(schemes=="covariate"){
      es <- fit.MNB(star=star,formula=formula,
                    dataSet=dataSet,tab=FALSE)
      parPert <- c(es$par,rep(0,nID))
      mhessW <- numDeriv::hessian(cova.pertu,x=parPert,formula=formula,
                        dataSet=dataSet,cova=cova)
      delta <- mhessW[1:length(es$par),(length(es$par)+1):(nID+length(es$par))]
      FF <- t(delta)%*%solve(es$hessian)%*%delta
      Ci <- diag(FF)
      Cdmax <- eigen(FF)$values[nID]
      vect <- eigen(FF)$vector[,nID]

    }

    if (schemes == "dispersion") {
        es <- fit.MNB(star = star, formula = formula,
            dataSet = dataSet,tab=FALSE)
        parPert <- c(es$par, rep(1, nID))
        mhessW <- numDeriv::hessian(dispersion, x = parPert,
            formula = formula, dataSet = dataSet)
        delta <- mhessW[1:length(es$par), (length(es$par) +
            1):(nID + length(es$par))]
        FF <- t(delta) %*% solve(es$hessian) %*%
            delta
        Ci <- diag(FF)
        Cdmax <- eigen(FF)$values[nID]
        vect <- eigen(FF)$vector[, nID]
    }

    result <- list()

    if (schemes == "weight") {
        result$local <- vect
        result$localTotal <- Ci
        if (plot==TRUE) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
          graphics::plot(abs(vect), main = "", ylab = "|dmax|",
            cex.axis = 1.2, cex.lab = 1.2, pch = 15,
            cex = 0.6, bg = 5)
          graphics::plot(abs(Ci), main = "", ylab = "|Ci|", cex.axis = 1.2,
            cex.lab = 1.2, pch = 15, cex = 0.6, bg = 5)}
    }

    if (schemes == "weight.obs") {
        result$local <- vect
        result$localTotal <- Ci
        if (plot==TRUE) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
          graphics::plot(abs(vect), main = "", ylab = "|dmax|",
            cex.axis = 1.2, cex.lab = 1.2, pch = 15,
            cex = 0.6, bg = 5)
          graphics::plot(abs(Ci), main = "", ylab = "|Ci|", cex.axis = 1.2,
            cex.lab = 1.2, pch = 15, cex = 0.6, bg = 5)}
    }

    if (schemes == "covariate") {
        result$local <- vect
        result$localTotal <- Ci
        if (plot==TRUE) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
          graphics::plot(abs(vect), main = "", ylab = "|dmax|",
            cex.axis = 1.2, cex.lab = 1.2, pch = 15,
            cex = 0.6, bg = 5)
          graphics::plot(abs(Ci), main = "", ylab = "|Ci|", cex.axis = 1.2,
            cex.lab = 1.2, pch = 15, cex = 0.6, bg = 5)}
    }

    if (schemes == "dispersion") {
        result$local <- vect
        result$localTotal <- Ci
        if (plot==TRUE) {
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))
          graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
          graphics::plot(abs(vect), main = "", ylab = "|dmax|",
            cex.axis = 1.2, cex.lab = 1.2, pch = 15,
            cex = 0.6, bg = 5)
          graphics::plot(abs(Ci), main = "", ylab = "|Ci|", cex.axis = 1.2,
            cex.lab = 1.2, pch = 15, cex = 0.6, bg = 5)}
    }

    return(result)

}
