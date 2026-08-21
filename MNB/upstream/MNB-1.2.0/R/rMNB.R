#' Generating Multivariate Negative Binomial Data
#'
#' @description It simulates a multivariate response variable, \eqn{Y_{ij}}, that is jth measurement
#' taken on the ith subject or cluster, i = 1,...,n and j= 1,...,mi.
#' @param n Length of the sample.
#' @param mi replicates on the ith subject or cluster.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones)
#' @param p.fix Vector of theoretical regression parameters of length p.
#' @return Generated response (\eqn{Y_{ij}})
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @examples
#'
#' n <- 100
#' mi <- 3
#' x1 <- rep(rnorm(n,0,1),each=mi)
#' x2 <- rep(c(0,1),each=150)
#' p.fix <- c(10,2.0,0.5,1)
#'
#' #generating a sample
#' sample.ex <- rMNB(n=n,mi=mi,formula=~x1+x2, p.fix=p.fix)
#' head(sample.ex)
#'
#' @export
#' @import flexsurv
#'
rMNB <- function(n, mi, formula, p.fix) {
    N = n * mi
    Y <- numeric(N)
    X <- stats::model.matrix(formula)
    p <- dim(X)[2]
    sigma <- (p.fix[1])^(-0.5)

    temp <- flexsurv::rgengamma(n = n, mu = 0, sigma = sigma,
        Q = sigma)
    ui <- log(temp)

    uij <- rep(ui, each = mi)
    eta <- X %*% (p.fix[2:(p + 1)])
    zij <- exp(eta + uij)
    Y <- stats::rpois(N, zij)
    ind = rep(1:n, each = mi)
    dat <- data.frame(Y = Y, X = X, ind = ind)
    return(dat)
}


