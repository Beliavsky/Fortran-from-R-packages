#' Residual analysis
#'
#' @description Weighted, standardized weighted, Pearson, standardized Pearson and
#' standardized deviance component residuals are available to assess possible departures
#' from the multivariate negative binomial model for fitting correlated data with overdispersion.
#' @param star Initial values for the parameters to be optimized over.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @details Similarly to GLMs theory (Agresti, 2015; Faraway, 2016), weighted and
#' the standardized weighted residuals are deduced trough Fisher scoring iterative process.
#' Based in the Pearson residual, Fabio (2017) suggest the standardized Pearson residuals
#' for the multivariate model from the random intercept Poisson-GLG model. In addition,
#' it is available the standardized deviance component residual for the ith subject
#' (Fabio et al., 2012).
#' @return Residuals
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Agresti, A. (2015). Foundations of Linear and Generalized Linear Models. Wiley.
#' \item Faraway, F. (2016). Extending the Linear Model with R: Generalized Linear,
#' Mixed Effects and nonparametric regression models. Taylor & Francis, New York.
#' \item Fabio, L., Paula, G. A., and de Castro, M. (2012). A Poisson mixed model
#' with nonormal random effect distribution. Computational Statistics and
#' Data Analysis, 56, 1499-1510.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#' }
#' @examples
#'
#' data(seizures)
#' head(seizures)
#'
#' star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)
#'
#' r <- re.MNB(formula=Y ~ trt + period + trt:period +
#' offset(weeks),star=star,dataSet=seizures)
#'
#' plot(r$ij.Sweighted.residual,cex.axis = 1.2, cex.lab = 1.2,
#' pch = 15,cex = 0.6, bg = 5,ylab="weighted.residual")
#'
#' abline(h=c(-3,0,3),lwd = 2, lty = 2)
#'
#' data(alzheimer)
#' head(alzheimer)
#'
#' star <- list(phi=10,beta1=2, beta2=0.2)
#' r <- re.MNB(formula = Y ~ trat,star=star,dataSet=alzheimer)
#' names(r)
#'
#' @export


re.MNB <- function(star, formula, dataSet) {

    op <- fit.MNB(star = star, formula = formula,
                  dataSet = dataSet,tab=FALSE)

    Y <- stats::model.response(data = stats::model.frame(formula,
                                                         dataSet))
    X <- stats::model.matrix(formula, dataSet)
    off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

    dataSet.ind <- split(dataSet, f = dataSet$ind)
    n <- length(dataSet.ind)
    p <- dim(X)[2]
    mi <- dim(dataSet.ind[[1]])[1]
    N <- n * mi

    ind <- dataSet$ind
    phi_hat <- op$par[1]
    beta_hat <- op$par[2:(p + 1)]

    eta_hat <- X %*% beta_hat
    if(methods::is(off,"NULL")){
      mu_hat <- exp(X %*% beta_hat)}else{mu_hat <- exp(X %*% beta_hat+off)}

    YI <- numeric(n)
    for (i in 1:n) {
      YI[i] <- sum(Y[ind == i])
    }

    mu_hati <- numeric(n)
    for (i in 1:n) {
      mu_hati[i] <- sum(mu_hat[ind == i])
    }

    W  <- diag(as.numeric(mu_hat)) + (rep((phi_hat+mu_hati)^(-1),each=mi)*mu_hat) %*% t(mu_hat)

    wi <- numeric(n)
    for (i in 1:n) {
      wi[i] <- sum(diag(W)[ind == i])
    }

    D <- matrix(NA,n,mi)
    i <- 1
    cont <- 1
    limit <- N-mi+2
    while(i<limit){
      l <- i
      c <- i+mi-1
      A <- sqrt(W[l:c,l:c])%*%X[l:c,1:ncol(X)]
      inv=solve(t(X)%*%W%*%X)
      D[cont,] <- diag(A%*%inv%*%t(A))
      i <- i+mi
      cont <- cont+1
    }

    h.ij <- c(t(D))


    # Residuals

    #-------------------------------------
    # ij- subject
    # ....................................
    # weighted residual - ij
    # ....................................

    a.i <- rep((phi_hat+YI)/(phi_hat+mu_hati),each=mi)
    wr.ij <- (Y-a.i*mu_hat)/sqrt(diag(W))

    # ....................................
    # Standardized weighted residual -ij
    # ....................................

    swr.ij <- wr.ij/sqrt(1-h.ij)

    # ......................................
    # Pearson.residual -ij
    # ......................................

    nump = phi_hat^(1/2)*(Y - mu_hat)
    denp = sqrt(mu_hat*(phi_hat+mu_hat))
    p.ij = nump/denp


    # ....................................
    # Standardized.Pearson.residual - ij
    # ...................................

    sp.ij = p.ij/sqrt(1 - h.ij)



    # ...........................
    # deviance residuals - i
    # ...........................

    aux.Y <- matrix(Y,ncol=mi,nrow=n,byrow = T)
    aux.mu_hat <- matrix(mu_hat,ncol=mi,nrow=n,byrow = T)
    aux1 <- numeric()
    for (i in 1:n) {
      temp1 <- numeric()
      for(j in 1:mi){
        if(aux.Y[i,j]==0){temp1[j] <- 0}else{
          auxN <- aux.Y[i,j] * (phi_hat + mu_hati[i])
          auxD <- aux.mu_hat[i,j] * (phi_hat + YI[i])
          auxLog <- log(auxN/auxD)
          temp1[j] <- aux.Y[i,j]*auxLog
        }
      }
      aux1[i] <- sum(temp1)
    }
    aux2Log <- log((phi_hat + mu_hati)/(phi_hat + YI))
    aux2 <- phi_hat*aux2Log
    d2 = 2*(aux2 + aux1)

    hii <- numeric(n)
    for (i in 1:n) {
      hii[i] <- mean(h.ij[ind == i])
    }

    numd = sign(YI - mu_hati)*sqrt(d2)
    dend = sqrt(1 - hii)

    d.i = numd/dend

    # ..................................................
    # Exit
    # .................................................

    result <- list()
    result$ij.weighted.residual <- as.numeric(wr.ij)
    result$ij.Sweighted.residual <- as.numeric(swr.ij)
    result$ij.pearson.residual <- as.numeric(p.ij)
    result$ij.Spearson.residual <- as.numeric(sp.ij)
    result$Deviance.residual <- as.numeric(d.i)


    return(result)
  }
