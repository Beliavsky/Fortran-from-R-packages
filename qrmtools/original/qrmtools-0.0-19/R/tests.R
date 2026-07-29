### Statistical tests ##########################################################

### Tests of multivariate normality ############################################

##' @title Normality Test based on Distribution of Squared Mahalanobis Distances
##' @param x (n, d) data matrix
##' @param type the type of test to be used
##' @param dist distribution to check
##' @param ... additional arguments passed to the underlying type
##' @return htest object as returned by ad.test() or ks.test()
##' @author Marius Hofert
maha2_test <- function(x, type = c("ad.test", "ks.test"), dist = c("chi2", "beta"),
                       ...)
{
    if(!is.matrix(x))
        x <- as.matrix(x)
    type <- match.arg(type)
    dist <- match.arg(dist)
    ## Build squared Mahalanobis distances
    Xbar <- colMeans(x)
    S <- cov(x)
    D2 <- mahalanobis(x, center = Xbar, cov = S)
    ## Address different types
    d <- ncol(x)
    switch(type,
           "ad.test" = {
               switch(dist,
                      "chi2" = {
                          ad.test(D2, distr.fun = pchisq, df = d, ...)
                      },
                      "beta" = {
                          n <- nrow(x)
                          ad.test(D2 * n * (n-1)^2, distr.fun = pbeta,
                                  shape1 = d/2, shape2 = (n-d-1)/2, ...)
                      },
                      stop("Wrong 'dist'"))
           },
           "ks.test" = {
               switch(dist,
                      "chi2" = {
                          ks.test(D2, y = "pchisq", df = d, ...)
                      },
                      "beta" = {
                          n <- nrow(x)
                          ks.test(D2 * n * (n-1)^2, y = "pbeta",
                                  shape1 = d/2, shape2 = (n-d-1)/2, ...)
                      },
                      stop("Wrong 'dist'"))
           },
           stop("Wrong 'type'"))
}

##' @title Mardia's Test of Multivariate Normality
##' @param x (n, d) data matrix
##' @param type the type of test to be used ("skewness" is based on Mahalanobis angles,
##'        "kurtosis" on Mahalanobis distances).
##' @param method method for computing the Mahalanobis angles or distances
##' @return htest object
##' @author Marius Hofert
mardia_test <- function(x, type = c("kurtosis", "skewness"), method = c("direct", "chol"))
{
    ## Basics
    if(!is.matrix(x))
        x <- as.matrix(x)
    n <- nrow(x)
    d <- ncol(x)
    type <- match.arg(type)
    method <- match.arg(method)

    ## Main
    switch(type,
           "kurtosis" = {
               ## Build squared Mahalanobis distances
               D2 <- mahalanobis(x, center = colMeans(x), cov = cov(x))

               ## Compute the test results
               k <- mean(D2^2)
               T <- (k - d * (d + 2)) / sqrt(8 * d * (d + 2) / n)
               p.val <- pnorm(abs(T), lower.tail = FALSE)
               alternative <- "two-sided"
               strng <- paste0("Mardia's kurtosis test (computed with method = '",method,"')")
           },
           "skewness" = {
               ## Build Mahalanobis angles
               X.centered <- sweep(x, MARGIN = 2, STATS = colMeans(x)) # X_i - bar(X)
               S <- cov(x) # sample covariance matrix S
               D <- switch(method,
                           "direct" = {
                               S.inv <- solve(S) # S^{-1}; as in mahalanobis()
                               X.centered %*% S.inv %*% t(X.centered) # (n, d) %*% (d, d) %*% (d, n)
                           },
                           "chol" = { # similar to QRM::MardiaTest but more efficient
                               R <- chol(S)
                               R.inv <- solve(R)
                               Z <- X.centered %*% R.inv # (n, d)-matrix
                               Z %*% t(Z) # (n, d) %*% (d, n)
                           },
                           stop("Wrong 'method'"))

               ## Compute the test statistic results
               b <- mean(D^3)
               T <- (n/6) * b
               p.val <- pchisq(T, d * (d + 1) * (d + 2) / 6, lower.tail = FALSE)
               alternative <- "one-sided"
               strng <- paste0("Mardia's skewness test (computed with method = '",method,"')")
           },
           stop("Wrong 'type'"))

    ## Build and return htest object
    structure(class = "htest",
              list(statistic = c(statistic = T), # have to give it a name so that it shows up in output
                   p.value = p.val, # shows up with 'p-value'
                   alternative = alternative,
                   method = strng,
                   data.name = deparse(substitute(x)))) # name of the provided data object
}
