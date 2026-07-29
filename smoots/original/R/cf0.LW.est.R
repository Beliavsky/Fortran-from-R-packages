#'
#'@importFrom stats arima
#'

# Function for selecting M and estimating cf nonparametrically-----------------

# removed and now written in C++

#------------------------------------------------------------------------------


# Function for estimating cf0 from AR(0) to AR(p.max), e.g. p.max = 5----------

cf0.AR.est <- function(Xt) {
  p.max <- 5
  BIC <- (0:p.max) * 0
  n = length(Xt)
  for (l in 0:p.max) {
    AR.l <- arima(Xt, order = c(l, 0, 0))
    BIC[l + 1] <- -2 * AR.l$loglik + log(n) * (l + 2) - log(n) * 2
  }
  p.BIC <- (0:p.max)[BIC == min(BIC)]
  AR.BIC <- arima(Xt, order = c(p.BIC, 0, 0))
  sc.AR <- 0
  if (p.BIC > 0) {sc.AR <- sum(AR.BIC$coef[1:p.BIC])}
  cf0.AR <- (1 / (1 - sc.AR)) ** 2 * AR.BIC$sigma2

  results <- list(cf0.AR = cf0.AR, AR.BIC = AR.BIC, p.BIC = p.BIC)

  results
}
#### end of the cf0.AR.est function

### Estimation of optimal ARMA model and the resulting variance factor

cf0.ARMA.est <- function(Xt, p.max = 5, q.max = 5) {
  p.arma <- rep(0:p.max, times = q.max + 1)
  q.arma <- rep(0:q.max, each = p.max + 1)
  bic <- mapply(function(x, p0, q0) {
    arma.nobs <- length(x)
    arma.est <- arima(x, order = c(p0, 0, q0))
    -2 * arma.est$loglik + log(arma.nobs) * (p0 + q0)
  }, p0 = p.arma, q0 = q.arma,
  MoreArgs = list(x = Xt))
  opt <- which(bic == min(bic))
  p.BIC <- p.arma[opt]
  q.BIC <- q.arma[opt]
  ARMA.BIC <- arima(Xt, order = c(p.BIC, 0, q.BIC))
  if (p.BIC != 0) {
    sc.AR <- sum(ARMA.BIC$coef[1:p.BIC])
  } else {
    sc.AR <- 0
  }
  if (q.BIC != 0) {
    sc.MA <- sum(ARMA.BIC$coef[(p.BIC + 1):(p.BIC + q.BIC)])
  } else {
    sc.MA <- 0
  }
  cf0.ARMA <- ((1 + sc.MA) / (1 - sc.AR)) ^ 2 * ARMA.BIC$sigma2
  results <- list(cf0.ARMA = cf0.ARMA, ARMA.BIC <- ARMA.BIC,
                  p.BIC = p.BIC, q.BIC = q.BIC)

  results
}

### Estimation of the optimal MA model and the resulting variance factor

cf0.MA.est <- function(Xt, q.max = 5) {
  q.arma <- 0:q.max
  bic <- sapply(q.arma, FUN = function(x, q0) {
    ma.nobs <- length(x)
    ma.est <- arima(x, order = c(0, 0, q0))
    -2 * ma.est$loglik + log(ma.nobs) * q0
  }, x = Xt)
  opt <- which(bic == min(bic))
  q.BIC <- q.arma[opt]
  MA.BIC <- arima(Xt, order = c(0, 0, q.BIC))
  if (q.BIC != 0) {
    sc.MA <- sum(MA.BIC$coef[1:q.BIC])
  } else {
    sc.MA <- 0
  }
  cf0.MA <- (1 + sc.MA) ^ 2 * MA.BIC$sigma2
  results <- list(cf0.MA = cf0.MA, MA.BIC <- MA.BIC,
                  q.BIC = q.BIC)

  results
}
