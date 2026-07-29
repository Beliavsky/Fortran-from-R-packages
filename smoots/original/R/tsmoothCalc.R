tsmoothCalc <- function(y, p = c(1, 3), mu = c(0, 1, 2, 3),
                    Mcf = c("NP", "ARMA", "AR", "MA"),
                    InfR = c("Opt", "Nai", "Var"), bStart = 0.15,
                    bvc = c("Y", "N"), bb = c(0, 1), cb = 0.05,
                    method = c("lpr", "kr")) {

  # using the knn idea, bb=1, or not
  # cb=x, 0, 0.025 or 0.05, x*n obs at each end not be used for choosing b
  # In this code ccf and bStart will be provided by the user
  # Input parameters
  n <- length(y)
  k <- p + 1
  pd <- p + 2
  runc <- 1
  n1 <- trunc(n * cb)

  # New method for the kernel constants with p = 1 or 3
  # Kernel constants
  m <- 1000000  # for the numerical integral
  u <- (-m:m) / (m + 0.5)
  # For p=1, any kernel in the given c-mu class is possible
  if (p == 1) {
    wkp <- (1 - u^2)^(mu)  # a standardization is not necessary
  }

  # For p=3, the four 4-th order kernels (Table 5.7, Mueller, 1988)
  # table saved internally
  if (p == 3) {  # the constant factor does not play any role
    wkp <- lookup$p3_lookup[[mu + 1]](u)
  }
  Rp <- sum(wkp^2) / m
  mukp <- sum((u^k) * wkp) / m

  # Two constant in the bandwidth
  c1 <- factorial(k)^2 / (2 * k)
  c2 <- (1 - 2 * cb) * Rp / (mukp)^2

  steps <- rep(NA, 40)
  p.char <- as.character(p)
  bd_func <- lookup$InfR_lookup[[p.char, InfR]]
  bv_func <- lookup$bvc_lookup[[p.char, mu + 1]]
  expo1 <- lookup$expon1[[p.char]]
  expo2 <- lookup$expon2[[p.char]]

  # The main iteration------------------------------------------------------

  noi <- 40  # maximal number of iterations

  for (i in 1:noi) {
    if (runc == 1) {
      if (i > 1) {bold1 <- bold}
      if (i == 1) {bold <- bStart} else {bold <- bopt}

      # Look up the EIM inflation rate in the internal list 'lookup'
      bd <- bd_func(bold)

      if (bd >= 0.49) {bd <- 0.49}
      yed <- c(gsmoothCalcCpp(y, k, pd, mu, bd, bb))
      I2 <- sum(yed[(n1 + 1):(n - n1)]^2) / (n - 2 * n1)

      # Use an enlarged bandwidth for estimating cf or not (Feng/Heiler, 2009)
      if (bvc == "Y") {

        # Look up the enlarged bandwidth
        bv <- bv_func(bold)

      } else {bv <- bold}

      if (bv >= 0.49) {bv <- 0.49}

      ye <- c(gsmoothCalcCpp(y, 0, p, mu, bv, bb))

      # The optimal bandwidth-----------------------------------------------

      # Estimating the variance factor
      yd <- y - ye

      # Lag-Window for estimating the variance
      if (Mcf == "NP") {
        cf0.LW.est <- cf0Cpp(yd)  #[(n1+1):(n-n1)])
        cf0.LW <- cf0.LW.est$cf0.LW
        cf0 <- cf0.LW
        cf0.AR <- NA
        cf0.MA <- NA
        cf0.ARMA <- NA
        L0.opt <- cf0.LW.est$L0.opt
        ARMA.BIC <- NA
        AR.BIC <- NA
        MA.BIC <- NA
        p.BIC <- NA
        q.BIC <- NA
      } else if (Mcf == "ARMA") {
        cf0.ARMA.est <- cf0.ARMA.est(yd)
        cf0.ARMA <- cf0.ARMA.est$cf0.ARMA
        cf0.AR <- NA
        cf0.MA <- NA
        cf0 <- cf0.ARMA
        AR.BIC <- NA
        MA.BIC <- NA
        ARMA.BIC <- cf0.ARMA.est$ARMA.BIC
        p.BIC <- cf0.ARMA.est$p.BIC
        q.BIC <- cf0.ARMA.est$q.BIC
        cf0.LW <- NA
        L0.opt <- NA
      } else if (Mcf == "AR") {
        cf0.AR.est <- cf0.AR.est(yd)  # [(n1+1):(n-n1)])
        cf0.AR <- cf0.AR.est$cf0.AR
        cf0 <- cf0.AR
        cf0.LW <- NA
        cf0.MA <- NA
        cf0.ARMA <- NA
        AR.BIC <- cf0.AR.est$AR.BIC
        p.BIC <- cf0.AR.est$p.BIC
        ARMA.BIC <- NA
        MA.BIC <- NA
        q.BIC <- NA
        L0.opt <- NA
      } else if (Mcf == "MA") {
        cf0.MA.est <- cf0.MA.est(yd)
        cf0.MA <- cf0.MA.est$cf0.MA
        cf0 <- cf0.MA
        cf0.LW <- NA
        cf0.AR <- NA
        cf0.ARMA <- NA
        MA.BIC <- cf0.MA.est$MA.BIC
        q.BIC <- cf0.MA.est$q.BIC
        ARMA.BIC <- NA
        AR.BIC <- NA
        p.BIC <- NA
        L0.opt <- NA
      }

      c3 <- cf0 / I2

      bopt <- (c1 * c2 * c3)^expo1 * n^(-expo1)
      if (bopt < n^expo2) bopt <- n^expo2
      if (bopt > 0.49) {bopt = 0.49}
      steps[i] <- bopt
      if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc <- 0}
      if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n) {
        bopt <- (bold + bopt) / 2
        runc <- 0
      }
    }
  }

  # Smooth with the selected bandwidth--------------------------------------

  if (bopt < n^expo2) bopt <- n^expo2
  if (bopt > 0.49) {bopt <- 0.49}

  results <- list(b0 = bopt, cf0 = cf0, I2 = I2,
                  cf0.AR = cf0.AR, cf0.MA = cf0.MA, cf0.ARMA = cf0.ARMA,
                  cf0.LW = cf0.LW, L0.opt = L0.opt,
                  AR.BIC = AR.BIC, MA.BIC = MA.BIC, ARMA.BIC = ARMA.BIC,
                  p.BIC = p.BIC, q.BIC = q.BIC, n = n,
                  niterations = length(steps[!is.na(steps)]),
                  orig = y, iterations = steps[!is.na(steps)],
                  p = p, mu = mu, Mcf = Mcf, InfR = InfR, bStart = bStart,
                  bvc = bvc, bb = bb, cb = cb, v = 0)
  if (method == "lpr") {
    est.opt <- gsmoothCalc2Cpp(y, 0, p, mu, bopt, bb)
    ye <- c(est.opt$ye)
    attributes(ye) <- attributes(y)
    results[["ws"]] <- est.opt$ws
    res <- y - ye
    attr(results, "method") = "lpr"
  } else if (method == "kr") {
    est.opt <- knsmooth(y, mu, bopt, bb)
    res <- est.opt$res
    ye <- est.opt$ye
    attr(results, "method") = "kr"
  }
  results[["ye"]] <- ye
  results[["res"]] <- res
  results
}
