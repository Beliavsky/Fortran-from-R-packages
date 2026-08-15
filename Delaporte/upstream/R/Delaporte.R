# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

ddelap <- function(x, alpha, beta, lambda, log = FALSE) {
  x <- as.double(x)
  xvalid <- x[!(is.nan(x) | is.na(x))]
  if (any(xvalid > floor(xvalid))) {
    warning("Non-integers passed to ddelap. These will have 0 probability.")
  }
  if (log) log_f <- 1L else log_f <- 0L
  .Call(ddelap_C, x, as.double(alpha), as.double(beta), as.double(lambda),
        log_f, getDelapThreads())
}

pdelap <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  # These interrupts throw errors even using expect_error. Excluding for now
  # nocov start
  if (any(q[is.finite(q)] >= 2 ^ 63)) {
    stop("Function cannot handle values >= 2^63.")
  }
  if (any(q[is.finite(q)] >= 2 ^ 15)) {
    cat("There are values >= 32768.",
        "This may take minutes if not hours to compute. Are you sure?\n")
    resp <- readline("Press 'y' to continue.\n")
    if (tolower(resp) != "y") {
      cat("Stopping\n")
      return(invisible(NULL))
    }
  }
  # nocov end
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  .Call(pdelap_C, as.double(q), as.double(alpha), as.double(beta),
        as.double(lambda), lt_f, lp_f, getDelapThreads())
}

qdelap <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE,
                   exact = TRUE) {
  p <- as.double(p)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  if (exact) {
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f,
                   getDelapThreads())
  } else if (length(alpha) > 1L || length(beta) > 1L || length(lambda) > 1L ||
          anyNA(p)) {
    warning("Quantile approximation relies on pooling and is not accurate when",
            " passed vector-valued parameters, NaNs, or NAs. Using exact",
            " version instead.")
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f,
                   getDelapThreads())
  } else if (alpha <= 0 || beta <= 0 || lambda <= 0) {
    QDLAP <- rep.int(NaN, length(p))
  } else {
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      pValid <- p[p > 0 & p < 1]
      pNeg <- p[p < 0]
      p0 <- p[p == 0]
      pInf <- p[p >= 1]
      n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 5), 1e7)
      shiftedGammas <- rgamma(n, shape = alpha, scale = beta)
      DP <- rpois(n, lambda = (shiftedGammas + lambda))
      qValid <- as.vector(quantile(DP, pValid, na.rm = TRUE, type = 8))
      qNeg <- rep.int(NaN, times = length(pNeg))
      q0 <- rep.int(0, times = length(p0))
      qInf <- rep.int(Inf, times = length(pInf))
      QDLAP <- as.vector(c(qNeg, q0, qValid, qInf), mode = "double")
  }
  if (any(is.nan(QDLAP))) warning("NaNs produced")
  QDLAP
}

rdelap <- function(n, alpha, beta, lambda, exact = TRUE) {
  if (n < 0) {
    stop("invalid arguments")
  }
  n <- as.integer(n)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  if (exact) {
    RDLAP <- .Call(rdelap_C, n, alpha, beta, lambda, getDelapThreads())
  } else if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0)) {
      RDLAP <- (rep.int(NaN, n))
  } else {
    shiftedGammas <- rgamma(n, shape = alpha, scale = beta)
    RDLAP <- rpois(n, lambda = (shiftedGammas + lambda))
  }
  if (any(is.nan(RDLAP))) warning("NaNs produced")
  RDLAP
}

MoMdelap <- function(x, type = 2L) { # nolint object_name_linter
  type <- as.integer(type)
  if (!(type %in% c(1L, 2L, 3L))) stop("Skew type must be one of 1, 2, or 3.")
  moMDLAP <- .Call(MoMdelap_C, as.double(x), type)
  if (any(moMDLAP <= 0)) {
    stop("Method of moments not appropriate for this data; results include ",
         "non-positive parameters.")
  }
  moMDLAP
}
