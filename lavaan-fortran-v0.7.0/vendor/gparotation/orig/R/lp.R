# Iterative re-weighted least squares for Lp rotation.
# Functions provided by Xinyi Liu, with random start integration by Coen Bernaards (February 2025).
#
# Each rotation uses an inner loop (GPA) and an outer loop (IRWLS).
# Random starts add a third outer loop.
#
# Note: lpQ and lpT exist as separate random-start wrappers rather than using
# .GPA_RS_engine because the Lp rotation calls GPFoblq.lp and GPForth.lp
# rather than GPFoblq and GPForth. Integrating with .GPA_RS_engine would
# require changes to the core GP algorithms, which was intentionally avoided.
#
# 4 functions:
#   lpQ        - oblique Lp rotation wrapper (with optional random starts)
#   GPFoblq.lp - main oblique Lp GP function
#   lpT        - orthogonal Lp rotation wrapper (with optional random starts)
#   GPForth.lp - main orthogonal Lp GP function

# --- OBLIQUE Lp ROTATION ---

lpQ <- function(A, Tmat = diag(ncol(A)), p = 1, normalize = FALSE,
                eps = 1e-05, maxit = 2000, randomStarts = 0, gpaiter = 5) {
  if (normalize)
    message("Normalization is not recommended for Lp rotation. Use may have unexpected results.")

  if (randomStarts > 0)
    Tmat <- Random.Start(ncol(A))

  r <- GPFoblq.lp(A, Tmat = Tmat, p, normalize = normalize, eps = eps,
                  maxit = maxit, gpaiter = gpaiter)

  if (randomStarts > 1) {
    Qvalues   <- sum(abs(r$loadings)^p) / nrow(A)
    Qmin      <- Qvalues
    Qconverged <- r$convergence

    for (inum in 2:randomStarts) {
      gpout <- GPFoblq.lp(A, Tmat = Random.Start(ncol(A)), p,
                          normalize = normalize, eps = eps,
                          maxit = maxit, gpaiter = gpaiter)
      Qval       <- sum(abs(gpout$loadings)^p) / nrow(A)
      Qvalues    <- c(Qvalues, Qval)
      Qconverged <- c(Qconverged, gpout$convergence)
      if (Qval < Qmin) {
        r    <- gpout
        Qmin <- Qval
      }
    }

    Qmin      <- eps * round(Qmin    * 1 / eps, 0)
    Qvalues   <- eps * round(Qvalues * 1 / eps, 0)
    randStartChar <- c(randomStarts      = randomStarts,
                       Converged         = sum(Qconverged),
                       atMinimum         = sum(Qvalues == Qmin),
                       localMins         = length(unique(Qvalues)))

    r <- list(loadings    = r$loadings,
              Phi         = r$Phi,
              Th          = r$Th,
              Table       = r$Table,
              method      = r$method,
              orthogonal  = FALSE,
              convergence = r$convergence,
              Gq          = r$Gq,
              randStartChar = randStartChar,
              Qvalues     = Qvalues)
    class(r) <- "GPArotation"
  }

  dimnames(r$Phi) <- list(colnames(r$loadings), colnames(r$loadings))
  r
}

GPFoblq.lp <- function(A, Tmat = diag(ncol(A)), p = 1, normalize = FALSE,
                        eps = 1e-05, maxit = 10000, gpaiter = 5) {
  # Oblique Lp gradient projection via iterative re-weighted least squares.
  # Inner loop: GPA with weighted least squares criterion (vgQ.lp.wls)
  # Outer loop: update weights until loadings converge

  j_sum      <- 0
  start_time <- proc.time()
  convergence <- FALSE

  # Corrected initial weights based on current rotated loadings
  W <- ((A %*% t(solve(Tmat)))^2 + eps)^(p / 2 - 1)

  for (it in 1:maxit) {
    r     <- suppressWarnings(
               GPFoblq(A, Tmat = Tmat, normalize = normalize, eps = eps,
                       maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W)))
    T_new <- r$Th
    L     <- r$loadings
    k     <- nrow(r$Table)
    j     <- r$Table[k, 1]
    j_sum <- j_sum + j

    W <- (L^2 + eps)^(p / 2 - 1)

    if (max(abs(L - A %*% solve(t(Tmat)))) < eps) {
      Tmat        <- T_new
      convergence <- TRUE
      break
    }
    Tmat <- T_new
  }

  Phi <- t(Tmat) %*% Tmat
  dimnames(Phi) <- list(colnames(A), colnames(A))

  Table <- data.frame(
    iter = it,
    f    = sum(abs(L)^p) / nrow(L),
    time = proc.time()[3] - start_time[3]
  )

  r <- list(
    loadings    = L,
    Phi         = Phi,
    Th          = Tmat,
    Table       = Table,
    method      = paste0("Lp rotation, p=", p),
    orthogonal  = FALSE,
    convergence = convergence
  )
  class(r) <- "GPArotation"
  r
}

# --- ORTHOGONAL Lp ROTATION ---

lpT <- function(A, Tmat = diag(ncol(A)), p = 1, normalize = FALSE,
                eps = 1e-05, maxit = 2000, randomStarts = 0, gpaiter = 5) {
  if (normalize)
    message("Normalization is not recommended for Lp rotation. Use may have unexpected results.")

  if (randomStarts > 0)
    Tmat <- Random.Start(ncol(A))

  r <- GPForth.lp(A, Tmat = Tmat, p, normalize = normalize, eps = eps,
                  maxit = maxit, gpaiter = gpaiter)

  if (randomStarts > 1) {
    Qvalues    <- sum(abs(r$loadings)^p) / nrow(A)
    Qmin       <- Qvalues
    Qconverged <- r$convergence

    for (inum in 2:randomStarts) {
      gpout <- GPForth.lp(A, Tmat = Random.Start(ncol(A)), p,
                          normalize = normalize, eps = eps,
                          maxit = maxit, gpaiter = gpaiter)
      Qval       <- sum(abs(gpout$loadings)^p) / nrow(A)
      Qvalues    <- c(Qvalues, Qval)
      Qconverged <- c(Qconverged, gpout$convergence)
      if (Qval < Qmin) {
        r    <- gpout
        Qmin <- Qval
      }
    }

    Qmin    <- eps * round(Qmin    * 1 / eps, 0)
    Qvalues <- eps * round(Qvalues * 1 / eps, 0)
    randStartChar <- c(randomStarts      = randomStarts,
                       Converged         = sum(Qconverged),
                       atMinimum         = sum(Qvalues == Qmin),
                       localMins         = length(unique(Qvalues)))

    r <- list(loadings    = r$loadings,
              Phi         = r$Phi,
              Th          = r$Th,
              Table       = r$Table,
              method      = r$method,
              orthogonal  = TRUE,
              convergence = r$convergence,
              Gq          = r$Gq,
              randStartChar = randStartChar,
              Qvalues     = Qvalues)
    class(r) <- "GPArotation"
  }

  dimnames(r$Phi) <- list(colnames(r$loadings), colnames(r$loadings))
  r
}

GPForth.lp <- function(A, Tmat = diag(ncol(A)), p = 1, normalize = FALSE,
                        eps = 1e-05, maxit = 10000, gpaiter = 5) {
  # Orthogonal Lp gradient projection via iterative re-weighted least squares.
  # Inner loop: GPA with weighted least squares criterion (vgQ.lp.wls)
  # Outer loop: update weights until loadings converge

  j_sum      <- 0
  start_time <- proc.time()
  convergence <- FALSE

  # Corrected initial weights based on current rotated loadings
  W <- ((A %*% t(solve(Tmat)))^2 + eps)^(p / 2 - 1)

  for (it in 1:maxit) {
    r     <- suppressWarnings(
               GPForth(A, Tmat = Tmat, normalize = normalize, eps = eps,
                       maxit = gpaiter, method = "lp.wls", methodArgs = list(W = W)))
    T_new <- r$Th
    L     <- r$loadings
    k     <- nrow(r$Table)
    j     <- r$Table[k, 1]
    j_sum <- j_sum + j

    W <- (L^2 + eps)^(p / 2 - 1)

    if (max(abs(L - A %*% solve(t(Tmat)))) < eps) {
      Tmat        <- T_new
      convergence <- TRUE
      break
    }
    Tmat <- T_new
  }

  Phi <- t(Tmat) %*% Tmat
  dimnames(Phi) <- list(colnames(A), colnames(A))

  Table <- data.frame(
    iter = it,
    f    = sum(abs(L)^p) / nrow(L),
    time = proc.time()[3] - start_time[3]
  )

  r <- list(
    loadings    = L,
    Phi         = Phi,
    Th          = Tmat,
    Table       = Table,
    method      = paste0("Lp rotation, p=", p),
    orthogonal  = TRUE,
    convergence = convergence
  )
  class(r) <- "GPArotation"
  r
}