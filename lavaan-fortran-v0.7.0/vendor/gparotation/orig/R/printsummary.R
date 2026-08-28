# S3 class functions for GPArotation objects:
#   .sortGPALoadings             --- internal helper, not exported
#   print.GPArotation
#   summary.GPArotation
#   print.summary.GPArotation

.sortGPALoadings <- function(x) {
  # Sort and sign-correct a GPArotation object.
  # Factors are reordered by descending variance explained and signs
  # are adjusted so that the sum of loadings per factor is positive.
  # Adapted from factanal function sortLoadings (R Core Team).
  # Called internally by print.GPArotation and summary.GPArotation.
  #
  # Args:
  #   x : a GPArotation object
  #
  # Returns:
  #   x with loadings, Th, and Phi updated consistently

  # Preserve original column names for restoration after reordering
  cln <- colnames(x$loadings)

  # Variance explained per factor:
  # orthogonal: sum of squared loadings per column
  # oblique: diagonal of Phi %*% t(L) %*% L (accounts for factor correlations)
  vx <- if (x$orthogonal)
    colSums(x$loadings^2)
  else
    diag(x$Phi %*% crossprod(x$loadings))

  # Save covariance attribute before subsetting strips it
  cov_attr <- attr(x$loadings, "covariance")

  # --- Legacy compatibility ---
  # Objects created before 2027 have no covariance attribute.
  # Reconstruct it from Phi (oblique) or identity (orthogonal).
  if (is.null(cov_attr)) {
    cov_attr <- if (x$orthogonal)
      diag(ncol(x$loadings))
    else if (!is.null(x$Phi))
      x$Phi
    else
      diag(ncol(x$loadings))
  }
  
  # Reorder factors from largest to smallest variance explained
  io.ssq     <- order(vx, decreasing = TRUE)
  x$loadings <- x$loadings[, io.ssq, drop = FALSE]

  # Sign correction: flip any factor column whose loadings sum to negative
  neg               <- colSums(x$loadings) < 0
  x$loadings[, neg] <- -x$loadings[, neg]

  # Restore column names after reordering
  colnames(x$loadings) <- cln

  # unit vector: 1 for columns kept as-is, -1 for sign-flipped columns
  unit <- c(1, -1)[neg + 1L]

  # Restore class and covariance attribute after subsetting
  class(x$loadings) <- "loadings"
  attr(x$loadings, "covariance") <- if (x$orthogonal)
    cov_attr                                       # identity - unaffected by reordering
  else
    diag(unit) %*% cov_attr[io.ssq, io.ssq] %*% diag(unit)  # reorder and sign-correct Phi

  # Apply reordering and sign correction to rotation matrix Th
  x$Th <- x$Th[, io.ssq] %*% diag(unit)
  
  # We must treat the starting matrix exactly like the final one  
  if (!is.null(x$Tinit)) {
    x$Tinit <- x$Tinit[, io.ssq, drop = FALSE] %*% diag(unit)
  }

  # Apply reordering and sign correction to Phi (oblique only)
  if (!x$orthogonal && "Phi" %in% names(x))
    x$Phi <- diag(unit) %*% x$Phi[io.ssq, io.ssq] %*% diag(unit)

  x
}
print.GPArotation <- function(x, digits = 3L, sortLoadings = TRUE, cutoff = .1,
                              rotateMat = FALSE, Table = FALSE, ...) {
  # --- Sort factors by descending variance explained ---
  if (sortLoadings) x <- .sortGPALoadings(x)

  # --- Convergence header ---
  cat(if (x$orthogonal) "Orthogonal" else "Oblique")
  cat(" rotation method", x$method)
  cat(if (!x$convergence) " NOT")
  cat(" converged")
  cat(if ("randStartChar" %in% names(x)) " at lowest minimum.\n" else ".\n")

  # Only show algorithm/fwindow if stored and non-default
  if (!is.null(x$algorithm) && !is.null(x$fwindow)) {
    if (!identical(x$algorithm, "bb") || x$fwindow != 10L) {
      cat("Algorithm:", x$algorithm, " |  fwindow:", x$fwindow, "\n")
    }
  }
    
  # --- Random start diagnostics ---
  if ("randStartChar" %in% names(x)) {
    rs <- x$randStartChar
    cat("Of ", rs[1], " random starts ",
        round(100 * rs[2] / rs[1]), "% converged, ",
        round(100 * rs[3] / rs[1]), "% at the same lowest minimum.\n", sep = "")
    if (rs[4] > 1)
      cat("Random starts converged to ", rs[4], " different local minima.\n", sep = "")
  }

  # --- AUC and FSI --- must be computed AFTER sorting, BEFORE printing ---
  auc <- calc_AUC(x$loadings)
  fsi <- calc_FSI(x$loadings)

  # --- Loadings + stats block ---
  L         <- x$loadings
  L_rounded <- round(L, digits)
  Lambda    <- unclass(L_rounded)
  p         <- nrow(Lambda)

  # Print loadings matrix
  cat(if("randStartChar" %in% names(x)) 
    "Loadings at lowest minimum:\n" 
    else "Loadings:\n")
  fx <- setNames(format(round(Lambda, digits)), NULL)
  nc <- nchar(fx[1L], type = "c")
  fx[abs(Lambda) < cutoff] <- strrep(" ", nc)
  print(fx, quote = FALSE)

  # Build varex block
  Phi <- attr(L, "covariance")
  vx  <- if (is.null(Phi) || identical(Phi, diag(ncol(Lambda))))
    colSums(Lambda^2)
  else
    diag(Phi %*% crossprod(Lambda))
  names(vx) <- colnames(Lambda)

  varex <- rbind(`SS loadings` = vx)

  if (x$orthogonal) {
    varex <- rbind(varex, `Proportion Var` = vx / p)
    if (ncol(Lambda) > 1)
      varex <- rbind(varex, `Cumulative Var` = cumsum(vx / p))
  }

  # Append AUC and FSI
  auc <- calc_AUC(x$loadings)
  fsi <- calc_FSI(x$loadings)
  varex <- rbind(varex, AUC = auc$AUC, FSI = fsi$FSI)

  cat("\n")
  print(round(varex, digits))
  
  # --- Phi (oblique only) ---
  if (!x$orthogonal) {
    dimnames(x$Phi) <- list(colnames(x$loadings), colnames(x$loadings))
    cat("\nPhi:\n")
    print(round(x$Phi, digits))
  }

  # --- Optional rotation matrix ---
  if (rotateMat) {
    cat("\nRotating matrix:\n")
    print(t(solve(x$Th)), digits = digits)
  }

  # --- Optional iteration table ---
  if (Table) {
    cat("\nIteration table:\n")
    print(x$Table, digits = digits)
  }

  invisible(x)
}


summary.GPArotation <- function(object, digits = 3L, Structure = TRUE, ...) {

  object <- .sortGPALoadings(object)

  # Preserve covariance attribute through the copy
  L <- object$loadings
  attr(L, "covariance") <- attr(object$loadings, "covariance")

  r <- list(
    loadings      = L,
    Phi           = object$Phi,
    method        = object$method,
    orthogonal    = object$orthogonal,
    convergence   = object$convergence,
    iters         = rev(object$Table[, 1])[1],
    Structure     = Structure,
    digits        = digits,
    algorithm     = object$algorithm,
    fwindow       = object$fwindow,
    randStartChar = object$randStartChar   # NULL if not present, harmless
  )
  class(r) <- "summary.GPArotation"
  r
}

print.summary.GPArotation <- function(x, cutoff = 0.1, ...) {

  # --- Convergence header ---
  cat(if (x$orthogonal) "Orthogonal" else "Oblique")
  cat(" rotation method", x$method)
  if (!x$convergence) cat(" NOT")
  cat(" converged in ", x$iters, " iterations.\n", sep = "")

 # Always show in summary if available
  if (!is.null(x$algorithm) && !is.null(x$fwindow)) {
   cat("Algorithm:", x$algorithm, " |  fwindow:", x$fwindow, "\n")
  ##    " |  Iterations:", nrow(x$Table) - 1, "\n")
  }
     
  # --- Random start diagnostics ---
  if (!is.null(x$randStartChar)) {
    rs <- x$randStartChar
    cat("Of ", rs[1], " random starts ",
        round(100 * rs[2] / rs[1]), "% converged, ",
        round(100 * rs[3] / rs[1]), "% at the same lowest minimum.\n", sep = "")
    if (rs[4] > 1)
      cat("Random starts converged to ", rs[4], " different local minima.\n", sep = "")
  }

  # --- AUC and FSI --- computed after sorting, before printing ---
  auc <- calc_AUC(x$loadings)
  fsi <- calc_FSI(x$loadings)

  # --- Pattern matrix header ---
  prstr <- !x$orthogonal && x$Structure
  if (prstr) cat("Pattern (loadings) - see Structure (correlations) below:\n")

  # --- Loadings + varex block (same approach as print.GPArotation) ---
  L         <- x$loadings
  L_rounded <- round(L, x$digits)
  Lambda    <- unclass(L_rounded)
  p         <- nrow(Lambda)

  cat("\nLoadings:\n")
  fx <- setNames(format(round(Lambda, x$digits)), NULL)
  nc <- nchar(fx[1L], type = "c")
  fx[abs(Lambda) < cutoff] <- strrep(" ", nc)
  print(fx, quote = FALSE)

  # Build varex block
  Phi_cov <- attr(L, "covariance")
  vx <- if (is.null(Phi_cov) || identical(Phi_cov, diag(ncol(Lambda))))
    colSums(Lambda^2)
  else
    diag(Phi_cov %*% crossprod(Lambda))
  names(vx) <- colnames(Lambda)

  varex <- rbind(`SS loadings` = vx)

  if (x$orthogonal) {
    varex <- rbind(varex, `Proportion Var` = vx / p)
    if (ncol(Lambda) > 1)
      varex <- rbind(varex, `Cumulative Var` = cumsum(vx / p))
  }

  varex <- rbind(varex, AUC = auc$AUC, FSI = fsi$FSI)

  cat("\n")
  print(round(varex, x$digits))

  # --- Structure matrix (oblique only) ---
  if (prstr) {
    S <- x$loadings %*% x$Phi
    dimnames(S) <- dimnames(x$loadings)
    S_rounded <- round(S, x$digits)
    class(S_rounded)              <- "loadings"
    attr(S_rounded, "covariance") <- diag(ncol(S))

    out <- capture.output(print(S_rounded, cutoff = cutoff, digits = x$digits))
    out[grep("^Loadings", out)] <- "Correlations:"
    out <- out[!grepl("^SS loadings", out)]
    out <- out[!grepl("^Proportion",  out)]
    out <- out[!grepl("^Cumulative",  out)]
    while (length(out) > 0 && grepl("^\\s*(Factor\\d+\\s*)+$", trimws(out[length(out)])))
      out <- out[-length(out)]
    while (length(out) > 0 && nchar(trimws(out[length(out)])) == 0)
      out <- out[-length(out)]
    cat(paste(out, collapse = "\n"), "\n")
  }

  # --- Simple structure summary ---
  ssi <- calc_simplicity(x$loadings)
  hp  <- calc_hyperplane(x$loadings)

  #cat("\nMean AUC:", auc$AUC_mean, " Min AUC:", auc$AUC_min, "\n")
  #cat("Mean FSI:", fsi$FSI_mean, " Min FSI:", fsi$FSI_min, "\n")
  cat("\nHyperplane total:", hp$HP_total, "of", hp$HP_max,
      "(", hp$HP_pct, "%) at cutoff", hp$cutoff, "\n")

  cat("\nPost-Hoc Simplicity Suite (overall solution):\n")
  cat("Hoffman Index:   ", format(round(ssi$Hoffman, x$digits), nsmall = x$digits), "\n")
  cat("Gini Coeficient: ", format(round(ssi$Gini,    x$digits), nsmall = x$digits), "\n")
  cat("Bentler Index:   ", format(round(ssi$Bentler, x$digits), nsmall = x$digits), "\n")

  # --- Phi (oblique only) ---
  if (!x$orthogonal) {
    cat("\nPhi:\n")
    Phi <- x$Phi
    dimnames(Phi) <- list(colnames(x$loadings), colnames(x$loadings))
    print(round(Phi, x$digits))
  }

  invisible(x)
}

###############################################################
### S3 functions
###############################################################

residuals.GPArotation <- function(object, ...) {
  R <- object$correlation
  if (is.null(R))
    stop("'correlation' not stored. Rotate a factanal object with n.obs specified.",
         call. = FALSE)
  L   <- unclass(object$loadings)
  Phi <- if (is.null(object$Phi)) diag(ncol(L)) else object$Phi
  R_hat       <- L %*% Phi %*% t(L)
  diag(R_hat) <- 1
  Resid       <- R - R_hat
  diag(Resid) <- 0
  Resid
}

