eiv <- function(A, identity = seq(ncol(A)), ...) {
  # Errors-in-variables rotation (eiv).
  # Selects rows specified by 'identity' to define the reference frame,
  # rotating so those rows form an identity matrix.
  #
  # Args:
  #   A        : Unrotated factor loading matrix (p x k), or a factanal object
  #   identity : Row indices defining the identity set (default: first k rows)
  #
  # Returns:
  #   A GPArotation object

  mc <- match.call()

  # --- Accept factanal object ---
  extracted   <- .extract_A(A)
  A           <- extracted$A
  n.obs       <- extracted$n.obs
  correlation <- extracted$correlation

  A1 <- A[identity, , drop = FALSE]
  g  <- solve(A1)

  if (max(abs(diag(nrow(A1)) - A1 %*% g)) > 1e-14)
    warning("Inverse is not well conditioned. ",
            "Consider setting identity to select different rows.")

  L              <- matrix(NA, nrow(A), ncol(A))
  L[identity, ]  <- diag(nrow(A1))
  L[-identity, ] <- A[-identity, , drop = FALSE] %*% g
  dimnames(L)    <- list(rownames(A), paste("factor", seq(ncol(A))))

  # --- covariance attribute for print.loadings ---
  Phi <- tcrossprod(A1)
  class(L) <- "loadings"
  attr(L, "covariance") <- Phi

  r <- list(
    loadings    = L,
    Phi         = Phi,
    Th          = t(A1),
    method      = "eiv",
    orthogonal  = FALSE,
    convergence = TRUE,
    Gq          = NULL,
    Table       = matrix(NA, 1, 4,
                    dimnames = list(NULL, c("iter", "f", "log10(s)", "alpha"))),
    n.obs       = n.obs,
    correlation = correlation,
    call        = mc
  )
  attr(r, "A_unrotated") <- A
  class(r) <- "GPArotation"
  r
}


#########################
#
#eiv <- function(L,  identity=seq(NCOL(L)), ...){
#   A1 <- L[ identity, , drop=FALSE]
#   g <- solve(A1)
#   if(1e-14 < max(abs(diag(1, length(identity)) - A1 %*% g)))
#      warning("Inverse is not well conditioned. Consider setting identity to select different rows.")
#   B <- array(NA, dim(L))
#   B[identity, ] <- diag(1, length(identity))
#   B[-identity,] <- L[-identity,, drop=FALSE] %*% g
#   dimnames(B) <- list(dimnames(L)[[1]], paste("factor", seq(NCOL(L))))
#   r <- list(loadings=B, Th=t(A1), method="eiv", orthogonal=FALSE, convergence=TRUE,
#        Phi= tcrossprod(A1))
#   class(r) <- "GPArotation"
#   r
#  }
