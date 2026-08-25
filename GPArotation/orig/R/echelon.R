echelon <- function(A, reference = seq(ncol(A)), ...) {
  # Echelon rotation.
  # Rotates the reference rows to lower triangular (Cholesky) form,
  # producing a unique solution for the factor loading matrix.
  #
  # Args:
  #   A         : Unrotated factor loading matrix (p x k), or a factanal object
  #   reference : Row indices defining the reference variables (default: first k rows)
  #
  # Returns:
  #   A GPArotation object

  mc <- match.call()

  # --- Accept factanal object ---
  extracted   <- .extract_A(A)
  A           <- extracted$A
  n.obs       <- extracted$n.obs
  correlation <- extracted$correlation

  A1 <- A[reference, , drop = FALSE]

  # Cholesky factorization of A1 %*% t(A1)
  # Note: ill-conditioning is a real danger - no singularity check performed
  B1   <- t(chol(A1 %*% t(A1)))
  Tmat <- solve(A1, B1)

  # Assemble rotated solution
  L              <- matrix(0, nrow(A), ncol(A))
  L[reference, ] <- B1
  L[-reference,] <- A[-reference, , drop = FALSE] %*% Tmat
  dimnames(L)    <- list(rownames(A), paste("factor", seq(ncol(A))))

  # --- covariance attribute for print.loadings ---
  class(L) <- "loadings"
  attr(L, "covariance") <- diag(ncol(L))  # orthogonal - plain SS

  r <- list(
    loadings    = L,
    Th          = Tmat,
    method      = "echelon",
    orthogonal  = TRUE,
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

###########################

#echelon <- function(L, reference=seq(NCOL(L)), ...) {

#   # Split L in reference part and the rest
#   A1 <- L[reference,, drop=FALSE]
#   #A2 is L[-reference,]

#   # Compute the part of A Phi A' corresponding to the reference variables
#   # Compute cholesky rot = rotated reference part
#   # No check or error message for singularity. Exact singularity is rare in
#   # practice but ill-conditioning is a real danger.

#   #  now assuming orthogonal (Phi=I)
#   #newPhi <- if (is.null(Phi)) A1 %*% t(A1) else A1 %*% Phi %*% t(A1)
#   #B1 <- t(chol(newPhi))
#   
#   B1 <- t(chol(A1 %*% t(A1)))

#   # Transformation matrix: B1 = A1 * Tmat
#   # Rotated solution for non-reference part: B2 = A2 * Tmat
#   Tmat <- solve(A1, B1)
#   
#   # Assemble rotated solution
#   B <- matrix(0, NROW(L), NCOL(L))
#   B[reference,]  <- B1
#   B[-reference,] <- L[-reference,, drop=FALSE] %*% Tmat

#   dimnames(B) <- list(dimnames(L)[[1]], paste("factor", seq(NCOL(L))))
#   r <- list(loadings=B, Th=Tmat, method="echelon", orthogonal=TRUE, 
#       convergence=TRUE)
#   class(r) <- "GPArotation"
#   r
#}
