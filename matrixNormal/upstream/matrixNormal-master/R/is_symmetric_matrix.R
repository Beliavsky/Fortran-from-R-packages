#' Is a matrix symmetric or positive-definite?
#'
#' @name is.symmetric.matrix
#' @family statistics
#' @keywords matrix
#'
#' @description Determines if a matrix is square, symmetric, positive-definite, or positive semi-definite.
#' @details
#' A tolerance is added to indicate if a matrix \emph{A} is approximately symmetric. If \emph{A} is not symmetric, a message and first few rows of the matrix is printed. If \emph{A} has any missing values, NA is returned.
#' \itemize{
#' \item \code{is.symmetric.matrix} returns TRUE if \emph{A} is a numeric, square and symmetric matrix; otherwise, returns FALSE.  A matrix is symmetric if the absolute difference between \emph{A} and its transpose is less than \code{tol}.
#' \item \code{is.positive.semi.definite} returns TRUE if a real, square, and symmetric matrix \emph{A} is positive semi-definite.  A matrix is positive semi-definite if its smallest eigenvalue is greater than or equal to zero.
#' \item \code{is.positive.definite} returns TRUE if a real, square, and symmetric matrix  \emph{A} is positive-definite.  A matrix is positive-definite if its smallest eigenvalue is greater than zero.
#' }

#' @note
#' Functions are adapted from Frederick Novomestky's \pkg{matrixcalc} package in order to implement the \code{rmatnorm} function. The following changes are made: \itemize{
#' \item  I changed argument x to A to reflect usual matrix notation.
#' \item For \code{is.symmetric}, I added a tolerance so that \emph{A} is symmetric even provided small differences between \emph{A} and its transpose. This is useful for \code{rmatnorm} function, which was used repeatedly to generate matrixNormal random variates in a Markov chain.
#' \item For \code{is.positive.semi.definite} and \code{is.positive.definite},  I also saved time by avoiding a \code{for-loop} and instead calculating the minimum of eigenvalues.
#' }


#' @param A A numeric matrix.
#' @param tol A numeric tolerance level used to check if a matrix is symmetric. That is, a matrix is symmetric if the difference between the matrix and its transpose is between -\code{tol} and \code{tol}.

# May want to produce a warning instead of stopping.

#' @examples
#' ## Example 0: Not square matrix
#' B <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#' B
#' is.square.matrix(B)
#'
#' ## Example 1: Not a matrix. should get an error.
#' df <- as.data.frame(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE))
#' df
#' \dontrun{
#' is.square.matrix(df)
#' }
#'
#' ## Example 2: Not symmetric & compare against matrixcalc
#' F <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
#' F
#' is.square.matrix(F)
#' is.symmetric.matrix(F) # should be FALSE
#' if (!requireNamespace("matrixcalc", quietly = TRUE)) {
#'   matrixcalc::is.symmetric.matrix(F)
#' } else {
#'   message("you need to install the package matrixcalc to compare this example")
#' }
#'
#' ## Example 3: Symmetric but negative-definite. The functions are same.
#' # eigenvalues are  3 -1
#' G <- matrix(c(1, 2, 2, 1), nrow = 2, byrow = TRUE)
#' G
#' is.symmetric.matrix(G)
#' if (!requireNamespace("matrixcalc", quietly = TRUE)) {
#'   matrixcalc::is.symmetric.matrix(G)
#' } else {
#'   message("you need to install the package matrixcalc to compare this example.")
#' }
#' isSymmetric.matrix(G)
#' is.positive.definite(G) # FALSE
#' is.positive.semi.definite(G) # FALSE
#'
#' ## Example 3b: A missing value in G
#' G[1, 1] <- NA
#' is.symmetric.matrix(G) # NA
#' is.positive.definite(G) # NA
#'
#' ## Example 4: positive definite matrix
#' # eigenvalues are 3.4142136 2.0000000 0.585786
#' Q <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, byrow = TRUE)
#' is.symmetric.matrix(Q)
#' is.positive.definite(Q)
#'
#' ## Example 5: identity matrix is always positive definite
#' I <- diag(1, 3)
#' is.square.matrix(I) # TRUE
#' is.symmetric.matrix(I) # TRUE
#' is.positive.definite(I) # TRUE
#' @importFrom utils head

#  # Another Symmetric Test found in base. Because of this, is.symmetric() may not be needed
#  isSymmetric.matrix(F)



#' @rdname is.symmetric.matrix
#' @export is.square.matrix
# Is matrix square? (A must be a matrix. ).
# Adapted from matrixcalc::is.square.matrix. Argument name changed to be consistent.
is.square.matrix <- function(A) {
  if (!is.numeric(A)) stop(sprintf("% is not a numeric matrix", A))
  if (!is.matrix(A)) stop(sprintf("% is not a matrix", A))
  is.square <- nrow(A) == ncol(A)
  return(is.square)
}

#' @rdname is.symmetric.matrix
#' @export is.symmetric.matrix
#  @description Is matrix is symmetric? A must be a numeric square matrix with no missing values.
is.symmetric.matrix <- function(A, tol = .Machine$double.eps^0.5) {
  if (anyNA(A)) {
    return(NA)
  }
  if (!is.square.matrix(A)) {
    warning(sprintf("% is not a square matrix", A))
    return(FALSE)
  }

  # Is A and t(A) equal within a tolerance?    # Edited from matrixcalc (PH)
  total.abs <- sum(abs(A - t(A)))
  if (total.abs < tol) { # to avoid being exactly equal
    okay <- TRUE
  } else {
    print("A is not symmetric. Top of the matrix: ")
    print(utils::head(A))
    okay <- FALSE
  }
  # cat("sum( abs(A - t(A)) : ", total.abs, "\n")
  # cat("Total Absolute Difference between Matrix & It's Transpose:", total.abs, "\n")
  return(okay)
}

# Find eigenvalues of a symmetric matrix A with no missing values.
# Paul Hargarten added this function.
find.eval <- function(A, tol = .Machine$double.eps^0.5) {
  # Check if A is symmetric.
  is.symm <- is.symmetric.matrix(A, tol)
  if (is.na(is.symm)) {
    eigenvalues <- NA
  } else if (!is.symm) { # Pass a Negative Number so it is returned false
    # stop(sprintf("%s is not symmetric so no eigenvalues are imputed", A))
    eigenvalues <- -100
  } else if (is.symm) {
    # If A is symmetric, find eigenvalues.
    eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    # Adjust small eigenvalues to be 0  #(Edited from for loop)
    eigenvalues <- ifelse(abs(eigenvalues) < tol, 0, eigenvalues)
  }
  return(eigenvalues)
}

#' @rdname is.symmetric.matrix
#' @export
is.positive.semi.definite <- function(A, tol = .Machine$double.eps^0.5) {
  # Positive semi-definite matrix have non-negative eigenvalues.
  eigenvalues <- find.eval(A, tol)
  if (anyNA(eigenvalues)) {
    return(NA)
  }
  pos.semi <- if (min(eigenvalues) >= 0) TRUE else FALSE
  return(pos.semi)
}

#' @rdname is.symmetric.matrix
#' @export
is.positive.definite <- function(A, tol = .Machine$double.eps^0.5) {
  # Positive definite matrix have positive e-values.
  eigenvalues <- find.eval(A, tol)
  pos <- if (anyNA(eigenvalues)) {
    NA
  } else if (min(eigenvalues) > 0) {
    TRUE
  } else {
    FALSE
  }
  return(pos)
}


# Edited check from tmvtnorm so that the symmetric test includes eigenvalues instead of the determinant.
# But not used in this function.
checkSymmetricPositiveDefinite <- function(x, name = "sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  if (any(diag(x) <= 0)) {
    stop(sprintf(
      "%s all diagonal elements must be positive",
      name
    ))
  }
  min.eval <- min(eigen(x)$values) # edited from tmvnorm function
  if (det(x) <= 0 & min.eval(x) < 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}
