#' The Matrix Normal Distribution

#' @family distribution
#' @keywords matrixNormal distribution

#' @name matrixNormal_Distribution
NULL
#' @description Computes the density (\code{dmatnorm}), calculates the cumulative distribution function (CDF, \code{pmatnorm}), and generates 1 random number (\code{rmatnorm}) from the matrix normal: \deqn{A \sim MatNorm_{n,p}(M, U, V) }.
#'
#' @details
#' These functions rely heavily on this following property of matrix normal distribution. Let \code{koch()} refer to the Kronecker product of a matrix. For a n x p matrix \emph{A},
#' if  \deqn{A \sim MatNorm(M, U, V),} then
#' \deqn{ vec(A) \sim MVN_{np} (M, Sigma = koch(V,U) ).}
#'
#' Thus, the probability of \code{Lower} < \code{A} < \code{Upper} in the matrix normal can be found by using the CDF of vec(A), which is given by \code{\link[mvtnorm]{pmvnorm}} function in \pkg{mvtnorm}.  See \code{\link[mvtnorm]{algorithms}} and \code{\link[mvtnorm]{pmvnorm}} for more information.
#'
#' Also, we can simulate a random matrix \emph{A} from a matrix normal by sampling \emph{vec(A)} from \code{\link[mvtnorm]{rmvnorm}} function in \pkg{mvtnorm}. This matrix \emph{A} takes the rownames from \emph{U} and the colnames from \emph{V}.

#' @section Calculating Matrix Normal Probabilities:
#'
#'  From the \code{mvtnorm} package, three algorithms are available for evaluating normal probabilities: \itemize{
#'  \item The default is the randomized Quasi-Monte-Carlo procedure by Genz (1992, 1993) and Genz and Bretz (2002) applicable to arbitrary covariance structures and dimensions up to 1000.
#'   \item For smaller dimensions (up to 20) and non-singular covariance matrices, the algorithm by Miwa et al. (2003) can be used as well.
#'   \item For two- and three-dimensional problems and semi-infinite integration region, TVPACK implements an interface to the methods described by Genz (2004).
#'     }
#'  The \code{...} arguments define the hyper-parameters for GenzBertz algorithm:
#'  \describe{
#' \item{maxpts}{maximum number of function values as integer. The internal FORTRAN code always uses a minimum number depending on the dimension.Default 25000.}
#' \item{abseps}{absolute error tolerance.}
#' \item{releps}{relative error tolerance as double.}
#' }

#' @note
#' Ideally, both scale matrices are positive-definite. If they do not appear to be symmetric, the tolerance should be increased. Since symmetry is checked, the `checkSymmetry` arguments in `mvtnorm::rmvnorm()` are set to FALSE.

# @note Not implemented yet.
#  To calculate the determinant more quickly for large matrices, the logdet function of a matrix \emph{A} from LaplacesDemon
#  is used instead of log(det(A)). Author: Statisticat, LLC. software@bayesian-inference.com.
#  #Example: > dmatnorm(A, M, U, V )
#  [1] -6367.045
#  > dmatnorm(A, M, U, V, log = FALSE)
#  [1] 0


#' @references
#'
#' Pocuca, N., Gallaugher, M.P., Clark, K.M., & McNicholas, P.D. (2019). Assessing and Visualizing Matrix Variate Normality. Methodology. <https://arxiv.org/abs/1910.02859>
#'
#' Gupta, A. K. and D. K. Nagar (1999). Matrix Variate Distributions. Boca Raton: Chapman & Hall/CRC Press.

#' @param A The numeric n x p matrix that follows the matrix-normal. Value used to calculate the density.
#' @param M  The mean n x p matrix that is numeric and real. Must contain non-missing values. Parameter of matrix Normal.
#' @param U  The individual scale n x n real positive-definite matrix (rows). Must contain non-missing values. Parameter of matrix Normal.
#' @param V  The parameter scale p x p  real positive-definite matrix (columns). Must contain non-missing values. Parameter of matrix Normal.
#' @inheritParams is.symmetric.matrix
#' @param log Logical; if TRUE, the logarithm of the density is returned.

#' @examples
#' # Data Used
#' # if( !requireNamespace("datasets", quietly = TRUE)) { install.packages("datasets")} #part of baseR.
#' A <- datasets::CO2[1:10, 4:5]
#' M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))
#' V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
#' V # Right covariance matrix (2 x 2), say the covariance between parameters.
#' U <- I(10) # Block of left-covariance matrix ( 84 x 84), say the covariance between subjects.
#'
#' # PDF
#' dmatnorm(A, M, U, V)
#' dmatnorm(A, M, U, V, log = FALSE)
#'
#' # Generating Probability Lower and Upper Bounds (They're matrices )
#' Lower <- matrix(rep(-1, 20), ncol = 2)
#' Upper <- matrix(rep(3, 20), ncol = 2)
#' Lower
#' Upper
#' # The probablity that a randomly chosen matrix A is between Lower and Upper
#' pmatnorm(Lower, Upper, M, U, V)
#'
#' # CDF
#' pmatnorm(Lower = -Inf, Upper, M, U, V)
#' # entire domain = 1
#' pmatnorm(Lower = -Inf, Upper = Inf, M, U, V)
#'
#' # Random generation
#' set.seed(123)
#' M <- cbind(rnorm(3, 435, 296), rnorm(3, 27, 11))
#' U <- diag(1, 3)
#' V <- matrix(c(10, 5, 5, 3), nrow = 2)
#' rmatnorm(1, M, U, V)
#'
#' # M has a different sample size than U; will return an error.
#' \dontrun{
#' M <- cbind(rnorm(4, 435, 296), rnorm(4, 27, 11))
#' rmatnorm(M, U, V)
#' }
#'
#' @rdname matrixNormal_Distribution
#' @importFrom mvtnorm pmvnorm
#' @export dmatnorm
dmatnorm <- function(
  A,
                     M,
                     U,
                     V,
                     tol = .Machine$double.eps^0.5,
                     log = TRUE
                     ) {
  n <- nrow(A)
  p <- ncol(A)

  # Checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (sum(dim(A) == dim(M)) != 2) stop("M must have same dimensions as A.")
  check_matnorm(s = 1, M, U, V, tol)

  # The Log Density
  log.dens <- (-n * p / 2) * log(2 * pi) - p / 2 * log(det(U)) - n / 2 * log(det(V)) +
    -1 / 2 * tr(solve(U) %*% (A - M) %*% solve(V) %*% t(A - M))

  # Return
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}

# #'@importFrom LaplacesDemon logdet() .  This is the Function used:
logdet <- function(x) {
  2 * sum(log(diag(chol(x))))
}
# Uses logdet(U) instead of log(det(U)) which could calculate .
# The log(det(U)) and log(det(V)) terms have same terms but are positive.

# Uses logdet function. Any difference using LaplaceNormal?
dmatnorm.logdet <- function(A, M, U, V,
                            tol = .Machine$double.eps^0.5,
                            log = TRUE) {
  n <- nrow(A)
  p <- ncol(A)

  # Checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (sum(dim(A) == dim(M)) != 2) stop("M must have same dimensions as A.")
  check_matnorm(s = 1, M, U, V, tol)

  # The Log Density
  log.dens <- (-n * p / 2) * log(2 * pi) - p / 2 * logdet(U)
  -n / 2 * logdet(V) +
    -1 / 2 * tr(solve(U) %*% (A - M) %*% solve(V) %*% t(A - M))

  # Return
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}



#' @rdname matrixNormal_Distribution
#' @param Lower	 The n x p matrix of lower limits for CDF.
#' @param Upper	 The n x p matrix of upper limits for CDF.
#' @inheritParams mvtnorm::pmvnorm
#' @export pmatnorm
pmatnorm <- function(
  Lower = -Inf,
                     Upper = Inf,
                     M,
                     U,
                     V,
                     tol = .Machine$double.eps^0.5,
                     keepAttr = TRUE,
                     algorithm = mvtnorm::GenzBretz(),
                     ...
                     ) {
  if (utils::packageVersion("mvtnorm") < "1.1-2") {
    warning("New argument added to `mvtnorm v. 1.1-2`. Please upgrade to avoid error when passing `keepAttr`.")
  }

  n <- nrow(M)
  p <- ncol(M)

  # Checks
  check_matnorm(s = 1, M, U, V, tol)

  # Convert the matrices to lower
  if (is.matrix(Lower)) {
    lower <- vec(Lower)
  } else {
    if (is.vector(Lower) & Lower == -Inf) {
      lower <- -Inf
    } else {
      stop("The lower limit must be a numeric matrix or -Inf.")
    }
  }

  if (is.matrix(Upper)) {
    upper <- vec(Upper)
  } else {
    if (is.vector(Upper) & Upper == Inf) {
      upper <- Inf
    } else {
      stop("The upper limit must be a numeric matrix or Inf.")
    }
  }
  # Calculating the probablity
  prob <- mvtnorm::pmvnorm(
    lower, upper,
    mean = vec(M),
    corr = NULL,
    sigma = kronecker(U, V),
    algorithm = algorithm,
    ...,
    keepAttr = keepAttr
  )
  warning("The covariance matrix is standardized. ")

  return(prob)
}

#' @rdname matrixNormal_Distribution
#' @param s The number of observations desired to simulate from the matrix normal. Defaults to 1. Currently has no effect but acts as a placeholder in future releases.
# #'@inheritParams mvtnorm::rmvnorm
#' @param method String specifying the matrix decomposition used to determine the matrix root of the Kronecker product of U and V in \code{rmatnorm}. Possible methods are eigenvalue decomposition ("eigen"), singular value decomposition ("svd"), and Cholesky decomposition ("chol"). The Cholesky (the default) is typically fastest, but not by much though. Passed to **mvtnorm**::\code{\link[mvtnorm]{rmvnorm}}.

#' @import mvtnorm
#' @export rmatnorm

rmatnorm <- function(
    s = 1,
    M,
    U,
    V,
    tol = .Machine$double.eps^0.5,
    method = "chol"
) {
  if (utils::packageVersion("mvtnorm") < "1.1-2") {
    warning("New argument added to `mvtnorm v. 1.1-2`. Please upgrade to avoid error.")
  }

  # Convert all to matrices -- added 5/2/20
  M <- as.matrix(M)
  U <- as.matrix(U)
  V <- as.matrix(V)

  n <- nrow(M)
  p <- ncol(M)

  # Checks
  # (Symmetry is already checked)
  check_matnorm(s, M, U, V, tol)

  # Vectorizing and sampling from rmvnorm
  if (utils::packageVersion("matrixNormal") <= "0.0.5") {
    warning("The construction of sigma has been found to be incorrect. Please upgrade to new version.")
  }
  # Sigma <- kronecker(U,V) #incorrect -- thanks @prockenschaub, https://github.com/phargarten2/matrixNormal/issues/1
  Sigma <- kronecker(V, U)

  vec.X <- mvtnorm::rmvnorm(1, vec(M), Sigma, method = method, checkSymmetry = FALSE)
  # pre0.9_9994 = FALSE #uses later version of package.

  # Reputting back into a matrix
  X <- matrix(vec.X,
              nrow = n, ncol = p,
              dimnames = list(rownames(U), colnames(V))
  )

  return(X)
}

# cov(vec(A))  #should be 1

# Check to make sure the parameters in MatrixNormal match.
check_matnorm <- function(s, M, U, V, tol) {
  if (!(s > 0)) stop("s must be > 0. s = ", s, call. = FALSE)
  if (anyNA(M)) {
    stop("M contains missing values.", call. = FALSE)
  }
  if (anyNA(U)) {
    stop("U contains missing values.")
  }
  if (anyNA(V)) {
    stop("V contains missing values.")
  }
  if (nrow(M) != nrow(U)) {
    stop("The mean matrix M has different sample size than the scale sample size
         matrix U. M has ", dim(M)[[1]], "rows, and U has ", dim(U)[[1]], ".")
  }
  if (ncol(M) != nrow(V)) {
    stop("The mean matrix M has different number of parameters than scale
         parameter matrix V: M  -- ", dim(M)[2], "; V -- ", dim(V)[1], ".")
  }
  if (!is.positive.definite(U, tol)) {
    stop("U is not positive definite. Calculation may not be accurate.
         Possibly raise tolerance.")
  }
  if (!is.positive.definite(V, tol)) {
    stop("V is not positive definite. Calculation may not be accurate.
         Possibly raise tolerance.")
  }
  return(invisible())
}

# Unsure if should add
# if (!is.matrix(M)) {  M <- matrix(M) }
