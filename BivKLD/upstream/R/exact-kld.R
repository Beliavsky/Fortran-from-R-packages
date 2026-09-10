.positive_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(sprintf("'%s' must be one finite positive number.", name),
             call. = FALSE)
    }
    as.numeric(x)
}

.probability_array <- function(x, name, normalize, tolerance) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0) || length(x) == 0L) {
        stop(sprintf("'%s' must contain finite non-negative probabilities.",
                     name), call. = FALSE)
    }
    total <- sum(x)
    if (total <= 0) {
        stop(sprintf("'%s' must have a positive sum.", name), call. = FALSE)
    }
    if (normalize) {
        x <- x / total
    } else if (abs(total - 1) > tolerance) {
        stop(sprintf("'%s' must sum to one (or use normalize = TRUE).", name),
             call. = FALSE)
    }
    x
}

#' Exact Bivariate KL Divergence for Discrete Distributions
#'
#' Computes the directed KL divergence between two probability tables. The
#' inputs can be matrices, arrays, or vectors, but must have identical shapes.
#'
#' @param p,q Probability tables for the first and second distributions.
#' @param normalize If `TRUE`, normalize each input to sum to one.
#' @param tolerance Allowed absolute error in a probability sum.
#'
#' @return A non-negative numeric scalar, possibly `Inf` when `q` is zero at a
#'   cell to which `p` assigns positive probability.
#'
#' @examples
#' p <- matrix(c(0.2, 0.3, 0.1, 0.4), 2)
#' q <- matrix(c(0.1, 0.4, 0.2, 0.3), 2)
#' biv_kld_discrete(p, q)
#'
#' @export
biv_kld_discrete <- function(p, q, normalize = FALSE,
                             tolerance = sqrt(.Machine$double.eps)) {
    if (!identical(dim(p), dim(q)) || length(p) != length(q)) {
        stop("'p' and 'q' must have identical dimensions and lengths.",
             call. = FALSE)
    }
    p <- .probability_array(p, "p", normalize, tolerance)
    q <- .probability_array(q, "q", normalize, tolerance)
    positive <- p > 0
    if (any(positive & q == 0)) {
        return(Inf)
    }
    sum(p[positive] * (log(p[positive]) - log(q[positive])))
}

.covariance_matrix <- function(x, name) {
    if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(2L, 2L)) ||
        any(!is.finite(x)) || max(abs(x - t(x))) > sqrt(.Machine$double.eps) ||
        any(eigen(x, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        stop(sprintf("'%s' must be a symmetric positive-definite 2 by 2 matrix.",
                     name), call. = FALSE)
    }
    x
}

#' Exact KL Divergence Between Bivariate Normal Distributions
#'
#' Computes `D(N(mean1, sigma1) || N(mean2, sigma2))` analytically.
#'
#' @param mean1,mean2 Numeric vectors of length two.
#' @param sigma1,sigma2 Symmetric positive-definite two-by-two covariance
#'   matrices.
#'
#' @return A non-negative numeric scalar.
#'
#' @examples
#' biv_kld_normal(c(-2, 2), matrix(c(1, 1, 1, 2), 2),
#'                c(-2, 2), matrix(c(3, 3, 3, 6), 2))
#'
#' @export
biv_kld_normal <- function(mean1, sigma1, mean2, sigma2) {
    valid_mean <- function(z, name) {
        if (!is.numeric(z) || length(z) != 2L || any(!is.finite(z))) {
            stop(sprintf("'%s' must contain two finite numbers.", name),
                 call. = FALSE)
        }
        as.numeric(z)
    }
    mean1 <- valid_mean(mean1, "mean1")
    mean2 <- valid_mean(mean2, "mean2")
    sigma1 <- .covariance_matrix(sigma1, "sigma1")
    sigma2 <- .covariance_matrix(sigma2, "sigma2")
    delta <- mean2 - mean1
    inverse2 <- solve(sigma2)
    value <- 0.5 * (determinant(sigma2, logarithm = TRUE)$modulus -
                   determinant(sigma1, logarithm = TRUE)$modulus - 2 +
                   sum(diag(inverse2 %*% sigma1)) +
                   drop(t(delta) %*% inverse2 %*% delta))
    as.numeric(value)
}

#' Exact KL Divergence for Two Models from the Reference Paper
#'
#' `biv_kld_pareto2()` evaluates the divergence between bivariate Pareto type
#' II models with survival functions `(1 + x1 + x2)^(-alpha)` and
#' `(1 + x1 + x2)^(-beta)`. `biv_kld_independent_weibull()` evaluates the
#' divergence between models with survival functions
#' `exp(-alpha[1] * x1 - alpha[2] * x2)` and the corresponding expression in
#' `beta`.
#'
#' @param alpha,beta Positive model parameters. Scalars for
#'   `biv_kld_pareto2()` and vectors of length two for
#'   `biv_kld_independent_weibull()`.
#'
#' @return A non-negative numeric scalar.
#'
#' @references
#' Chackochan, R., Sankaran, P. G., & Unnikrishnan Nair, N. (2026).
#' Bivariate Kullback-Leibler divergence. *Communications in Statistics -
#' Theory and Methods*, 55(1), 292-312.
#' \doi{10.1080/03610926.2025.2496687}
#'
#' @examples
#' biv_kld_pareto2(1, 2)
#' biv_kld_independent_weibull(c(1, 2), c(2, 3))
#'
#' @name model_kld
NULL

#' @rdname model_kld
#' @export
biv_kld_pareto2 <- function(alpha, beta) {
    alpha <- .positive_scalar(alpha, "alpha")
    beta <- .positive_scalar(beta, "beta")
    (beta - alpha) * (1 + 2 * alpha) / (alpha * (1 + alpha)) +
        log(alpha * (1 + alpha) / (beta * (1 + beta)))
}

#' @rdname model_kld
#' @export
biv_kld_independent_weibull <- function(alpha, beta) {
    if (!is.numeric(alpha) || length(alpha) != 2L ||
        any(!is.finite(alpha)) || any(alpha <= 0)) {
        stop("'alpha' must contain two finite positive numbers.", call. = FALSE)
    }
    if (!is.numeric(beta) || length(beta) != 2L ||
        any(!is.finite(beta)) || any(beta <= 0)) {
        stop("'beta' must contain two finite positive numbers.", call. = FALSE)
    }
    sum(beta / alpha) - 2 + log(prod(alpha) / prod(beta))
}

