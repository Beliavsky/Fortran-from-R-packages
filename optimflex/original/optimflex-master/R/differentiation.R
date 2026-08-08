#' @title Internal Epsilon Calculation for Finite Differences
#' @description Calculates an optimal step size (h) for numerical differentiation.
#' @param x_i Numeric. The parameter value.
#' @param type String. Type of derivative (grad_fwd, grad_cnt, hess_fwd, hess_cnt).
#' @noRd
get_eps <- function(x_i, type = c("grad_fwd", "grad_cnt", "hess_fwd", "hess_cnt")) {
  m_eps <- .Machine$double.eps
  scale <- max(abs(x_i), 1.0)
  
  h <- switch(match.arg(type),
              "grad_fwd" = sqrt(m_eps) * scale,
              "grad_cnt" = (m_eps)^(1/3) * scale,
              "hess_fwd" = (m_eps)^(1/4) * scale,
              "hess_cnt" = (m_eps)^(1/4) * scale)
  return(h)
}

#' Fast Numerical Gradient
#'
#' @description
#' Provides a high-speed numerical gradient using forward or central differences.
#'
#' @param f Function. The objective function.
#' @param x Numeric vector. Parameters at which to evaluate the gradient.
#' @param diff_method String. Differentiation method: "forward" or "central".
#' @param ... Additional arguments passed to \code{f}.
#'
#' @return A numeric vector of gradients.
#' @export
fast_grad <- function(f, x, diff_method = c("forward", "central"), ...) {
  diff_method <- match.arg(diff_method)
  n <- length(x)
  g <- numeric(n)
  
  if (diff_method == "forward") {
    f0 <- Re(as.numeric(f(x, ...))[1])
    for (i in seq_len(n)) {
      h <- get_eps(x[i], "grad_fwd")
      x_up <- replace(x, i, x[i] + h)
      g[i] <- (Re(as.numeric(f(x_up, ...))[1]) - f0) / h
    }
  } else {
    for (i in seq_len(n)) {
      h <- get_eps(x[i], "grad_cnt")
      x_u <- replace(x, i, x[i] + h)
      x_d <- replace(x, i, x[i] - h)
      g[i] <- (Re(as.numeric(f(x_u, ...))[1]) - Re(as.numeric(f(x_d, ...))[1])) / (2 * h)
    }
  }
  return(g)
}

#' Fast Numerical Hessian
#'
#' @description
#' High-speed numerical Hessian calculation using finite differences.
#'
#' @param f Function. The objective function.
#' @param x Numeric vector. Parameters at which to evaluate the Hessian.
#' @param diff_method String. Differentiation method: "forward" or "central".
#' @param ... Additional arguments passed to \code{f}.
#'
#' @return A symmetric Hessian matrix.
#' @export
fast_hess <- function(f, x, diff_method = c("forward", "central"), ...) {
  diff_method <- match.arg(diff_method)
  n <- length(x)
  H <- matrix(0, n, n)
  f0 <- Re(as.numeric(f(x, ...))[1])
  
  if (diff_method == "forward") {
    h_vec <- sapply(x, get_eps, type = "hess_fwd")
    f_grads <- numeric(n)
    for (i in seq_len(n)) f_grads[i] <- Re(as.numeric(f(replace(x, i, x[i] + h_vec[i]), ...))[1])
    
    for (i in seq_len(n)) {
      for (j in i:n) {
        if (i == j) {
          x_ii <- replace(x, i, x[i] + 2 * h_vec[i])
          H[i, i] <- (Re(as.numeric(f(x_ii, ...))[1]) - 2 * f_grads[i] + f0) / (h_vec[i]^2)
        } else {
          x_ij <- replace(replace(x, i, x[i] + h_vec[i]), j, x[j] + h_vec[j])
          val <- (Re(as.numeric(f(x_ij, ...))[1]) - f_grads[i] - f_grads[j] + f0) / (h_vec[i] * h_vec[j])
          H[i, j] <- H[j, i] <- val
        }
      }
    }
  } else {
    h_vec <- sapply(x, get_eps, type = "hess_cnt")
    for (i in seq_len(n)) {
      for (j in i:n) {
        if (i == j) {
          f_u <- Re(as.numeric(f(replace(x, i, x[i] + h_vec[i]), ...))[1])
          f_d <- Re(as.numeric(f(replace(x, i, x[i] - h_vec[i]), ...))[1])
          H[i, i] <- (f_u - 2 * f0 + f_d) / (h_vec[i]^2)
        } else {
          f_pp <- Re(as.numeric(f(replace(replace(x, i, x[i] + h_vec[i]), j, x[j] + h_vec[j]), ...))[1])
          f_pn <- Re(as.numeric(f(replace(replace(x, i, x[i] + h_vec[i]), j, x[j] - h_vec[j]), ...))[1])
          f_np <- Re(as.numeric(f(replace(replace(x, i, x[i] - h_vec[i]), j, x[j] + h_vec[j]), ...))[1])
          f_nn <- Re(as.numeric(f(replace(replace(x, i, x[i] - h_vec[i]), j, x[j] - h_vec[j]), ...))[1])
          H[i, j] <- H[j, i] <- (f_pp - f_pn - f_np + f_nn) / (4 * h_vec[i] * h_vec[j])
        }
      }
    }
  }
  return(H)
}

#' Fast Numerical Jacobian
#'
#' @description
#' Calculates the Jacobian matrix for a vector-valued function (e.g., residuals).
#'
#' @param f_res Function. A function returning a vector of residuals.
#' @param x Numeric vector. Parameters at which to evaluate the Jacobian.
#' @param diff_method String. Differentiation method: "forward" or "central".
#' @param ... Additional arguments passed to \code{f_res}.
#'
#' @return A Jacobian matrix of dimension (m x n).
#' @export
fast_jac <- function(f_res, x, diff_method = c("forward", "central"), ...) {
  diff_method <- match.arg(diff_method)
  r0 <- Re(as.numeric(f_res(x, ...)))
  m <- length(r0)
  n <- length(x)
  J <- matrix(0, m, n)
  
  if (diff_method == "forward") {
    for (j in seq_len(n)) {
      h <- get_eps(x[j], "grad_fwd")
      J[, j] <- (Re(as.numeric(f_res(replace(x, j, x[j] + h), ...))) - r0) / h
    }
  } else {
    for (j in seq_len(n)) {
      h <- get_eps(x[j], "grad_cnt")
      J[, j] <- (Re(as.numeric(f_res(replace(x, j, x[j] + h), ...))) - 
                   Re(as.numeric(f_res(replace(x, j, x[j] - h), ...)))) / (2 * h)
    }
  }
  return(J)
}