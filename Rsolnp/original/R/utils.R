# auto_scale_factors <- function(pars, fn, gr, min_scale = 1e-12, max_scale = 1e4, eps = .Machine$double.eps^0.25, ...)
# {
#   if (!is.null(gr)) {
#     g <- gr(pars, ...)
#     H_approx <- tcrossprod(g)  # g %*% t(g), outer product
#     if (!corpcor::is.positive.definite(H_approx)) H_approx <- corpcor::make.positive.definite(H_approx)
#     V_approx <- solve(H_approx)
#     diag_V <- abs(diag(V_approx))
#     scale <- sqrt(pmax(diag_V, eps))
#     D <- pmin(pmax(scale, min_scale), max_scale)
#   } else {
#     g <- grad(func = fn, x = pars, method = "Richardson", side = NULL, method.args = list(eps = eps), ...)
#     H_approx <- tcrossprod(g)  # g %*% t(g), outer product
#     if (!corpcor::is.positive.definite(H_approx)) H_approx <- corpcor::make.positive.definite(H_approx)
#     V_approx <- solve(H_approx)
#     diag_V <- abs(diag(V_approx))
#     scale <- sqrt(pmax(diag_V, eps))
#     D <- pmin(pmax(scale, min_scale), max_scale)
#   }
#   return(D)
# }


solnp_problem_setup <- function(pars, fn, gr, eq_fn, eq_b, eq_jac, ineq_fn, ineq_lower, ineq_upper, ineq_jac, lower, upper, ...)
{
  # index of function indicators
  # [1] length of pars
  # [2] has function gradient?
  # [3] 0 [empty for now]
  # [4] has ineq?
  # [5] ineq length
  # [6] has jacobian (inequality)
  # [7] has eq?
  # [8] eq length
  # [9] has jacobian (equality)
  # [10] has upper / lower bounds
  # [11] has either lower/upper bounds or ineq

  setup_index <- rep(0, 11)
  if (any(is.na(pars))) stop("\nNA found in initial pars.")
  n <- length(pars)
  setup_index[1] <- n

  # check objective function (fn)
  initial_fn_value <- try(fn(pars, ...), silent = TRUE)
  if (inherits(initial_fn_value, 'try-error')) stop("\nobjective function (fn) returned an error on evaluation with initial pars.")
  if (length(initial_fn_value) != 1L) stop("\nobjective function (fn) must return a scalar value.")
  if (is.na(initial_fn_value)) stop("\nobjective function (fn) returned NA with initial pars.")

  # check gradient function (gr)
  if (!is.null(gr)) {
    initial_gr_value <- try(gr(pars, ...), silent = TRUE)
    if (inherits(initial_gr_value, 'try-error')) stop("\ngradient (gr) returned an error on evaluation with initial pars.")
    if (length(initial_gr_value) != n) stop("\ngradient (gr) must return a vector equal to length(pars).")
    if (any(is.na(initial_gr_value))) stop("\nNA's detected in gradient (gr) evaluation with initial pars.")
    if (any(!is.finite(initial_gr_value))) stop("\nnon finite values detected in gradient (gr) evaluation with initial pars.")
    setup_index[2] <- 1
  }

  # check inequality function (ineq_fn) and jacobian
  if (!is.null(ineq_fn)) {
    initial_ineq_value <- try(ineq_fn(pars, ...), silent = TRUE)
    if (inherits(initial_ineq_value, 'try-error')) stop("\nineq_fn returned an error on evaluation with initial pars.")
    nineq <- length(initial_ineq_value)
    setup_index[4] <- 1
    setup_index[5] <- nineq
    if (is.null(ineq_lower) | is.null(ineq_upper)) stop("\nineq_lower and ineq_upper cannot be NULL when ineq_fn provided.")
    if (length(ineq_lower) != nineq) stop(paste0("\nineq_lower must be of length : ", nineq))
    if (length(ineq_upper) != nineq) stop(paste0("\nineq_upper must be of length : ", nineq))
    if (any(is.na(ineq_upper)) | any(!is.finite(ineq_upper))) stop("\nineq_upper must not contain any NA or non finite values.")
    if (any(is.na(ineq_lower)) | any(!is.finite(ineq_lower))) stop("\nineq_lower must not contain any NA or non finite values.")
    if (any((ineq_upper - initial_ineq_value) < 0)) warning("\nupper inequality values violated with initial pars.")
    if (any((initial_ineq_value - ineq_lower) < 0)) warning("\nlower inequality values violated with initial pars.")
    if (!is.null(ineq_jac)) {
      initial_ineq_jac_value <- try(ineq_jac(pars, ...), silent = TRUE)
      if (inherits(initial_ineq_jac_value, 'try-error')) stop("\nineq_jac returned an error on evaluation with initial pars.")
      if (!is.matrix(initial_ineq_jac_value)) stop(paste0("\nineq_jac must return a matrix of dimensions : ", nineq, " x ", n))
      if (any(is.na(initial_ineq_jac_value))) stop("\nNA's detected in ineq_jac evaluation with initial pars.")
      if (any(!is.finite(initial_ineq_jac_value))) stop("\nnon-finite values detected in ineq_jac evaluation with initial pars.")
      setup_index[6] <- 1
    }
  }

  # check equality function (eq_fn) and jacobian
  if (!is.null(eq_fn)) {
    initial_eq_value <- try(eq_fn(pars, ...), silent = TRUE)
    if (inherits(initial_eq_value, 'try-error')) stop("\neq_fn returned an error on evaluation with initial pars.")
    neq <- length(initial_eq_value)
    setup_index[7] <- 1
    setup_index[8] <- neq
    if (is.null(eq_b)) stop("\neq_b cannot be NULL with eq_fn provided.")
    if (any(is.na(eq_b)) | any(!is.finite(eq_b))) stop("\neq_b must not contain any NA or non finite values.")
    if (!is.null(eq_jac)) {
      initial_eq_jac_value <- try(eq_jac(pars, ...), silent = TRUE)
      if (inherits(initial_eq_jac_value, 'try-error')) stop("\neq_jac returned an error on evaluation with initial pars.")
      if (!is.matrix(initial_eq_jac_value)) stop(paste0("\neq_jac must return a matrix of dimensions : ", neq, " x ", n))
      if (any(is.na(initial_eq_jac_value))) stop("\nNA's detected in eq_jac evaluation with initial pars.")
      if (any(!is.finite(initial_eq_jac_value))) stop("\nnon-finite values detected in eq_jac evaluation with initial pars.")
      setup_index[9] <- 1
    }
  }

  # check bounds

  if (!is.null(lower) || !is.null(upper)) {
    setup_index[10] <- 1
    if (length(lower) != n) stop(paste0("\nlower must be of length : ", n))
    if (length(upper) != n) stop(paste0("\nupper must be of length : ", n))
    if (any(is.na(lower)) | any(!is.finite(lower))) stop("\nlower must not contain any NA or non finite values.")
    if (any(is.na(upper)) | any(!is.finite(upper))) stop("\nupper must not contain any NA or non finite values.")
    if (any((upper - pars) < 0)) warning("\nupper bound values violated with initial pars.")
    if (any((pars - lower) < 0)) warning("\nlower bound values violated with initial pars.")
  }
  setup_index[11] <- max(setup_index[10], setup_index[4])
  return(setup_index)
}


solnp_objective_wrapper <- function(fn, gr, n, ...)
{
  f <- function(pars)
  {
    if (any(is.na(pars))) stop("\nNA's detected in parameters")
    out <- try(fn(pars, ...), silent = TRUE)
    if (inherits(out, 'try-error')) {
      warning("\nerror in objective function. Returning 1e10")
      return(1e10)
    } else {
      if (is.na(out) | !is.finite(out)) return(1e10)
      return(out)
    }
  }

  if (!is.null(gr)) {
    g <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      out <- try(gr(pars, ...), silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in gradient function. Returning 1e10")
        return(1e10)
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) return(rep(1e5, n))
        return(out)
      }
    }
  } else {
    g <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      out <- try(grad(func = f, x = pars), silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in gradient function. Returning 1e10")
        return(1e10)
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) return(rep(1e5, n))
        return(out)
      }
    }
  }
  return(list(objective_fun = f, gradient_fun = g))
}


# solnp_scaled_objective_wrapper <- function(fn, gr, n, scale_vector, ...)
# {
#   f <- function(pars)
#   {
#     if (any(is.na(pars))) stop("\nNA's detected in parameters")
#     out <- try(fn(pars * scale_vector, ...), silent = TRUE)
#     if (inherits(out, 'try-error')) {
#       warning("\nerror in objective function. Returning 1e10")
#       return(1e10)
#     } else {
#       if (is.na(out) | !is.finite(out)) return(1e10)
#       return(out)
#     }
#   }
#
#   if (!is.null(gr)) {
#     g <- function(pars)
#     {
#       if (any(is.na(pars))) stop("\nNA's detected in parameters")
#       out <- try(gr(pars * scale_vector, ...) * scale_vector, silent = TRUE)
#       if (inherits(out, 'try-error')) {
#         warning("\nerror in gradient function. Returning 1e10")
#         return(1e10)
#       } else {
#         if (any(is.na(out)) | any(!is.finite(out))) return(rep(1e5, n))
#         return(out)
#       }
#     }
#   } else {
#     g <- function(pars)
#     {
#       if (any(is.na(pars))) stop("\nNA's detected in parameters")
#       out <- try(grad(func = f, x = pars * scaled_vector) * scale_vector, silent = TRUE)
#       if (inherits(out, 'try-error')) {
#         warning("\nerror in gradient function. Returning 1e10")
#         return(1e10)
#       } else {
#         if (any(is.na(out)) | any(!is.finite(out))) return(rep(1e5, n))
#         return(out)
#       }
#     }
#   }
#   return(list(objective_fun = f, gradient_fun = g))
# }


solnp_ineq_wrapper <- function(ineq_fn, ineq_jac, n_ineq, n, ...)
{
  f <- j <- NULL
  if (!is.null(ineq_fn)) {
    f <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      # Apply scaling factor if needed
      out <- try(ineq_fn(pars, ...), silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in inequality function evaluation. Returning large violation.")
        return(rep(1e10, n_ineq))
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) {
          out[is.na(out)] <- 1e10
          out[!is.finite(out)] <- 1e10
        }
        return(out)
      }
    }

    j <- NULL
    if (!is.null(ineq_jac)) {
      # If user provides an analytical Jacobian
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to inequality Jacobian function wrapper.")
        }
        jac_out <- try(ineq_jac(pars, ...), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in user-provided analytical inequality Jacobian function. Returning zeros.")
          return(.zeros(n_ineq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in user-provided analytical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    } else {
      # If no analytical Jacobian is provided, use numerical differentiation (numDeriv)
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to numerical inequality Jacobian function wrapper.")
        }
        jac_out <- try(jacobian(func = f, x = pars, method = "Richardson"), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in numerical inequality Jacobian calculation. Returning zeros.")
          return(.zeros(n_ineq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in numerical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    }
  } else {
    f <- function(pars) return(NULL)
    j <- function(pars) return(NULL)
  }
  return(list(ineq_f = f, ineq_j = j))
}

solnp_scaled_ineq_wrapper <- function(ineq_fn, ineq_jac, n_ineq, n, scale_vector, ...)
{
  f <- j <- NULL
  if (!is.null(ineq_fn)) {
    f <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      # Apply scaling factor if needed
      out <- try(ineq_fn(pars * scale_vector, ...), silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in inequality function evaluation. Returning large violation.")
        return(rep(1e10, n_ineq))
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) {
          out[is.na(out)] <- 1e10
          out[!is.finite(out)] <- 1e10
        }
        return(out)
      }
    }

    j <- NULL
    if (!is.null(ineq_jac)) {
      # If user provides an analytical Jacobian
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to inequality Jacobian function wrapper.")
        }
        jac_out <- try(ineq_jac(pars * scale_vector, ...) %*% diag(scale_vector), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in user-provided analytical inequality Jacobian function. Returning zeros.")
          return(.zeros(n_ineq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in user-provided analytical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    } else {
      # If no analytical Jacobian is provided, use numerical differentiation (numDeriv)
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to numerical inequality Jacobian function wrapper.")
        }
        jac_out <- try(jacobian(func = f, x = pars * scale_vector, method = "Richardson") %*% diag(scale_vector), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in numerical inequality Jacobian calculation. Returning zeros.")
          return(.zeros(n_ineq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in numerical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    }
  } else {
    f <- function(pars) return(NULL)
    j <- function(pars) return(NULL)
  }
  return(list(ineq_f = f, ineq_j = j))
}


solnp_eq_wrapper <- function(eq_fn, eq_b, eq_jac, n_eq, n , ...)
{
  f <- j <- NULL
  if (!is.null(eq_fn)) {
    f <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      # Apply scaling factor if needed
      out <- try(eq_fn(pars, ...) - eq_b, silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in equality function evaluation. Returning large violation.")
        return(rep(1e10, n_eq))
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) {
          out[is.na(out)] <- 1e10
          out[!is.finite(out)] <- 1e10
        }
        return(out)
      }
    }
    j <- NULL
    if (!is.null(eq_jac)) {
      # If user provides an analytical Jacobian
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to equality Jacobian function wrapper.")
        }
        jac_out <- try(eq_jac(pars, ...), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in user-provided analytical inequality Jacobian function. Returning zeros.")
          return(.zeros(n_eq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in user-provided analytical equality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    } else {
      # If no analytical Jacobian is provided, use numerical differentiation (numDeriv)
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to numerical inequality Jacobian function wrapper.")
        }
        jac_out <- try(jacobian(func = f, x = pars, method = "Richardson"), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in numerical inequality Jacobian calculation. Returning zeros.")
          return(.zeros(n_eq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in numerical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    }
  }  else {
    f <- function(pars) return(NULL)
    j <- function(pars) return(NULL)
  }
  return(list(eq_f = f, eq_j = j))
}

solnp_scaled_eq_wrapper <- function(eq_fn, eq_b, eq_jac, n_eq, n, scale_vector, ...)
{
  f <- j <- NULL
  if (!is.null(eq_fn)) {
    f <- function(pars)
    {
      if (any(is.na(pars))) stop("\nNA's detected in parameters")
      # Apply scaling factor if needed
      out <- try(eq_fn(pars * scale_vector, ...) - eq_b, silent = TRUE)
      if (inherits(out, 'try-error')) {
        warning("\nerror in equality function evaluation. Returning large violation.")
        return(rep(1e10, n_eq))
      } else {
        if (any(is.na(out)) | any(!is.finite(out))) {
          out[is.na(out)] <- 1e10
          out[!is.finite(out)] <- 1e10
        }
        return(out)
      }
    }
    j <- NULL
    if (!is.null(eq_jac)) {
      # If user provides an analytical Jacobian
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to equality Jacobian function wrapper.")
        }
        jac_out <- try(eq_jac(pars * scale_vector, ...) %*% diag(scale_vector), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in user-provided analytical inequality Jacobian function. Returning zeros.")
          return(.zeros(n_eq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in user-provided analytical equality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    } else {
      # If no analytical Jacobian is provided, use numerical differentiation (numDeriv)
      j <- function(pars) {
        if (any(is.na(pars))) {
          stop("\nNA's detected in parameters passed to numerical inequality Jacobian function wrapper.")
        }
        jac_out <- try(jacobian(func = f, x = pars * scale_vector, method = "Richardson") %*% diag(scale_vector), silent = TRUE)
        if (inherits(jac_out, 'try-error')) {
          warning("\nError in numerical inequality Jacobian calculation. Returning zeros.")
          return(.zeros(n_eq, n)) # Return a matrix of zeros on error
        } else {
          if (any(is.na(jac_out)) | any(!is.finite(jac_out))) {
            warning("\nNA/Inf detected in numerical inequality Jacobian. Converting to zeros.")
            jac_out[is.na(jac_out)] <- 0
            jac_out[!is.finite(jac_out)] <- 0
          }
          return(jac_out)
        }
      }
    }
  }  else {
    f <- function(pars) return(NULL)
    j <- function(pars) return(NULL)
  }
  return(list(eq_f = f, eq_j = j))
}


augmented_function_gradient <- function(gr, n, n_ineq)
{
  f <- function(x) {
    # Extract original parameters from the augmented vector (augmented_pars = [s, P])
    original_parameters <- x[(n_ineq + 1):(n_ineq + n)]
    # Call the user's gradient function
    user_grad <- gr(original_parameters)
    # Construct the full augmented gradient: [ 0 | d(obj)/d(P) ]
    # Objective function does not depend on slack variables, so components for s are zero.
    augmented_grad <- c(.zeros(n_ineq, 1), user_grad)
    return(augmented_grad)
  }
  return(f)
}

augmented_equality_jacobian <- function(eq_jac, n_eq, n_ineq, n)
{
  f <- function(x) {
    original_parameters <- x[(n_ineq + 1):(n_ineq + n)]
    user_eq_jac <- eq_jac(original_parameters)
    # Construct the full augmented Jacobian: [ 0 | d(eq)/d(P) ]
    # Equality constraints do not depend on slack variables, so columns for s are zero.
    augmented_jac <- cbind(.zeros(n_eq, n_ineq), user_eq_jac)
    return(augmented_jac)
  }
  return(f)
}

augmented_inequality_jacobian <- function(ineq_jac, n_ineq, n)
{
  f <- function(x) {
    original_parameters <- x[(n_ineq + 1):(n_ineq + n)]
    user_ineq_jac <- ineq_jac(original_parameters)
    # Construct the full augmented Jacobian: [ d(ineq)/d(P) | d(ineq)/d(s) ]
    # where d(ineq)/d(s) is -I because the internal constraint is G_j(P,s) = ineq_j(P) - s_j = 0
    augmented_jac <- cbind(user_ineq_jac, -diag(n_ineq))
    return(augmented_jac)
  }
  return(f)
}


.solnp_ctrl <- function(control_list) {
  # Set default control parameters
  default_control <- list(
    rho = 1,        # Penalty parameter
    max_iter = 400,     # Maximum number of major iterations
    min_iter = 800,     # Maximum number of minor iterations
    tol = 1.0e-8,   # Tolerance on feasibility and optimality
    trace = FALSE   # Trace optimization progress
  )
  # Update defaults with user-provided control parameters
  for (name in names(control_list)) {
    if (name %in% names(default_control)) {
      default_control[[name]] <- control_list[[name]]
    } else {
      warning(paste("Unknown control parameter:", name))
    }
  }
  return(default_control)
}

.solver_warnings <- function(message_indicator)
{
  m1 <- paste0("\nsolnp: Redundant constraints were found.")
  m2 <- paste0("\nLinearized problem has no feasible solution.The problem may not be feasible.")
  m3 <- paste0("\nMinor optimization routine did not converge in the specified number of minor iterations.")
  ans <- switch(message_indicator,
               "M1" = m1,
               "M2" = m2,
               "M3" = m3)
  warning(ans)
}

