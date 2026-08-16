quick_stationary_feasible <- function(x0, lower, upper, ineq_fn = NULL, ineq_lower = NULL, ineq_upper = NULL,
                                      eq_fn = NULL, eq_b = NULL, penalty = 1e4, maxit = 10)
{
    # Inequality constraints (may be empty)
    standard_form_ineq_fn <- function(x) {
        if (is.null(ineq_fn) || is.null(ineq_lower) || is.null(ineq_upper)) return(numeric(0))
        h <- ineq_fn(x)
        c(ineq_lower - h, h - ineq_upper)
    }
    # Equality constraints (may be empty)
    standard_form_eq_fn <- function(x) {
        if (is.null(eq_fn) || is.null(eq_b)) return(numeric(0))
        eq_fn(x) - eq_b
    }
    penalized <- function(x) {
        ineq_val <- standard_form_ineq_fn(x)
        eq_val <- standard_form_eq_fn(x)
        penalty_ineq <- if (length(ineq_val) == 0) 0 else sum(pmax(ineq_val, 0)^2)
        penalty_eq   <- if (length(eq_val)   == 0) 0 else sum(eq_val^2)
        norm_term <- sum((x - x0)^2)
        norm_term + penalty * (penalty_ineq + penalty_eq)
    }
    res <- optim(x0, penalized, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = maxit))
    res$par
}

sample_box_feasible <- function(lower, upper, eps = 1e-4) {
    effective_lower <- lower + eps * (upper > lower + 2*eps)
    effective_upper <- upper - eps * (upper > lower + 2*eps)
    is_too_narrow <- effective_lower >= effective_upper
    res <- numeric(length(lower))
    res[is_too_narrow] <- (lower[is_too_narrow] + upper[is_too_narrow]) / 2
    res[!is_too_narrow] <- effective_lower[!is_too_narrow] +
        (effective_upper[!is_too_narrow] - effective_lower[!is_too_narrow]) * runif(sum(!is_too_narrow))
    res
}

generate_feasible_starts <- function(n, fn = NULL, lower, upper, ineq_fn, ineq_lower, ineq_upper, eq_fn, eq_b, maxit = 100, penalty = 1e4, eps = 1e-4, seed = NULL)
{
    if (!is.null(seed)) set.seed(seed)
    candidates <- vector("list", n)
    for (i in seq_len(n)) {
        x0 <- sample_box_feasible(lower, upper, eps = eps)
        if (!is.null(ineq_fn) | !is.null(eq_fn)) {
            feasible <- quick_stationary_feasible(x0, lower, upper, ineq_fn, ineq_lower, ineq_upper, eq_fn, eq_b, maxit = maxit, penalty = penalty)
        } else {
            feasible <- x0
        }
        candidates[[i]] <- feasible
    }
    # Optionally, rank by objective value
    if (!is.null(fn)) {
        obj_vals <- sapply(candidates, fn)
        ord <- order(obj_vals, decreasing = FALSE)
        candidates <- candidates[ord]
    }
    out <- do.call(rbind, candidates)
    if (!is.null(names(lower))) colnames(out) <- names(lower)
    return(out)
}

