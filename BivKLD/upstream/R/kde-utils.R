.as_bivariate_sample <- function(x, name = deparse(substitute(x))) {
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 2L) {
        stop(sprintf("'%s' must be a numeric matrix or data frame with exactly two columns.",
                     name), call. = FALSE)
    }
    storage.mode(x) <- "double"
    if (nrow(x) < 3L) {
        stop(sprintf("'%s' must contain at least three observations.", name),
             call. = FALSE)
    }
    if (any(!is.finite(x))) {
        stop(sprintf("'%s' must contain only finite values.", name),
             call. = FALSE)
    }
    x
}

.validate_grid_size <- function(grid_size) {
    if (length(grid_size) == 1L) {
        grid_size <- rep(grid_size, 2L)
    }
    if (length(grid_size) != 2L || any(!is.finite(grid_size)) ||
        any(grid_size < 10) || any(grid_size != as.integer(grid_size))) {
        stop("'grid_size' must contain one or two integers of at least 10.",
             call. = FALSE)
    }
    as.integer(grid_size)
}

.validate_bandwidth <- function(H, name) {
    if (!is.matrix(H) || !is.numeric(H) || !identical(dim(H), c(2L, 2L)) ||
        any(!is.finite(H))) {
        stop(sprintf("'%s' must be a finite numeric 2 by 2 matrix.", name),
             call. = FALSE)
    }
    if (max(abs(H - t(H))) > sqrt(.Machine$double.eps) ||
        any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        stop(sprintf("'%s' must be symmetric and positive definite.", name),
             call. = FALSE)
    }
    (H + t(H)) / 2
}

.select_bandwidth <- function(x, H, bandwidth, name) {
    if (!is.null(H)) {
        return(.validate_bandwidth(H, name))
    }
    selector <- switch(bandwidth,
                       scv = ks::Hscv,
                       normal = ks::Hns)
    .validate_bandwidth(selector(x = x), name)
}

.standardize_samples <- function(samples, standardize) {
    if (standardize == "none") {
        return(samples)
    }
    standardize_one <- function(z, center, spread) {
        if (any(!is.finite(spread)) || any(spread <= 0)) {
            stop("Cannot standardize a coordinate with zero variance.",
                 call. = FALSE)
        }
        sweep(sweep(z, 2L, center, "-"), 2L, spread, "/")
    }
    if (standardize == "pooled") {
        pooled <- do.call(rbind, samples)
        center <- colMeans(pooled)
        spread <- apply(pooled, 2L, stats::sd)
        return(lapply(samples, standardize_one,
                      center = center, spread = spread))
    }
    lapply(samples, function(z) {
        standardize_one(z, colMeans(z), apply(z, 2L, stats::sd))
    })
}

.validate_range <- function(range) {
    if (is.list(range) && length(range) == 2L) {
        range <- do.call(rbind, range)
    } else if (is.numeric(range) && length(range) == 4L && is.null(dim(range))) {
        range <- matrix(range, nrow = 2L, byrow = TRUE)
    }
    if (!is.matrix(range) || !is.numeric(range) ||
        !identical(dim(range), c(2L, 2L)) || any(!is.finite(range)) ||
        any(range[, 1L] >= range[, 2L])) {
        stop(paste0("'range' must be a list of two increasing ranges, a 2 by 2 ",
                    "matrix, or c(xmin, xmax, ymin, ymax)."), call. = FALSE)
    }
    range
}

.integration_grid <- function(samples, bandwidths, grid_size, range) {
    if (is.null(range)) {
        pooled <- do.call(rbind, samples)
        endpoints <- t(apply(pooled, 2L, base::range))
        kernel_sd <- vapply(seq_len(2L), function(j) {
            max(vapply(bandwidths, function(H) sqrt(H[j, j]), numeric(1L)))
        }, numeric(1L))
        observed_span <- endpoints[, 2L] - endpoints[, 1L]
        padding <- pmax(0.15 * observed_span, 3 * kernel_sd,
                        sqrt(.Machine$double.eps))
        range <- cbind(endpoints[, 1L] - padding,
                       endpoints[, 2L] + padding)
    } else {
        range <- .validate_range(range)
    }
    axes <- list(seq(range[1L, 1L], range[1L, 2L], length.out = grid_size[1L]),
                 seq(range[2L, 1L], range[2L, 2L], length.out = grid_size[2L]))
    list(points = expand.grid(axes[[1L]], axes[[2L]]),
         axes = axes,
         cell_area = diff(axes[[1L]][1:2]) * diff(axes[[2L]][1:2]),
         range = range)
}

.density_on_grid <- function(x, H, grid) {
    as.numeric(ks::kde(x = x, H = H, eval.points = grid$points)$estimate)
}

.integrated_kld <- function(p, q, cell_area) {
    positive <- p > 0
    if (any(positive & q <= 0)) {
        return(Inf)
    }
    value <- sum(p[positive] * (log(p[positive]) - log(q[positive]))) *
        cell_area
    if (is.finite(value) && value < 0 && abs(value) < 100 * .Machine$double.eps) {
        value <- 0
    }
    value
}
