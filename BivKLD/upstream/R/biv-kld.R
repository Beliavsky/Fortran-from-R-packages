#' Estimate Bivariate Kullback-Leibler Divergence
#'
#' Estimates the directed divergence from the distribution represented by
#' `x` to that represented by `y`. Bivariate Gaussian kernel densities are
#' evaluated on a common rectangular grid and numerically integrated.
#'
#' @param x,y Numeric matrices or data frames with exactly two columns.
#' @param Hx,Hy Optional symmetric positive-definite bandwidth matrices. When
#'   omitted, bandwidths are selected according to `bandwidth`.
#' @param bandwidth Either `"scv"` for unconstrained smoothed cross-validation
#'   or `"normal"` for the normal-scale selector.
#' @param grid_size One integer or two integers giving the number of grid points
#'   in each coordinate.
#' @param range Optional integration limits. See Details.
#' @param standardize Standardization applied before estimation: `"none"`,
#'   `"pooled"`, or `"separate"`.
#' @param details If `TRUE`, return densities, grid, and bandwidth matrices in
#'   addition to the estimate.
#'
#' @details
#' `range` may be a list of two increasing numeric ranges, a two-by-two matrix
#' whose rows contain coordinate limits, or
#' `c(x_min, x_max, y_min, y_max)`. By default the observed pooled range is
#' padded using the selected bandwidths.
#'
#' Pooled standardization applies the same affine transformation to both
#' samples and therefore leaves the population KL divergence unchanged.
#' Separate standardization, used in the motivating analyses, compares the
#' shapes after removing each sample's own location and scale; it does not
#' estimate the divergence between the original distributions.
#'
#' @return A non-negative numeric estimate when `details = FALSE`. Otherwise,
#'   an object of class `bivkld_estimate` containing the estimate, bandwidths,
#'   evaluation grid, densities, and settings.
#'
#' @references
#' Chackochan, R., Sankaran, P. G., & Unnikrishnan Nair, N. (2026).
#' Bivariate Kullback-Leibler divergence. *Communications in Statistics -
#' Theory and Methods*, 55(1), 292-312.
#' \doi{10.1080/03610926.2025.2496687}
#'
#' @examples
#' set.seed(2026)
#' x <- cbind(rnorm(60), rnorm(60))
#' y <- cbind(rnorm(60, 0.4), rnorm(60, -0.2))
#' biv_kld(x, y, bandwidth = "normal", grid_size = 40)
#'
#' @export
biv_kld <- function(x, y, Hx = NULL, Hy = NULL,
                    bandwidth = c("scv", "normal"), grid_size = 100L,
                    range = NULL,
                    standardize = c("none", "pooled", "separate"),
                    details = FALSE) {
    bandwidth <- match.arg(bandwidth)
    standardize <- match.arg(standardize)
    grid_size <- .validate_grid_size(grid_size)
    samples <- list(.as_bivariate_sample(x, "x"),
                    .as_bivariate_sample(y, "y"))
    samples <- .standardize_samples(samples, standardize)
    bandwidths <- list(
        .select_bandwidth(samples[[1L]], Hx, bandwidth, "Hx"),
        .select_bandwidth(samples[[2L]], Hy, bandwidth, "Hy")
    )
    grid <- .integration_grid(samples, bandwidths, grid_size, range)
    densities <- list(
        x = .density_on_grid(samples[[1L]], bandwidths[[1L]], grid),
        y = .density_on_grid(samples[[2L]], bandwidths[[2L]], grid)
    )
    estimate <- .integrated_kld(densities$x, densities$y, grid$cell_area)
    if (!isTRUE(details)) {
        return(estimate)
    }
    structure(list(
        estimate = estimate,
        bandwidth = list(x = bandwidths[[1L]], y = bandwidths[[2L]]),
        grid = grid,
        density = densities,
        standardize = standardize,
        bandwidth_selector = bandwidth
    ), class = "bivkld_estimate")
}

#' @export
print.bivkld_estimate <- function(x, ...) {
    cat("Bivariate kernel KL divergence\n")
    cat("  D(x || y):", format(x$estimate, digits = 7L), "\n")
    cat("  bandwidth selector:", x$bandwidth_selector, "\n")
    cat("  standardization:", x$standardize, "\n")
    invisible(x)
}

