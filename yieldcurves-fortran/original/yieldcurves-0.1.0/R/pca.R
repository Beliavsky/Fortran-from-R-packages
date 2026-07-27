#' Principal Component Analysis of Yield Curves
#'
#' Perform PCA on a time series of yield curves to extract the dominant
#' factors (level, slope, curvature) following Litterman and Scheinkman
#' (1991).
#'
#' @param curves_matrix Numeric matrix where each row is a yield curve
#'   observation (e.g., daily curves) and each column is a tenor.
#'   Column names should be maturity labels.
#' @param n_components Integer. Number of principal components to retain.
#'   Default is 3 (level, slope, curvature).
#' @param scale Logical. Whether to scale variables before PCA. Default
#'   is `FALSE` (use covariance matrix, standard in yield curve PCA).
#'
#' @return A `yc_pca` object (S3 class) with components:
#'   \describe{
#'     \item{loadings}{Matrix of factor loadings (tenors x components).}
#'     \item{scores}{Matrix of factor scores (observations x components).}
#'     \item{variance_explained}{Numeric vector of proportion of variance
#'       explained by each component.}
#'     \item{cumulative_variance}{Numeric vector of cumulative variance
#'       explained.}
#'     \item{sdev}{Standard deviations of each component.}
#'     \item{n_components}{Number of components retained.}
#'     \item{tenors}{Column names from the input matrix.}
#'   }
#'
#' @references
#' Litterman, R. and Scheinkman, J. (1991). Common Factors Affecting Bond
#' Returns. *The Journal of Fixed Income*, 1(1), 54--61.
#' \doi{10.3905/jfi.1991.692347}
#'
#' @export
#' @examples
#' # Simulate 100 days of yield curves at 5 tenors
#' set.seed(42)
#' n_days <- 100
#' tenors <- c(1, 2, 5, 10, 30)
#' base_rates <- c(0.045, 0.043, 0.042, 0.040, 0.043)
#' curves <- matrix(NA, n_days, length(tenors))
#' colnames(curves) <- paste0(tenors, "Y")
#' level <- cumsum(rnorm(n_days, 0, 0.001))
#' slope <- cumsum(rnorm(n_days, 0, 0.0005))
#' for (i in seq_len(n_days)) {
#'   curves[i, ] <- base_rates + level[i] + slope[i] * (tenors - mean(tenors)) / 30
#' }
#' pca_result <- yc_pca(curves)
#' pca_result
yc_pca <- function(curves_matrix, n_components = 3, scale = FALSE) {
  if (!is.matrix(curves_matrix) && !is.data.frame(curves_matrix)) {
    cli_abort("{.arg curves_matrix} must be a matrix or data frame.")
  }

  curves_matrix <- as.matrix(curves_matrix)

  if (any(is.na(curves_matrix))) {
    cli_abort("{.arg curves_matrix} must not contain NA values.")
  }

  if (nrow(curves_matrix) < 3) {
    cli_abort("{.arg curves_matrix} must have at least 3 rows (observations).")
  }

  if (ncol(curves_matrix) < 2) {
    cli_abort("{.arg curves_matrix} must have at least 2 columns (tenors).")
  }

  n_components <- min(n_components, ncol(curves_matrix), nrow(curves_matrix) - 1)

  pca <- prcomp(curves_matrix, center = TRUE, scale. = scale)

  var_explained <- pca$sdev^2 / sum(pca$sdev^2)

  tenors <- colnames(curves_matrix)
  if (is.null(tenors)) tenors <- paste0("V", seq_len(ncol(curves_matrix)))

  structure(
    list(
      loadings = pca$rotation[, seq_len(n_components), drop = FALSE],
      scores = pca$x[, seq_len(n_components), drop = FALSE],
      variance_explained = var_explained[seq_len(n_components)],
      cumulative_variance = cumsum(var_explained)[seq_len(n_components)],
      sdev = pca$sdev[seq_len(n_components)],
      center = pca$center,
      n_components = n_components,
      tenors = tenors
    ),
    class = "yc_pca"
  )
}

#' Print Method for Yield Curve PCA Objects
#'
#' @param x A `yc_pca` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#'
#' @export
print.yc_pca <- function(x, ...) {
  cli_h1("Yield Curve PCA")
  cli_bullets(c(
    "*" = "Components: {.val {x$n_components}}",
    "*" = "Tenors: {.val {length(x$tenors)}} ({paste(x$tenors, collapse = ', ')})"
  ))
  for (i in seq_len(x$n_components)) {
    labels <- c("Level", "Slope", "Curvature")
    label <- if (i <= 3) labels[i] else paste0("PC", i)
    cli_bullets(c(
      "*" = "PC{i} ({label}): {.val {round(x$variance_explained[i] * 100, 1)}}% variance"
    ))
  }
  cli_bullets(c(
    "*" = "Cumulative: {.val {round(x$cumulative_variance[x$n_components] * 100, 1)}}%"
  ))
  invisible(x)
}

#' Summary Method for Yield Curve PCA Objects
#'
#' @param object A `yc_pca` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.yc_pca <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' Plot Method for Yield Curve PCA Objects
#'
#' Plots the factor loadings for each principal component across tenors.
#'
#' @param x A `yc_pca` object.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return The input object, invisibly.
#'
#' @export
plot.yc_pca <- function(x, ...) {
  n_comp <- x$n_components
  labels <- c("Level", "Slope", "Curvature")
  cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  tenor_idx <- seq_along(x$tenors)

  plot(tenor_idx, x$loadings[, 1], type = "b", pch = 19,
       col = cols[1], ylim = range(x$loadings),
       xlab = "Tenor", ylab = "Loading",
       main = "PCA Factor Loadings", xaxt = "n", ...)
  axis(1, at = tenor_idx, labels = x$tenors)
  abline(h = 0, lty = 2, col = "grey60")
  grid()

  for (i in seq_len(min(n_comp, 5))[-1]) {
    lines(tenor_idx, x$loadings[, i], type = "b", pch = 19, col = cols[i])
  }

  legend_labels <- vapply(seq_len(min(n_comp, 5)), function(i) {
    lab <- if (i <= 3) labels[i] else paste0("PC", i)
    paste0(lab, " (", round(x$variance_explained[i] * 100, 1), "%)")
  }, character(1))
  legend("topright", legend = legend_labels,
         col = cols[seq_len(min(n_comp, 5))],
         lty = 1, pch = 19, bty = "n", cex = 0.8)

  invisible(x)
}
