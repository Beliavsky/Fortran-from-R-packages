#' Pairwise Bivariate Kullback-Leibler Divergence
#'
#' Computes all directed kernel KL divergences among two or more groups on one
#' common integration grid. Each group density and bandwidth is evaluated only
#' once.
#'
#' @param x Either a list of bivariate samples or a single numeric matrix or
#'   data frame with exactly two columns.
#' @param group A grouping vector when `x` is a single sample. Must be `NULL`
#'   when `x` is a list.
#' @param H Optional list of bandwidth matrices, one per group.
#' @inheritParams biv_kld
#'
#' @return A numeric matrix. Row `i`, column `j` is the estimated directed
#'   divergence from group `i` to group `j`. Bandwidths, grid information, and
#'   standardization are stored as attributes.
#'
#' @examples
#' measurements <- iris[, c("Sepal.Length", "Sepal.Width")]
#' kld <- biv_kld_matrix(measurements, iris$Species,
#'                       bandwidth = "normal", grid_size = 40)
#' round(kld, 3)
#'
#' @export
biv_kld_matrix <- function(x, group = NULL, H = NULL,
                           bandwidth = c("scv", "normal"),
                           grid_size = 100L, range = NULL,
                           standardize = c("none", "pooled", "separate")) {
    bandwidth <- match.arg(bandwidth)
    standardize <- match.arg(standardize)
    grid_size <- .validate_grid_size(grid_size)

    if (is.list(x) && !is.data.frame(x)) {
        if (!is.null(group)) {
            stop("'group' must be NULL when 'x' is a list.", call. = FALSE)
        }
        samples <- x
        group_names <- names(samples)
        if (is.null(group_names) || any(!nzchar(group_names))) {
            group_names <- paste0("Group", seq_along(samples))
        }
    } else {
        x <- .as_bivariate_sample(x, "x")
        if (is.null(group) || length(group) != nrow(x) || anyNA(group)) {
            stop("'group' must contain one non-missing value per row of 'x'.",
                 call. = FALSE)
        }
        rows <- split(seq_len(nrow(x)), group, drop = TRUE)
        group_names <- names(rows)
        samples <- lapply(rows, function(i) x[i, , drop = FALSE])
    }
    if (length(samples) < 2L) {
        stop("At least two groups are required.", call. = FALSE)
    }
    if (anyDuplicated(group_names)) {
        stop("Group names must be unique.", call. = FALSE)
    }
    samples <- Map(function(z, nm) .as_bivariate_sample(z, nm),
                   samples, group_names)
    samples <- .standardize_samples(samples, standardize)

    if (is.null(H)) {
        H <- rep(list(NULL), length(samples))
    }
    if (!is.list(H) || length(H) != length(samples)) {
        stop("'H' must be NULL or a list with one bandwidth per group.",
             call. = FALSE)
    }
    bandwidths <- Map(function(z, h, nm) {
        .select_bandwidth(z, h, bandwidth, paste0("H[[", nm, "]]"))
    }, samples, H, group_names)

    grid <- .integration_grid(samples, bandwidths, grid_size, range)
    densities <- Map(function(z, h) .density_on_grid(z, h, grid),
                     samples, bandwidths)
    ng <- length(samples)
    answer <- matrix(0, nrow = ng, ncol = ng,
                     dimnames = list(group_names, group_names))
    for (i in seq_len(ng)) {
        for (j in seq_len(ng)) {
            if (i != j) {
                answer[i, j] <- .integrated_kld(
                    densities[[i]], densities[[j]], grid$cell_area)
            }
        }
    }
    attr(answer, "bandwidth") <- stats::setNames(bandwidths, group_names)
    attr(answer, "grid") <- grid
    attr(answer, "standardize") <- standardize
    class(answer) <- c("bivkld_matrix", "matrix", "array")
    answer
}

#' @export
print.bivkld_matrix <- function(x, ...) {
    shown <- x
    attr(shown, "bandwidth") <- NULL
    attr(shown, "grid") <- NULL
    attr(shown, "standardize") <- NULL
    class(shown) <- c("matrix", "array")
    print(shown, ...)
    invisible(x)
}
