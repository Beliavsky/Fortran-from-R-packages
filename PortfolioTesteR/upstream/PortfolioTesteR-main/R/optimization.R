# R/optimization.R - Parameter Grid Optimization with Phase 1 & 2 optimizations


#' Calculate Sharpe Ratio with Frequency Detection
#'
#' @param bt Backtest result object with $returns and (optionally) $dates
#' @return Annualized Sharpe ratio
#' @export
metric_sharpe <- function(bt) {
  returns <- bt$returns
  if (length(returns) == 0 || all(is.na(returns))) return(NA_real_)

  freq_factor <- if (!is.null(bt$dates)) {
    .pt_detect_frequency(as.Date(bt$dates))
  } else {
    52
  }

  mu <- mean(returns, na.rm = TRUE)
  sdv <- stats::sd(returns, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(NA_real_)

  (mu / sdv) * sqrt(freq_factor)
}
#' Run Parameter Grid Optimization (safe + ergonomic)
#'
#' @param prices Data frame with Date + symbol columns
#' @param grid   Data frame (each row = a combo) OR a **named list** of vectors
#' @param builder Function(prices, params, ...) -> weights (Date + symbols)
#' @param metric  Scoring function(backtest) -> numeric. Defaults to metric_sharpe.
#' @param name_prefix String prefix for backtest names
#' @param verbose Logical
#' @param light_mode Logical: speed-ups in backtest
#' @param precompute_returns Logical: precompute log-returns once (light_mode only)
#' @param builder_args List of extra args forwarded to builder (e.g., caches)
#' @param n_cores Integer (kept for API compatibility; ignored here)
#' @param fixed Optional named list of **constant** parameters merged into every combo.
#'   If a name appears in both `grid` and `fixed`, the **fixed value wins** and that
#'   column is **pruned from the grid** (fewer duplicate combos; clearer counts).
#' @return param_grid_result
#' @export
run_param_grid <- function(prices,
                           grid,
                           builder,
                           metric = NULL,
                           name_prefix = "Strategy",
                           verbose = FALSE,
                           light_mode = TRUE,
                           precompute_returns = TRUE,
                           builder_args = list(),
                           n_cores = 1,
                           fixed = NULL) {

  # --- validate inputs -------------------------------------------------------
  if (is.null(metric)) metric <- metric_sharpe
  if (!is.data.frame(prices)) stop("'prices' must be a data frame", call. = FALSE)
  if (!is.function(builder))  stop("'builder' must be a function", call. = FALSE)
  if (!is.function(metric))   stop("'metric' must be a function", call. = FALSE)

  # normalize fixed
  if (!is.null(fixed)) {
    if (!is.list(fixed) || is.null(names(fixed)) || any(names(fixed) == "")) {
      stop("'fixed' must be a named list, e.g., list(top_n = 10)", call. = FALSE)
    }
  } else {
    fixed <- list()
  }

  # Guardrail: forbid cadence/timeframe knobs in the grid
  .pt_check_grid_params(grid, "parameter grid optimization")

  # --- prune collisions BEFORE expand.grid / counting combos -----------------
  # If `grid` is a named list: remove names that appear in `fixed`.
  # If `grid` is a data.frame: drop columns that appear in `fixed`.
  prune_collisions <- function(grid, fixed) {
    if (length(fixed) == 0) return(grid)
    if (is.list(grid) && !is.data.frame(grid)) {
      gnames <- names(grid)
      coll   <- intersect(gnames, names(fixed))
      if (length(coll)) {
        warning("Parameters in both grid and fixed will use the fixed values: ",
                paste(coll, collapse = ", "), call. = FALSE)
        grid[coll] <- NULL
      }
      grid
    } else if (is.data.frame(grid)) {
      coll <- intersect(names(grid), names(fixed))
      if (length(coll)) {
        warning("Parameters in both grid and fixed will use the fixed values: ",
                paste(coll, collapse = ", "), call. = FALSE)
        grid <- grid[, setdiff(names(grid), coll), drop = FALSE]
      }
      grid
    } else {
      stop("'grid' must be a data frame or a named list", call. = FALSE)
    }
  }

  grid <- prune_collisions(grid, fixed)

  # Accept grid as named list (normalize to data.frame)
  if (is.list(grid) && !is.data.frame(grid)) {
    if (is.null(names(grid)) || any(names(grid) == "")) {
      stop("'grid' list must be named, e.g., list(lookback = 1:4, n_top = 1:3)", call. = FALSE)
    }
    grid <- do.call(expand.grid, c(grid, stringsAsFactors = FALSE))
  }
  if (!is.data.frame(grid)) stop("'grid' must be a data frame or a named list", call. = FALSE)
  if (nrow(grid) == 0)      stop("'grid' must have at least one row after pruning", call. = FALSE)

  # Detect whether run_backtest supports 'track_transactions'
  rb_formals  <- tryCatch(formals(run_backtest), error = function(e) NULL)
  rb_has_track<- is.list(rb_formals) && ("track_transactions" %in% names(rb_formals))

  # Phase 1: optional precompute
  if (precompute_returns && light_mode) {
    P <- as.data.frame(prices)
    price_cols <- setdiff(names(P), "Date")
    if (length(price_cols) > 0L && nrow(P) > 1L) {
      price_matrix <- as.matrix(P[, price_cols, drop = FALSE])
      returns_matrix <- diff(log(price_matrix))
      builder_args$.precomputed_returns <- returns_matrix
    }
  }

  results <- vector("list", nrow(grid))
  scores  <- rep(NA_real_, nrow(grid))

  for (i in seq_len(nrow(grid))) {
    # merge fixed -> wins over row values (collisions already pruned above)
    params_row <- as.list(grid[i, , drop = FALSE])
    params     <- utils::modifyList(params_row, fixed)

    if (verbose) {
      cat(sprintf("Testing %d/%d: %s\n",
                  i, nrow(grid),
                  paste(sprintf("%s=%s", names(params), params), collapse = ", ")))
    }

    tryCatch({
      # 1) Build weights
      weights <- if (length(builder_args) > 0) {
        do.call(builder, c(list(prices, params), builder_args))
      } else {
        builder(prices, params)
      }

      # 2) OPTIONAL: validate/align weights (soft)
      if (exists(".pt_validate_weights", mode = "function", inherits = TRUE)) {
        weights2 <- tryCatch(
          .pt_validate_weights(prices, weights, context = sprintf("grid[%d]", i), allow_sparse = TRUE),
          error = function(e) {
            if (verbose) cat("  [validator skipped:", e$message, "]\n")
            NULL
          }
        )
        if (!is.null(weights2)) weights <- weights2
      }

      # 3) Basic sanity
      if (is.null(weights) || nrow(weights) == 0) {
        if (verbose) cat("  [weights empty -> NA score]\n")
        scores[i]   <- NA_real_
        results[[i]]<- list(params = params, score = NA_real_, backtest = NULL)
      } else {
        # 4) Backtest
        rb_args <- list(
          prices = prices,
          weights = weights,
          initial_capital = 100000,
          name = paste0(name_prefix, "_", i),
          verbose = FALSE
        )
        if (light_mode && rb_has_track) rb_args$track_transactions <- FALSE
        bt <- do.call(run_backtest, rb_args)

        # 5) Score
        score <- metric(bt)
        scores[i]    <- score
        results[[i]] <- list(params = params, score = score, backtest = bt)
        if (verbose && is.na(score)) cat("  [metric returned NA]\n")
      }
    }, error = function(e) {
      if (verbose) cat(sprintf("  Error: %s\n", e$message))
      scores[i]    <- NA_real_
      results[[i]] <- list(params = params, score = NA_real_, error = e$message)
    })
  }

  # NA-safe best
  valid_scores <- scores[!is.na(scores)]
  if (length(valid_scores) == 0) {
    best_idx <- NULL; best_params <- NULL; best_score <- NA_real_; best_backtest <- NULL
  } else {
    best_idx      <- which.max(replace(scores, is.na(scores), -Inf))[1]
    best_params   <- results[[best_idx]]$params
    best_score    <- scores[best_idx]
    best_backtest <- results[[best_idx]]$backtest
  }

  structure(
    list(
      best_params = best_params,
      best_idx = best_idx,
      optimization_score = best_score,
      best_backtest = best_backtest,
      all_results = results,
      all_scores  = scores,
      grid        = grid,       # pruned grid actually evaluated
      fixed       = fixed,      # carried for printing
      n_valid     = sum(!is.na(scores)),
      n_tested    = nrow(grid)  # pruned-count -> override case will drop duplicates
    ),
    class = "param_grid_result"
  )
}



#' Print a param_grid_result
#'
#' @param x A `param_grid_result` object returned by [run_param_grid()].
#' @param ... Additional arguments passed to methods (ignored).
#' @return Invisibly returns `x`.
#' @method print param_grid_result
#' @export
print.param_grid_result <- function(x, ...) {
  cat("Parameter Grid Optimization Result\n")
  cat("==================================\n")
  cat("Tested:", x$n_tested, "combinations\n")
  cat("Valid:", x$n_valid,  "results\n")

  if (!is.null(x$fixed) && length(x$fixed)) {
    cat("\nFixed parameters:\n")
    for (nm in names(x$fixed)) {
      cat(sprintf("  %s: %s\n", nm, paste(x$fixed[[nm]], collapse = ",")))
    }
  }

  if (!is.null(x$best_params)) {
    cat("\nBest parameters:\n")
    for (nm in names(x$best_params)) {   # <- removed extra ')'
      cat(sprintf("  %s: %s\n", nm, paste(x$best_params[[nm]], collapse = ",")))
    }
    cat(sprintf("\nBest score: %.4f\n", x$optimization_score))
  } else {
    cat("\nNo valid results found\n")
  }

  invisible(x)
}




# =========================
# Plotting for param_grid_result (generic, robust)
# =========================

# ---- small utilities ----
#' Internal helper: restore par after plotting
#' @keywords internal
#' @noRd
.pt_with_par <- function(expr) {
  old_par <- par(no.readonly = TRUE)
  on.exit(try(par(old_par), silent = TRUE), add = TRUE)
  force(expr)
}

#' Internal helper: default HCL palette
#' @keywords internal
#' @noRd
.hcl_pal <- function(n = 100) grDevices::hcl.colors(n, "YlOrRd", rev = TRUE)

#' Internal helper: match facet values by type
#' @keywords internal
#' @noRd
.match_facet_values <- function(vec, requested) {
  if (is.factor(vec)) {
    chr_req <- as.character(requested)
    chr_req[chr_req %in% levels(vec)]
  } else if (is.numeric(vec)) {
    num_req <- suppressWarnings(as.numeric(requested))
    num_req[!is.na(num_req) & num_req %in% unique(vec)]
  } else {
    chr_req <- as.character(requested)
    chr_req[chr_req %in% as.character(unique(vec))]
  }
}

#' Internal helper: union z-limits across multiple heatmaps
#' @keywords internal
#' @noRd
.zlim_union <- function(hm_list, clip = c(0.02, 0.98)) {
  z_flat <- unlist(lapply(hm_list, function(h) as.numeric(h$z)), use.names = FALSE)
  z_flat <- z_flat[is.finite(z_flat)]
  if (!length(z_flat)) return(c(0,1))
  qs <- stats::quantile(z_flat, probs = clip, na.rm = TRUE)
  range(qs, na.rm = TRUE)
}

#' Internal helper: compute panel grid (rows x cols)
#' @keywords internal
#' @noRd
.panel_dims <- function(k, ncol = NULL) {
  if (!is.null(ncol) && is.finite(ncol) && ncol > 0) {
    nc <- as.integer(ncol); nr <- ceiling(k / nc)
  } else {
    nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  }
  list(nr = nr, nc = nc)
}

# Flatten param_grid_result -> data.frame (params + score + result_idx)
.pg_flatten <- function(opt) {
  rows <- lapply(seq_along(opt$all_results), function(i) {
    r <- opt$all_results[[i]]
    if (!is.null(r$params) && !is.null(r$score)) {
      d <- as.data.frame(r$params, stringsAsFactors = FALSE)
      d$score <- r$score; d$result_idx <- i; d
    } else NULL
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) stop("No valid results to plot.")
  do.call(rbind, rows)
}

# Heatmap builder that accepts numeric/character/factor axes
.make_heatmap <- function(df, xcol, ycol, zcol,
                          agg = c("mean","median","max"),
                          clip = c(0.02, 0.98),
                          na.rm = TRUE) {
  agg <- match.arg(agg)
  f   <- switch(agg, mean = mean, median = median, max = max)

  xi <- if (is.list(df[[xcol]])) unlist(df[[xcol]]) else df[[xcol]]
  yi <- if (is.list(df[[ycol]])) unlist(df[[ycol]]) else df[[ycol]]
  zi <- if (is.list(df[[zcol]])) unlist(df[[zcol]]) else df[[zcol]]

  # aggregate onto y-x grid (rows = y, cols = x)
  z <- tapply(zi, list(yi, xi), f, na.rm = na.rm)
  z <- as.matrix(z)
  if (is.null(dim(z)) || any(dim(z) == 0)) {
    stop("No grid could be formed for ", xcol, " x ", ycol,
         ". Try `fixed=` to reduce parameters or check coverage.")
  }

  cx <- colnames(z); rx <- rownames(z)
  x_num <- suppressWarnings(as.numeric(cx))
  y_num <- suppressWarnings(as.numeric(rx))
  x_is_num <- length(cx) > 0 && all(is.finite(x_num))
  y_is_num <- length(rx) > 0 && all(is.finite(y_num))

  if (x_is_num) { xo <- order(x_num); x_vals <- x_num[xo]; z <- z[, xo, drop=FALSE]; x_labels <- NULL }
  else          { x_vals <- seq_len(ncol(z)); x_labels <- if (length(cx)) cx else as.character(seq_len(ncol(z))) }
  if (y_is_num) { yo <- order(y_num); y_vals <- y_num[yo]; z <- z[yo, , drop=FALSE]; y_labels <- NULL }
  else          { y_vals <- seq_len(nrow(z)); y_labels <- if (length(rx)) rx else as.character(seq_len(nrow(z))) }

  # winsorized zlim for stable colouring
  z_flat <- as.numeric(z); z_flat <- z_flat[is.finite(z_flat)]
  zlim <- if (length(z_flat)) {
    qs <- stats::quantile(z_flat, probs = clip, na.rm = TRUE); range(qs, na.rm = TRUE)
  } else c(0,1)

  list(
    x_vals = as.numeric(x_vals),
    y_vals = as.numeric(y_vals),
    x_labels = x_labels,
    y_labels = y_labels,
    z = z, zlim = zlim
  )
}

# Draw using index grid and transpose so x-axis = columns, y-axis = rows
# Draw using index grid and transpose so x-axis = columns, y-axis = rows
.draw_heatmap_grid <- function(hm,
                               palette = NULL, zlim = NULL,
                               xlab = "", ylab = "", main = NULL,
                               contour = TRUE, label_cells = FALSE, mark_best = TRUE) {
  pal <- if (is.null(palette)) .hcl_pal(100) else palette
  zl  <- if (is.null(zlim)) hm$zlim else zlim

  nx <- ncol(hm$z); ny <- nrow(hm$z)
  xidx <- seq_len(nx); yidx <- seq_len(ny)

  # Suppress base axes (xaxt/yaxt) so we only see our custom tick labels
  graphics::image(
    xidx, yidx, t(hm$z), useRaster = TRUE,
    col = pal, zlim = zl, xlab = xlab, ylab = ylab, main = main,
    axes = FALSE   # <- robustly suppress both axes
  )


  # Custom axis labels (numeric values or categorical labels)
  xlabs <- if (is.null(hm$x_labels)) hm$x_vals else hm$x_labels
  ylabs <- if (is.null(hm$y_labels)) hm$y_vals else hm$y_labels
  graphics::axis(1, at = xidx, labels = xlabs, las = 2)
  graphics::axis(2, at = yidx, labels = ylabs, las = 2)

  # Contours: match the transpose used above
  if (contour && ny >= 5 && nx >= 5) {
    graphics::contour(xidx, yidx, t(hm$z), add = TRUE, drawlabels = FALSE)
  }

  # Cell value labels: label the transposed matrix to match the plot orientation
  if (label_cells) {
    tx <- rep(xidx, each = ny)
    ty <- rep(yidx, times = nx)
    graphics::text(tx, ty, labels = sprintf("%.3f", as.vector(t(hm$z))), cex = 0.8)
  }

  # Mark best: swap i/j to account for transpose (you already did this correctly)
  if (mark_best && any(is.finite(hm$z))) {
    ij <- which(hm$z == max(hm$z, na.rm = TRUE), arr.ind = TRUE)[1, ]
    graphics::points(ij[2], ij[1], pch = 21, bg = "white", cex = 1.2)
  }

  graphics::box()
  invisible(zl)
}

#' Internal helper: draw legend colorbar on the right side
#' @keywords internal
#' @noRd
#' @importFrom graphics axis box image mtext
.add_colorbar <- function(zlim, palette, n = 100, main = NULL) {
  pal <- if (is.null(palette)) .hcl_pal(n) else palette
  zseq <- seq(zlim[1], zlim[2], length.out = n)
  par(mar = c(5, 1.5, if (is.null(main)) 4 else 5, 3))
  image(x = 1, y = zseq, z = matrix(zseq, nrow = 1),
        col = pal, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(4, las = 1); box()
  if (!is.null(main)) mtext(main, side = 3, line = 1)
}

# =========================
# plot() method
# =========================
# --- helpers specifically for surfaces (add ABOVE the plot method) ---
#' Internal helper: sanitize z for 3D surface
#' @keywords internal
#' @noRd
.pt_sanitize_z_for_surface <- function(z) {
  ok <- is.finite(z)
  if (!any(ok)) return(list(z = z, all_bad = TRUE))
  z_fill <- z
  z_fill[!ok] <- stats::median(z[ok], na.rm = TRUE)
  list(z = z_fill, all_bad = FALSE)
}

#' Internal helper: surface colors for persp
#' @keywords internal
#' @noRd
.pt_surface_colours <- function(z, zlim, palette) {
  pal <- if (is.null(palette)) grDevices::hcl.colors(100, "YlOrRd", rev = TRUE) else palette
  zf <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)]) / 4
  zf_clamped <- pmin(pmax(zf, zlim[1]), zlim[2])
  brks <- seq(zlim[1], zlim[2], length.out = length(pal) + 1)
  idx  <- cut(zf_clamped, breaks = brks, include.lowest = TRUE)
  pal[as.integer(idx)]
}

#' Plot Parameter Grid Results (1D/2D/3D and Facets)
#'
#' Generic plotter for objects returned by [run_param_grid()]. Supported types:
#' \itemize{
#'   \item "line": 1D metric vs one parameter.
#'   \item "heatmap": 2D heatmap over two parameters.
#'   \item "surface": 3D surface (persp) over two parameters.
#'   \item "slices": 2D heatmaps faceted by a third parameter.
#'   \item "surface_slices": 3D surfaces faceted by a third parameter.
#'   \item "auto": chosen from the above based on the number of parameters.
#' }
#'
#' @param x A \code{param_grid_result}.
#' @param y Ignored.
#' @param type One of "auto","line","heatmap","surface","slices","surface_slices".
#' @param metric Column name to plot (defaults to "score" if present).
#' @param params Character vector of parameter columns to use for axes/facets.
#'   If \code{NULL}, uses all parameter columns detected in \code{x$all_results}.
#' @param fixed Optional named list of parameter values to condition on.
#'   Rows not matching are dropped before plotting.
#' @param agg Aggregation when multiple rows map to a cell:
#'   "mean","median", or "max".
#' @param na.rm Logical; drop NA metric rows before plotting.
#' @param palette Optional color vector used for heatmaps/surfaces.
#'   Defaults to \code{grDevices::hcl.colors(..., "YlOrRd", rev = TRUE)}.
#' @param zlim Optional two-element numeric range for color scaling.
#' @param clip Two quantiles used to winsorize z-limits for stable coloring.
#' @param main,sub,xlab,ylab Base plotting annotations.
#' @param ... Additional options depending on \code{type}, e.g.:
#'   \code{legend}, \code{label_cells}, \code{contour}, \code{mark_best},
#'   \code{theta}, \code{phi}, \code{shade}, \code{expand}, \code{ticktype},
#'   \code{shared_zlim}, \code{facet}, \code{facet_values}, \code{ncol},
#'   \code{debug}, \code{impute_for_surface}.
#'
#' @return An invisible list describing the plot (kind/params/zlim/etc.).
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   b <- function(prices, params, ...) {
#'     weight_equally(filter_top_n(calc_momentum(prices, params$lookback), 10))
#'   }
#'   opt <- run_param_grid(
#'     prices  = sample_prices_weekly,
#'     grid    = list(lookback = c(8, 12, 26)),
#'     builder = b
#'   )
#'   plot(opt, type = "line", params = "lookback")
#' }
#'
#' @seealso [run_param_grid()], [print.param_grid_result()]
#'
#' @method plot param_grid_result
#' @export
#' @importFrom grDevices hcl.colors
#' @importFrom graphics axis box image layout mtext persp plot.new points title
#' @importFrom stats median quantile
#' @importFrom utils head
plot.param_grid_result <- function(
    x, y = NULL,
    type = c("auto","line","heatmap","surface","slices","surface_slices"),
    metric = NULL,
    params = NULL,
    fixed  = list(),
    agg    = c("mean","median","max"),
    na.rm  = TRUE,
    palette = NULL,
    zlim    = NULL,
    clip    = c(0.02, 0.98),
    main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ...
) {
  type <- match.arg(type); agg <- match.arg(agg)
  dots <- list(...)
  dbg  <- isTRUE(dots$debug)

  # ---- flatten grid results ----
  results_list <- lapply(seq_along(x$all_results), function(i) {
    r <- x$all_results[[i]]
    if (!is.null(r$params) && !is.null(r$score)) {
      d <- as.data.frame(r$params, stringsAsFactors = FALSE)
      d$score <- r$score; d$result_idx <- i; d
    } else NULL
  })
  results_list <- Filter(Negate(is.null), results_list)
  if (!length(results_list)) stop("No valid results to plot.")
  df <- do.call(rbind, results_list)

  if (dbg) {
    message("DEBUG: head(flat)"); print(utils::head(df, 10))
    message("DEBUG: unique params"); print(lapply(setdiff(names(df), c("score","result_idx")), unique))
  }

  # ---- metric column ----
  if (is.null(metric)) metric <- if ("score" %in% names(df)) "score" else names(df)[ncol(df)]
  if (!metric %in% names(df)) stop("Metric not found: ", metric)
  if (na.rm) df <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (!nrow(df)) stop("No non-NA metric rows to plot.")

  # ---- candidate params & fixed filters ----
  param_cols_all <- setdiff(names(df), c(metric, "result_idx"))
  if (is.null(params)) params <- param_cols_all
  if (length(fixed)) {
    for (nm in names(fixed)) if (nm %in% names(df)) df <- df[df[[nm]] == fixed[[nm]], , drop = FALSE]
    if (!nrow(df)) stop("No rows match 'fixed' constraints.")
  }

  # ---- auto type ----
  if (type == "auto") {
    npar <- length(params)
    type <- if (npar <= 1) "line" else if (npar == 2) "heatmap" else "slices"
  }

  # small helper to always reset layout
  layout_reset <- function() on.exit(try(layout(1), silent = TRUE), add = TRUE)

  # ---------- LINE (1 param) ----------
  if (type == "line") {
    if (!length(params)) stop("Provide one parameter for type='line'.")
    xcol <- params[1L]
    if (!xcol %in% names(df)) stop("Parameter not found: ", xcol)
    f <- switch(agg, mean=mean, median=median, max=max)
    xi <- if (is.list(df[[xcol]])) unlist(df[[xcol]]) else df[[xcol]]
    yagg <- tapply(df[[metric]], xi, f, na.rm = TRUE)
    xv <- suppressWarnings(as.numeric(names(yagg)))
    x_is_num <- all(is.finite(xv))
    if (!x_is_num) xv <- seq_along(yagg)
    ord <- order(xv); xv <- xv[ord]; yv <- as.numeric(yagg)[ord]

    .pt_with_par({
      graphics::plot(xv, yv, type = "b", pch = 19,
                     xlab = xcol, ylab = metric,
                     main = if (is.null(main)) sprintf("%s vs %s (%s)", metric, xcol, agg) else main)
      if (!x_is_num) graphics::axis(1, at = xv, labels = names(yagg)[ord], las = 2)
      b <- which.max(yv); graphics::points(xv[b], yv[b], pch = 21, bg = "white", cex = 1.3)
      graphics::box()
    })
    return(invisible(list(kind="line", param=xcol, agg=agg)))
  }

  # ---------- HEATMAP (2 params) ----------
  if (type == "heatmap") {
    if (length(params) < 2) stop("Need two parameters for heatmap.")
    xcol <- params[1L]; ycol <- params[2L]
    if (!all(c(xcol, ycol) %in% names(df))) stop("Params not found in results.")

    hm  <- .make_heatmap(df, xcol, ycol, metric, agg = agg, clip = clip, na.rm = na.rm)
    pal <- if (is.null(palette)) .hcl_pal(100) else palette
    if (!is.null(zlim)) hm$zlim <- zlim

    nr <- nrow(hm$z); nc <- ncol(hm$z)
    draw_contour <- if (!is.null(dots$contour)) isTRUE(dots$contour) else (nr >= 5 && nc >= 5)
    add_legend   <- isTRUE(dots$legend)
    mark_best    <- if (is.null(dots$mark_best)) TRUE else isTRUE(dots$mark_best)
    label_cells  <- isTRUE(dots$label_cells)

    if (dbg) {
      message("DEBUG[heatmap]: dim(z) = ", paste(dim(hm$z), collapse="x"),
              " | x_vals=", paste(hm$x_vals, collapse=","),
              " | y_vals=", paste(hm$y_vals, collapse=","))
    }

    .pt_with_par({
      layout_reset()
      if (add_legend) {
        layout(matrix(c(1,2), 1, 2), widths = c(4,1))
        par(mar = c(5,5,4,2) + 0.1)
      }
      .draw_heatmap_grid(
        hm, palette = pal, zlim = hm$zlim,
        xlab = if (is.null(xlab)) xcol else xlab,
        ylab = if (is.null(ylab)) ycol else ylab,
        main = if (is.null(main)) sprintf("%s heatmap (%s)", metric, agg) else main,
        contour = draw_contour, label_cells = label_cells, mark_best = mark_best
      )
      if (add_legend) .add_colorbar(hm$zlim, palette = pal, main = "Score")
    })
    return(invisible(list(kind="heatmap", params=c(xcol,ycol), agg=agg,
                          grid=list(x=hm$x_vals, y=hm$y_vals, z=hm$z, zlim=hm$zlim))))
  }

  # ---------- SURFACE (2 params, 3D) ----------
  if (type == "surface") {
    if (length(params) < 2) stop("Need two parameters for surface.")
    xcol <- params[1L]; ycol <- params[2L]
    if (!all(c(xcol, ycol) %in% names(df))) stop("Params not found in results.")

    hm <- .make_heatmap(df, xcol, ycol, metric, agg = agg, clip = clip, na.rm = na.rm)
    if (nrow(hm$z) < 2 || ncol(hm$z) < 2) stop("Surface needs at least a 2x2 grid.")

    impute <- if (is.null(dots$impute_for_surface)) TRUE else isTRUE(dots$impute_for_surface)
    sani <- .pt_sanitize_z_for_surface(hm$z)
    if (!impute && (sani$all_bad || any(!is.finite(hm$z)))) {
      stop("Surface has non-finite cells. Set impute_for_surface=TRUE or adjust grid/data.")
    }
    z0 <- if (impute) sani$z else hm$z

    # persp() expects rows = length(x), cols = length(y)  -> transpose
    zt <- t(z0)

    zl_use  <- if (!is.null(zlim)) zlim else range(zt[is.finite(zt)], na.rm = TRUE)
    col_vec <- .pt_surface_colours(zt, zl_use, palette)

    if (dbg) {
      message("DEBUG[surface]: dim(z)=", paste(dim(hm$z), collapse="x"),
              " | dim(tz)=", paste(dim(zt), collapse="x"),
              " | len(x)=", length(hm$x_vals), " | len(y)=", length(hm$y_vals),
              " | zlim=", paste(round(zl_use,4), collapse=" .. "))
    }

    theta <- if (is.null(dots$theta)) 40 else dots$theta
    phi   <- if (is.null(dots$phi))   25 else dots$phi
    shade <- if (is.null(dots$shade)) 0.5 else dots$shade
    expand<- if (is.null(dots$expand))0.6 else dots$expand
    ticktype<- if (is.null(dots$ticktype)) "detailed" else dots$ticktype

    .pt_with_par({
      par(mar = c(4,4,4,2) + 0.1)
      res <- graphics::persp(x = hm$x_vals, y = hm$y_vals, z = zt,
                             theta = theta, phi = phi, expand = expand,
                             shade = shade, zlim = zl_use,
                             col = as.vector(col_vec), border = NA,
                             xlab = if (is.null(xlab)) xcol else xlab,
                             ylab = if (is.null(ylab)) ycol else ylab,
                             zlab = metric, ticktype = ticktype,
                             main = if (is.null(main)) sprintf("%s surface (%s)", metric, agg) else main)
      if (any(is.finite(zt))) {
        ij <- which(zt == max(zt, na.rm = TRUE), arr.ind = TRUE)[1, ]
        p  <- grDevices::trans3d(hm$x_vals[ij[1]], hm$y_vals[ij[2]], zt[ij[1], ij[2]], res)
        graphics::points(p, pch = 21, bg = "white", cex = 1.2)
      }
    })
    return(invisible(list(kind="surface", params=c(xcol,ycol), agg=agg, zlim=zl_use)))
  }

  # ---------- SLICES (3 params: X/Y + facet; >3 -> require fixed=) ----------
  if (type == "slices") {
    if (!length(params)) stop("Provide parameters for slices.")
    # Pick facet name without partial matching; allow default to 3rd param
    facet <- dots[["facet"]]
    if (length(params) >= 3L && is.null(facet)) facet <- params[3L]
    if (length(params) == 2L && is.null(facet)) {
      stop("Slices require 3 parameters (X, Y, facet).")
    }
    xcol <- params[1L]
    ycol <- params[2L]
    fcol <- as.character(facet)[1L]

    needed  <- c(xcol, ycol, fcol)
    missing <- setdiff(needed, names(df))
    if (length(missing)) {
      stop("Params not found: ", paste(missing, collapse = ", "),
           "\nAvailable: ", paste(names(df), collapse = ", "))
    }

    # Accept facet_values like c(5,20) and coerce to the facet column's type
    fvec <- df[[fcol]]
    fv <- if (!is.null(dots[["facet_values"]])) {
      .match_facet_values(fvec, dots[["facet_values"]])
    } else {
      unique(fvec)
    }
    if (!length(fv)) stop(sprintf("No facet levels to plot for %s after filtering.", fcol))


    hm_list <- lapply(fv, function(v) {
      dfi <- df[df[[fcol]] == v, , drop = FALSE]
      if (!nrow(dfi)) return(NULL)
      .make_heatmap(dfi, xcol, ycol, metric, agg = agg, clip = clip, na.rm = na.rm)
    })
    names(hm_list) <- as.character(fv)
    hm_list <- Filter(Negate(is.null), hm_list)
    if (!length(hm_list)) stop("No data to plot for any facet level.")

    pal        <- if (is.null(palette)) .hcl_pal(100) else palette
    shared     <- if (!is.null(dots$shared_zlim)) isTRUE(dots$shared_zlim) else TRUE
    add_legend <- isTRUE(dots$legend)
    draw_contour <- if (!is.null(dots$contour)) isTRUE(dots$contour) else {
      nr <- nrow(hm_list[[1]]$z); nc <- ncol(hm_list[[1]]$z); (nr >= 5 && nc >= 5)
    }
    mark_best   <- if (is.null(dots$mark_best)) TRUE else isTRUE(dots$mark_best)
    label_cells <- isTRUE(dots$label_cells)
    if (shared) zl <- if (is.null(zlim)) .zlim_union(hm_list, clip = clip) else zlim

    if (dbg) {
      message("DEBUG[slices]: facet levels -> ", paste(names(hm_list), collapse=", "))
    }

    dims <- .panel_dims(length(hm_list), ncol = if (!is.null(dots$ncol)) dots$ncol else NULL)
    nr <- dims$nr; nc <- dims$nc

    .pt_with_par({
      layout_reset()
      if (add_legend) {
        base <- matrix(seq_len(nr * nc), nrow = nr, ncol = nc, byrow = TRUE)
        mat  <- cbind(base, rep(nr * nc + 1, nr))
        layout(mat, widths = c(rep(4, nc), 1))
      } else {
        layout(matrix(seq_len(nr * nc), nrow = nr, ncol = nc, byrow = TRUE))
      }
      par(mar = c(4, 4, 3, 1) + 0.1)

      for (nm in names(hm_list)) {
        hm <- hm_list[[nm]]
        zl_use <- if (!is.null(zlim)) zlim else if (shared) zl else hm$zlim
        .draw_heatmap_grid(
          hm, palette = pal, zlim = zl_use,
          xlab = if (is.null(xlab)) xcol else xlab,
          ylab = if (is.null(ylab)) ycol else ylab,
          main = sprintf("%s = %s", fcol, nm),
          contour = draw_contour, label_cells = label_cells, mark_best = mark_best
        )
      }

      if (add_legend) .add_colorbar(if (!is.null(zlim)) zlim else if (shared) zl else hm_list[[1]]$zlim,
                                    palette = pal, main = if (is.null(main)) "Score" else main)
    })
    return(invisible(list(kind = "slices",
                          axes   = c(xcol, ycol),
                          facet  = fcol,
                          levels = names(hm_list))))
  }

  # ---------- SURFACE_SLICES (3 params: X/Y + facet, 3D panels) ----------
  if (type == "surface_slices") {
    if (!length(params)) stop("Provide parameters for surface_slices.")
    facet <- dots[["facet"]]
    if (is.null(facet)) {
      if (length(params) < 3L) stop("surface_slices needs 3 parameters (X, Y, facet).")
      facet <- params[3L]
    }
    xcol <- params[1L]
    ycol <- params[2L]
    fcol <- as.character(facet)[1L]

    needed  <- c(xcol, ycol, fcol)
    missing <- setdiff(needed, names(df))
    if (length(missing)) {
      stop("Params not found: ", paste(missing, collapse = ", "),
           "\nAvailable: ", paste(names(df), collapse = ", "))
    }

    fvec <- df[[fcol]]
    fv <- if (!is.null(dots[["facet_values"]])) {
      .match_facet_values(fvec, dots[["facet_values"]])
    } else {
      unique(fvec)
    }
    if (!length(fv)) stop(sprintf("No facet levels to plot for %s.", fcol))


    hm_list <- lapply(fv, function(v) {
      dfi <- df[df[[fcol]] == v, , drop = FALSE]
      if (!nrow(dfi)) return(NULL)
      .make_heatmap(dfi, xcol, ycol, metric, agg = agg, clip = clip, na.rm = na.rm)
    })
    names(hm_list) <- as.character(fv)
    hm_list <- Filter(Negate(is.null), hm_list)
    if (!length(hm_list)) stop("No data to plot for any facet level.")

    shared  <- if (!is.null(dots$shared_zlim)) isTRUE(dots$shared_zlim) else TRUE
    if (shared) {
      z_all <- unlist(lapply(hm_list, function(h) as.numeric(h$z)), use.names = FALSE)
      z_all <- z_all[is.finite(z_all)]
      zl_shared <- if (length(z_all)) range(z_all, na.rm = TRUE) else c(0,1)
    }

    theta <- if (is.null(dots$theta)) 40 else dots$theta
    phi   <- if (is.null(dots$phi))   25 else dots$phi
    shade <- if (is.null(dots$shade)) 0.5 else dots$shade
    expand<- if (is.null(dots$expand))0.6 else dots$expand
    ticktype<- if (is.null(dots$ticktype)) "detailed" else dots$ticktype

    dims <- .panel_dims(length(hm_list), ncol = if (!is.null(dots$ncol)) dots$ncol else NULL)
    nr <- dims$nr; nc <- dims$nc

    if (dbg) {
      message("DEBUG[surface_slices]: facet levels -> ", paste(names(hm_list), collapse=", "))
    }

    .pt_with_par({
      on.exit(try(layout(1), silent = TRUE), add = TRUE)
      layout(matrix(seq_len(nr * nc), nrow = nr, ncol = nc, byrow = TRUE))
      par(mar = c(2.5, 2.5, 3, 0.5) + 0.1)

      for (nm in names(hm_list)) {
        hm <- hm_list[[nm]]

        if (nrow(hm$z) < 2 || ncol(hm$z) < 2) {
          plot.new(); title(main = sprintf("%s = %s (insufficient grid)", fcol, nm))
          next
        }

        impute <- if (is.null(dots$impute_for_surface)) TRUE else isTRUE(dots$impute_for_surface)
        sani <- .pt_sanitize_z_for_surface(hm$z)
        if (!impute && (sani$all_bad || any(!is.finite(hm$z)))) {
          plot.new(); title(main = sprintf("%s = %s (NA/Inf in z)", fcol, nm))
          next
        }
        z0 <- if (impute) sani$z else hm$z

        # persp() expects rows = length(x), cols = length(y)  -> transpose
        zt <- t(z0)

        zl_use <- if (!is.null(zlim)) {
          zlim
        } else if (shared) {
          zl_shared
        } else {
          range(zt[is.finite(zt)], na.rm = TRUE)
        }

        col_vec <- .pt_surface_colours(zt, zl_use, palette)

        if (dbg) {
          message(sprintf("DEBUG[surface_slices][%s]: dim(z)=%sx%s, dim(tz)=%sx%s, x=%d, y=%d, zlim=%0.4f..%0.4f",
                          nm, nrow(hm$z), ncol(hm$z), nrow(zt), ncol(zt),
                          length(hm$x_vals), length(hm$y_vals), zl_use[1], zl_use[2]))
        }

        res <- graphics::persp(x = hm$x_vals, y = hm$y_vals, z = zt,
                               theta = theta, phi = phi, expand = expand,
                               shade = shade, zlim = zl_use,
                               col = as.vector(col_vec), border = NA,
                               xlab = xcol, ylab = ycol, zlab = metric, ticktype = ticktype,
                               main = sprintf("%s = %s", fcol, nm))
        if (any(is.finite(zt))) {
          ij <- which(zt == max(zt, na.rm = TRUE), arr.ind = TRUE)[1, ]
          p  <- grDevices::trans3d(hm$x_vals[ij[1]], hm$y_vals[ij[2]], zt[ij[1], ij[2]], res)
          graphics::points(p, pch = 21, bg = "white", cex = 1.1)
        }
      }
    })

    return(invisible(list(kind="surface_slices",
                          axes=c(xcol,ycol), facet=fcol, levels = names(hm_list))))
  }

  message(sprintf("Plot type '%s' not implemented here.", type))
  invisible(list(type=type, params=params))
}
