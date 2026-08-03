#' @keywords internal
#' @importFrom stats predict
NULL

# ===========================================
# R/ml_helper.R  --- minimal, general ML helpers
# Safe to source from R/ml.R during development
# ===========================================

# ---- tiny internals ---------------------------------------------------------
# Don't redefine if already present elsewhere in the package
#' Internal infix: A or B if A is NULL
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Internal utilities
#' @keywords internal
#' @noRd
.ml_symcols <- function(D) setdiff(names(D), "Date")
#' @keywords internal
#' @noRd
.ml_first  <- function(x, default = NULL) if (length(x)) x[[1]] else default


#' Panel-safe binary operation on aligned wide panels
#'
#' Applies an elementwise binary operator to two date-aligned wide panels
#' (first column `Date`, other columns are symbols), preserving the `Date`
#' column and a consistent symbol set. Supports intersection or union of
#' column sets; missing entries introduced by `how="union"` are filled.
#'
#' @param A,B Data frames with a `Date` column and one column per symbol.
#' @param op Binary function to apply elementwise (e.g., `*`, `/`, `+`).
#' @param how Character; `"intersect"` (default) or `"union"` for the set
#'   of symbol columns to operate on.
#' @param fill Numeric; value used to fill gaps when `how="union"`.
#'
#' @return A data.frame with `Date` and the operated symbol columns.
#' @examples
#' \dontrun{
#' out <- ml_panel_op(mom12_panel, vol_panel, op = `*`)
#' }
#' @export
ml_panel_op <- function(A, B, op = `*`, how = c("intersect","union"), fill = NA_real_) {
  # Always work with plain data.frames (avoid data.table's j-eval)
  A <- as.data.frame(A, check.names = FALSE)
  B <- as.data.frame(B, check.names = FALSE)

  # Harmonize Date types and check alignment
  if (!inherits(A$Date, "Date")) A$Date <- as.Date(A$Date)
  if (!inherits(B$Date, "Date")) B$Date <- as.Date(B$Date)
  stopifnot(identical(A$Date, B$Date))

  how   <- match.arg(how)
  symsA <- setdiff(names(A), "Date")
  symsB <- setdiff(names(B), "Date")
  syms  <- if (how == "intersect") intersect(symsA, symsB) else union(symsA, symsB)
  if (!length(syms)) stop("ml_panel_op(): no symbol columns after ", how, ".")

  # Build numeric matrices column-by-column (works even if some cols are missing)
  if (how == "intersect") {
    aM <- as.matrix(A[, syms, drop = FALSE])
    bM <- as.matrix(B[, syms, drop = FALSE])
  } else {  # union
    aM <- matrix(fill, nrow = nrow(A), ncol = length(syms), dimnames = list(NULL, syms))
    bM <- matrix(fill, nrow = nrow(B), ncol = length(syms), dimnames = list(NULL, syms))
    existA <- intersect(syms, symsA)
    existB <- intersect(syms, symsB)
    if (length(existA)) aM[, existA] <- as.matrix(A[, existA, drop = FALSE])
    if (length(existB)) bM[, existB] <- as.matrix(B[, existB, drop = FALSE])
  }

  out <- data.frame(Date = A$Date, check.names = FALSE)
  out[syms] <- op(aM, bM)
  out
}

#' Reduce multiple panels with a binary operator
#'
#' Folds a list of panels using `ml_panel_op()` across a set of named panels.
#'
#' @param features List of panels (each a wide data frame with `Date`).
#' @param panels Character vector of names in `features` to reduce (length \eqn{\ge}{>=} 2).
#' @param op Binary function to apply elementwise.
#' @param how Column-set policy passed to [ml_panel_op()].
#' @param fill Fill value for `how="union"`.
#'
#' @return A data.frame panel (wide) with the reduced result.
#' @examples
#' \dontrun{
#' # product of three panels
#' prod_panel <- ml_panel_reduce(X, c("mom12","vol","rsi14"), op = `*`)
#' }
#' @export
ml_panel_reduce <- function(features, panels, op = `*`, how = "intersect", fill = NA_real_) {
  stopifnot(length(panels) >= 2L)
  res <- features[[panels[[1]]]]
  for (k in 2:length(panels)) {
    nxt <- features[[panels[[k]]]]
    stopifnot(!is.null(res), !is.null(nxt))
    res <- ml_panel_op(res, nxt, op = op, how = how, fill = fill)
  }
  res
}


#' Add interaction panels to a feature list
#'
#' Builds new panels from existing ones via elementwise operations.
#' Specification accepts either a shorthand `list(new = c("A","B"))` (defaults to product),
#' or a structured form `list(new = list(panels=c("A","B","C"), op=`/`, how="intersect", fill=NA))`.
#'
#' @param features List of existing panels (wide data frames with `Date`).
#' @param interactions Named list describing interactions to add.
#'
#' @return The input `features` list with additional interaction panels.
#' @examples
#' \dontrun{
#' X2 <- ml_add_interactions(X, list(mom_vol = c("mom12","vol")))
#' }
#' @export
ml_add_interactions <- function(features, interactions) {
  for (nm in names(interactions)) {
    spec <- interactions[[nm]]
    if (is.character(spec) && length(spec) == 2L) {
      features[[nm]] <- ml_panel_reduce(features, spec, op = `*`)
    } else if (is.list(spec)) {
      panels <- spec$panels %||% .ml_first(spec)
      op     <- spec$op     %||% `*`
      how    <- spec$how    %||% "intersect"
      fill   <- if (!is.null(spec$fill)) spec$fill else NA_real_
      stopifnot(is.character(panels), length(panels) >= 2L)
      features[[nm]] <- ml_panel_reduce(features, panels, op = op, how = how, fill = fill)
    } else {
      stop(sprintf("Unsupported interaction spec for '%s'", nm))
    }
  }
  features
}




#' Prepare tabular features (weekly + aligned daily volatility)
#'
#' Constructs a minimal, leakage-safe set of tabular features commonly used
#' in cross-sectional ML: weekly momentum (12/26/52), RSI(14), distance from
#' 20-week MA, and daily rolling volatility aligned to weekly dates
#' (tokens `vol{N}d_walign`, e.g., `"vol20d_walign"`).
#'
#' All outputs are **lagged by one period** to avoid look-ahead in backtests.
#'
#' @param prices_weekly Wide panel with `Date` and symbol columns (weekly).
#' @param prices_daily Optional wide panel (daily) if `vol*d_walign` are included.
#' @param include Character vector of feature tokens to include.
#' @param interactions Optional named list passed to [ml_add_interactions()].
#'
#' @return A named list of panels (each a wide data.frame with `Date`).
#' @examples
#' \dontrun{
#' X <- ml_prepare_features(sample_prices_weekly, sample_prices_daily,
#'                          include = c("mom12","vol20d_walign","rsi14"))
#' }
#' @export
ml_prepare_features <- function(prices_weekly,
                                prices_daily = NULL,
                                include = c("mom12","mom26","mom52",
                                            "vol20d_walign","dist20","rsi14"),
                                interactions = NULL) {
  out <- list()

  # momentum (weekly panels)
  if ("mom12" %in% include) out$mom12 <- panel_lag(calc_momentum(prices_weekly, 12L), 1L)
  if ("mom26" %in% include) out$mom26 <- panel_lag(calc_momentum(prices_weekly, 26L), 1L)
  if ("mom52" %in% include) out$mom52 <- panel_lag(calc_momentum(prices_weekly, 52L), 1L)

  # daily volatility aligned to weekly dates
  vol_keys <- grep("^vol(\\d+)d_walign$", include, value = TRUE)
  if (length(vol_keys)) {
    stopifnot(!is.null(prices_daily))
    weekly_dates <- if (inherits(prices_weekly$Date, "Date")) prices_weekly$Date else as.Date(prices_weekly$Date)
    for (vk in vol_keys) {
      N <- as.integer(sub("^vol(\\d+)d_walign$", "\\1", vk))
      vol_d <- calc_rolling_volatility(prices_daily, N)
      if (!inherits(vol_d$Date, "Date")) vol_d$Date <- as.Date(vol_d$Date)
      aligned <- align_to_timeframe(vol_d, weekly_dates, "forward_fill")
      nm <- if (N == 20L) "vol" else paste0("vol", N, "d")
      out[[nm]] <- panel_lag(aligned, 1L)
    }
  }

  # distance to MA(20); RSI(14)
  if ("dist20" %in% include) {
    ma20 <- calc_moving_average(prices_weekly, 20L)
    out$dist20 <- panel_lag(calc_distance(prices_weekly, ma20), 1L)
  }
  if ("rsi14" %in% include) out$rsi14 <- panel_lag(calc_rsi(prices_weekly, 14L), 1L)

  # interactions
  if (!is.null(interactions) && length(interactions)) {
    out <- ml_add_interactions(out, interactions)
  }
  out
}



#' Model factory for tabular cross-sectional learners
#'
#' Returns a pair of closures `fit(X,y)` / `predict(model, X)` implementing
#' a chosen learner. Implementations are NA-aware and conservative:
#' `glmnet` ridge drops rows with any non-finite input; `ranger` and `xgboost`
#' keep NA in `X` as missing; the linear baseline uses `lm.fit`.
#'
#' @param type One of `"ridge"`, `"rf"`, `"xgboost"`, `"linear"`.
#' @param params List of model parameters (passed to backend; used by xgboost).
#' @param nrounds Integer boosting rounds (xgboost).
#' @param ... Additional arguments forwarded to the backend.
#'
#' @details
#' Optional dependencies: `glmnet` (ridge), `ranger` (rf), `xgboost` (xgboost).
#' If a backend is not available, use `"linear"` or install the package.
#'
#' @return A list with functions `fit` and `predict`.
#' @examples
#' \dontrun{
#' ridge <- ml_make_model("ridge")
#' m <- ridge$fit(X_is, y_is)
#' s <- ridge$predict(m, X_oos)
#' }
#' @export
ml_make_model <- function(type = c("ridge","rf","xgboost","linear"),
                          params = list(), nrounds = 200L, ...) {
  type <- match.arg(type)

  if (type == "ridge") {
    stopifnot(requireNamespace("glmnet", quietly = TRUE))
    fit <- function(X, y) {
      X <- as.matrix(X); y <- as.numeric(y)
      keep <- is.finite(y) & rowSums(!is.finite(X)) == 0
      glmnet::cv.glmnet(X[keep, , drop = FALSE], y[keep], alpha = 0, nfolds = 3, ...)
    }
    pred <- function(m, X) {
      X <- as.matrix(X); X[!is.finite(X)] <- 0
      as.numeric(predict(m, newx = X, s = "lambda.min"))
    }
    return(list(fit = fit, predict = pred))
  }

  if (type == "rf") {
    stopifnot(requireNamespace("ranger", quietly = TRUE))
    fit <- function(X, y) {
      X <- as.data.frame(X); y <- as.numeric(y)
      keep <- is.finite(y) & rowSums(!is.finite(as.matrix(X))) == 0
      ranger::ranger(y = y[keep], x = X[keep, , drop = FALSE], ...)
    }
    pred <- function(m, X) {
      X <- as.data.frame(X); X[!is.finite(as.matrix(X))] <- 0
      as.numeric(predict(m, data = X)$predictions)
    }
    return(list(fit = fit, predict = pred))
  }

  if (type == "xgboost") {
    stopifnot(requireNamespace("xgboost", quietly = TRUE))
    fit <- function(X, y) {
      X <- as.matrix(X); y <- as.numeric(y)
      keep <- is.finite(y)  # XGBoost can handle NA in X; keep NA in X as missing
      d <- xgboost::xgb.DMatrix(data = X[keep, , drop = FALSE], label = y[keep], missing = NA)
      xgboost::xgb.train(params = params, data = d, nrounds = as.integer(nrounds), verbose = 0)
    }
    pred <- function(m, X) {
      X <- as.matrix(X); X[!is.finite(X)] <- NA
      as.numeric(predict(m, newdata = X, missing = NA))
    }
    return(list(fit = fit, predict = pred))
  }

  # linear baseline
  fit <- function(X, y) {
    X <- as.matrix(X); y <- as.numeric(y)
    keep <- is.finite(y) & rowSums(!is.finite(X)) == 0
    Xc <- cbind(1, X[keep, , drop = FALSE])
    list(coef = stats::lm.fit(Xc, y[keep])$coefficients)
  }
  pred <- function(m, X) {
    X <- as.matrix(X); X[!is.finite(X)] <- 0
    as.numeric(cbind(1, X) %*% m$coef)
  }
  list(fit = fit, predict = pred)
}

#' NA-tolerant ensemble blender (row-wise)
#'
#' Creates an equal- or user-weighted blend of multiple model objects produced
#' by [ml_make_model()]. The returned `fit` trains each component; `predict`
#' combines component predictions with an NA-safe weighted average.
#'
#' @param ... Two or more model objects each with `$fit`/`$predict`.
#' @param weights Optional numeric vector of blend weights (recycled).
#'
#' @return A list with `$fit` and `$predict` closures for the ensemble.
#' @examples
#' \dontrun{
#' ens <- ml_make_ensemble(ml_make_model("ridge"),
#'                         ml_make_model("rf"),
#'                         ml_make_model("xgboost"))
#' }
#' @export
ml_make_ensemble <- function(..., weights = NULL) {
  parts <- list(...)
  stopifnot(length(parts) >= 2L)
  if (is.null(weights)) weights <- rep(1/length(parts), length(parts))

  fit <- function(X, y) lapply(parts, function(w) w$fit(X, y))

  pred <- function(fitted, X) {
    pred_list <- lapply(seq_along(parts), function(i) {
      as.numeric(parts[[i]]$predict(fitted[[i]], X))  # force numeric vector
    })
    n <- length(pred_list[[1]])
    # sanity: all same length
    if (any(vapply(pred_list, length, 1L) != n)) {
      stop("Ensemble models returned different-length predictions.")
    }
    # build P (n x m)
    P <- matrix(NA_real_, nrow = n, ncol = length(pred_list))
    for (i in seq_along(pred_list)) P[, i] <- pred_list[[i]]

    # weighted average with NA masking
    W <- matrix(weights, nrow = n, ncol = length(weights), byrow = TRUE)
    mask <- is.finite(P)
    W_eff <- W * mask
    num <- rowSums(P * W_eff, na.rm = TRUE)
    den <- rowSums(W_eff, na.rm = TRUE)
    out <- num / den
    out[den == 0] <- NA_real_
    out
  }

  list(fit = fit, predict = pred)
}




#' Deterministic sequence model factory (GRU/LSTM/CNN1D with linear fallback)
#'
#' Returns `fit/predict` closures for sequence models that consume flattened
#' tabular inputs (n \eqn{\times}{x} (steps \eqn{\times}{x} p)) and internally reshape to (n, steps, p).
#' If Keras/TensorFlow is unavailable, falls back to a linear baseline so
#' examples remain runnable on CPU-only machines.
#'
#' Determinism knobs: fixed seeds, `TF_DETERMINISTIC_OPS=1`, no shuffle,
#' `workers=1`, and a fixed `pred_batch_size` to minimise retracing.
#'
#' @param type One of `"linear"`, `"gru"`, `"lstm"`, `"cnn1d"`.
#' @param steps Integer sequence length (e.g., 26 for 6 months of weeks).
#' @param units Hidden units for GRU/LSTM or filters for CNN1D.
#' @param dense Optional integer vector of additional dense layers.
#' @param dropout Dropout rate for recurrent/CNN core.
#' @param epochs,batch_size Training settings.
#' @param lr Learning rate.
#' @param patience Early-stopping patience.
#' @param seed Integer seed.
#' @param deterministic Logical; set determinism knobs when TRUE.
#' @param pred_batch_size Fixed batch size used at prediction time.
#'
#' @details
#' Optional dependencies: `keras` and `tensorflow`. When not available,
#' the factory returns the linear fallback.
#'
#' @return A list with `$fit` and `$predict` closures.
#' @examples
#' \dontrun{
#' seq_gru <- ml_make_seq_model("gru", steps = 26, units = 16, epochs = 12)
#' }
#' @export
ml_make_seq_model <- function(type = c("linear","gru","lstm","cnn1d"),
                              steps = 26L, units = 16L, dense = NULL,
                              dropout = 0.1, epochs = 12L, batch_size = 128L,
                              lr = 0.01, patience = 2L, seed = 123L,
                              deterministic = TRUE, pred_batch_size = 256L) {
  type <- match.arg(type)
  use_keras <- (type != "linear") &&
    requireNamespace("keras", quietly = TRUE) &&
    isTRUE(try(keras::is_keras_available(), silent = TRUE))

  # ---- Linear fallback (no Keras required) ----
  if (!use_keras) {
    fit  <- function(X, y) {
      Xc <- cbind(1, as.matrix(X))
      list(type = "lin", coef = stats::lm.fit(Xc, as.numeric(y))$coefficients)
    }
    pred <- function(m, Xnew) {
      Xc <- cbind(1, as.matrix(Xnew))
      as.numeric(Xc %*% m$coef)
    }
    return(list(fit = fit, predict = pred))
  }

  # ---- Keras models (GRU/LSTM/CNN1D) ----
  fit <- function(X, y) {
    set.seed(seed)

    # Determinism knobs
    if (deterministic) {
      Sys.setenv(PYTHONHASHSEED = as.character(seed))
      Sys.setenv(TF_DETERMINISTIC_OPS = "1")
      if (requireNamespace("tensorflow", quietly = TRUE)) {
        tensorflow::tf$random$set_seed(seed)
        try(tensorflow::tf$config$experimental$enable_op_determinism(TRUE), silent = TRUE)
      }
      keras::backend()$clear_session()
    }

    X <- as.matrix(X); y <- as.numeric(y)
    X[!is.finite(X)] <- 0; y[!is.finite(y)] <- 0
    p <- max(1L, as.integer(ncol(X) / steps))
    X3 <- array(X, dim = c(nrow(X), steps, p))

    inp <- keras::layer_input(shape = c(steps, p))
    core <- switch(type,
                   gru   = keras::layer_gru(units = units, dropout = dropout)(inp),
                   lstm  = keras::layer_lstm(units = units, dropout = dropout)(inp),
                   cnn1d = {
                     x <- keras::layer_conv_1d(filters = units, kernel_size = 3,
                                               padding = "causal", activation = "relu")(inp)
                     keras::layer_global_average_pooling_1d()(x)
                   }
    )

    head <- core
    if (!is.null(dense)) for (h in dense)
      head <- keras::layer_dense(units = h, activation = "relu")(head)
    out  <- keras::layer_dense(units = 1)(head)

    mdl <- keras::keras_model(inp, out)
    keras::compile(mdl, optimizer = keras::optimizer_adam(learning_rate = lr), loss = "mse")

    cb <- keras::callback_early_stopping(patience = patience, restore_best_weights = TRUE)

    keras::fit(
      mdl, X3, y,
      batch_size = batch_size, epochs = epochs,
      validation_split = 0.1, verbose = 0,
      callbacks = list(cb),
      shuffle = if (!deterministic) TRUE else FALSE,
      use_multiprocessing = if (!deterministic) TRUE else FALSE,
      workers = if (deterministic) 1L else 2L
    )

    list(type = type, model = mdl, p = p, steps = steps)
  }


  pred <- function(m, Xnew) {
    X <- as.matrix(Xnew); X[!is.finite(X)] <- 0
    X3 <- array(X, dim = c(nrow(X), m$steps, m$p))
    # Fixed inference batch size to reduce TF retracing warnings
    as.numeric(predict(m$model, X3, batch_size = as.integer(pred_batch_size), verbose = 0))
  }

  list(fit = fit, predict = pred)
}




#' Keep only rows with at least one finite score (rebalance dates)
#'
#' Convenience filter used by IC diagnostics to align to formation dates.
#'
#' @param scores_dt Wide scores table (`Date` + symbols) or long format coerced.
#' @return A `data.table` subset containing only rows with any finite score.
#' @keywords internal
#' @noRd
.ml_score_dates <- function(scores_dt) {
  DT <- data.table::as.data.table(scores_dt)
  syms <- setdiff(names(DT), "Date")
  keep <- rowSums(is.finite(as.matrix(DT[, ..syms]))) > 0
  DT[keep]
}


#' Rank-IC series computed on score (rebalance) dates
#'
#' Wrapper around `ic_series()` that first filters scores to formation dates
#' using an internal filter that keeps rows where at least one score is finite.
#'
#' @param scores_dt Scores table (wide).
#' @param labels_dt Labels table (wide).
#' @param method Correlation method; `"spearman"` (default) or `"pearson"`.
#'
#' @return A data frame/data.table with `Date` and IC values.
#' @examples
#' \dontrun{
#' ic <- ml_ic_series_on_scores(res_xgb$scores, Y, method = "spearman")
#' }
#' @export
ml_ic_series_on_scores <- function(scores_dt, labels_dt, method = "spearman") {
  ic_series(.ml_score_dates(scores_dt), labels_dt, method = method)
}



#' Rolling rank-IC plot (rebalance dates; leakage-safe)
#'
#' Computes the IC time series via [ml_ic_series_on_scores()] and plots the
#' rolling mean IC over a specified window. Returns the rolling statistics
#' invisibly for further inspection.
#'
#' @param scores_dt Scores (wide).
#' @param labels_dt Labels (wide).
#' @param window Integer window length (default 26).
#' @param method Correlation method; `"spearman"` (default) or `"pearson"`.
#'
#' @return (Invisibly) a data frame with `Date`, `roll_mean`, `roll_sd`, `roll_ICIR`.
#' @examples
#' \dontrun{
#' ris <- ml_plot_ic_roll(res_xgb$scores, Y, window = 8L)
#' }
#' @export
ml_plot_ic_roll <- function(scores_dt, labels_dt, window = 26L, method = "spearman") {
  ics <- ml_ic_series_on_scores(scores_dt, labels_dt, method = method)
  ris <- roll_ic_stats(ics, window = window)  # returns roll_mean, roll_sd, roll_ICIR
  plot(ris$Date, ris$roll_mean, type = "l",
       main = paste0("Rolling IC (", window, "w)"),
       xlab = "", ylab = "IC mean")
  invisible(ris)
}



#' Run multi-horizon ML backtests (pooled or sector-neutral)
#'
#' Convenience wrapper around [PortfolioTesteR::ml_backtest()] that repeats the
#' same specification across multiple horizons, returning a named list of
#' backtest objects keyed as `"H1w"`, `"H4w"`, `"H13w"`, etc.
#'
#' @param features_list Named list of data.tables with factor scores (each with a
#'   `Date` column and one column per symbol). Typically from
#'   [PortfolioTesteR::ml_prepare_features()].
#' @param prices_weekly Wide price table (weekly) with `Date` + one column per
#'   symbol (adjusted prices). Used both to create labels and run the backtest.
#' @param horizons Integer vector of horizons in *weeks* (e.g., `c(1L,4L,13L)`).
#' @param fit_fn,predict_fn Model fit/predict closures as returned by
#'   [PortfolioTesteR::ml_make_model()].
#' @param schedule Walk-forward schedule list with elements `is`, `oos`, `step`.
#' @param transform Feature transform (default `"zscore"`). Passed to
#'   [PortfolioTesteR::ml_backtest()].
#' @param selection List describing selection rules (e.g., `list(top_k=20L, max_per_group=3L)`).
#' @param weighting List describing weighting rules (e.g., `list(method="softmax", temperature=12)`).
#' @param caps List with exposure caps (e.g., `list(max_per_symbol=0.10, max_per_group=NULL)`).
#' @param group_mode `"pooled"` or `"per_group"`. If `"per_group"`, you must pass `group_map`.
#' @param group_map A two-column table with columns `Symbol` and `Group` defining
#'   the grouping (e.g., sectors) for `"per_group"` mode.
#' @param initial_capital Numeric. Starting capital for the backtest (default `1e5`).
#' @param name_prefix Optional string prefixed to each backtest title.
#' @param seed Optional integer. If provided, the same seed is set before each
#'   horizon’s backtest call to ensure deterministic tie-breaks.
#' @param ... Additional arguments forwarded to [PortfolioTesteR::ml_backtest()]
#'   (kept for future extensibility; no effect if unused).
#'
#' @return A named list of backtest objects (as returned by
#'   [PortfolioTesteR::ml_backtest()]), with names like `"H1w"`, `"H4w"`, … .
#'
#' @details
#' This function does **not** change core behavior; it only removes boilerplate
#' when running identical specs across horizons and (optionally) grouping
#' regimes. It preserves all defaults you pass for selection, weighting,
#' transforms, caps, and schedule.
#'
#' @examples
#' \donttest{
#' library(PortfolioTesteR)
#' data(sample_prices_weekly, package = "PortfolioTesteR")
#'
#' # Minimal features for the example
#' X <- ml_prepare_features(
#'   prices_weekly = sample_prices_weekly,
#'   include = c("mom12","mom26")
#' )
#'
#' # Simple deterministic model
#' model <- ml_make_model("linear")
#' sched <- list(is = 156L, oos = 4L, step = 4L)
#'
#' set.seed(42)
#' bt_list <- ml_backtest_multi(
#'   features_list = X,
#'   prices_weekly = sample_prices_weekly,
#'   horizons = c(1L, 4L),
#'   fit_fn = model$fit,
#'   predict_fn = model$predict,
#'   schedule = sched,
#'   selection = list(top_k = 5L),
#'   weighting = list(method = "softmax", temperature = 12),
#'   caps = list(max_per_symbol = 0.10),
#'   group_mode = "pooled",
#'   name_prefix = "Demo ",
#'   seed = 42
#' )
#'
#' names(bt_list)   # "H1w" "H4w"
#' }
#'
#' @export
ml_backtest_multi <- function(features_list, prices_weekly, horizons,
                              fit_fn, predict_fn, schedule,
                              transform = "zscore",
                              selection = list(top_k = 20L),
                              weighting = list(method = "softmax", temperature = 12),
                              caps = list(max_per_symbol = 0.10),
                              group_mode = c("pooled","per_group"),
                              group_map = NULL,
                              initial_capital = 1e5,
                              name_prefix = "",
                              seed = NULL,
                              ...) {
  group_mode <- match.arg(group_mode)

  # name vector: c("H1w","H4w",...)
  keys <- paste0("H", horizons, "w")

  res <- stats::setNames(vector("list", length(horizons)), keys)

  for (i in seq_along(horizons)) {
    h <- horizons[i]
    if (!is.null(seed)) set.seed(seed)

    Yh <- PortfolioTesteR::make_labels(
      prices  = prices_weekly,
      horizon = h,
      type    = "log"
    )

    args_common <- list(
      features_list   = features_list,
      labels          = Yh,
      fit_fn          = fit_fn,
      predict_fn      = predict_fn,
      schedule        = schedule,
      transform       = transform,
      selection       = selection,
      weighting       = weighting,
      caps            = caps,
      prices          = prices_weekly,
      initial_capital = initial_capital
    )

    res[[i]] <-
      if (group_mode == "pooled") {
        do.call(PortfolioTesteR::ml_backtest,
                c(args_common,
                  list(group = "pooled",
                       name  = sprintf("%sPooled (%dw)", name_prefix, h)),
                  list(...)))
      } else {
        stopifnot(!is.null(group_map),
                  all(c("Symbol","Group") %in% names(group_map)))
        do.call(PortfolioTesteR::ml_backtest,
                c(args_common,
                  list(group = "per_group", group_map = group_map,
                       name  = sprintf("%sSector-neutral (%dw)", name_prefix, h)),
                  list(...)))
      }
  }

  # Attach a tiny config breadcrumb for reproducibility
  attr(res, "config") <- list(
    horizons = horizons,
    schedule = schedule,
    selection = selection,
    weighting = weighting,
    caps = caps,
    transform = transform,
    group_mode = group_mode,
    seed = seed,
    name_prefix = name_prefix
  )

  res
}


#' Collect diagnostics from two `ml_backtest_multi()` runs
#'
#' Builds a compact set of outputs—coverage, IC series, OOS-only rolling IC,
#' performance tables (Full/Pre/Post), turnover, and a cost sweep—given two
#' lists of backtests (pooled and per-group) produced by `ml_backtest_multi()`.
#'
#' @param bt_pooled_list Named list of backtests (keys like `H4w`, `H13w`)
#'   produced with `group_mode = "pooled"`.
#' @param bt_neutral_list Named list of backtests (same keys) produced with
#'   `group_mode = "per_group"`.
#' @param weekly_prices Deprecated alias for `prices`; kept for backwards compatibility.
#' @param horizons Integer vector of horizons (in weeks) expected in the lists.
#' @param split_date `Date` used to split performance into `Pre`/`Post`.
#' @param cost_bps Numeric vector of per-rebalance cost levels (in basis points)
#'   for the turnover-based cost sweep. Default `c(5, 10, 15, 20, 25, 50)`.
#' @param freq Integer frequency used by `perf_metrics()` (e.g., `52` for weekly).
#' @param prices Optional price table (preferred). If `NULL`, `weekly_prices` is used.
#' @param ic_roll_window Integer window length (weeks) for rolling IC on OOS decision dates.
#'   Default `26L`.
#' @param mask_scores_to_decision_dates Logical; if `TRUE` (default) scores are
#'   masked to OOS decision dates only (see [scores_oos_only()]).
#' @param cost_model Function `(turnover, bps)` returning per-period cost to subtract
#'   from returns in the sweep. Default scales linearly with turnover.
#'
#' @details
#' Both input lists must have identical horizon keys (`paste0("H", h, "w")`).
#' Coverage and IC series are computed from stored `scores`; rolling IC is built
#' on OOS decision dates only; performance is summarised for the full sample
#' and `Pre`/`Post` relative to `split_date`; turnover is derived from realised
#' sector-neutral weights; and a turnover-based cost sweep is evaluated on the
#' sector-neutral run across `cost_bps`.
#'
#' @return
#' A named list with one element per horizon, each containing:
#' \itemize{
#'   \item `bt_pooled`, `bt_neutral` — the input backtests;
#'   \item `coverage` — coverage by date for pooled/neutral;
#'   \item `ic_series` — raw IC series for pooled/neutral;
#'   \item `icroll_oos_26w` — rolling IC (OOS-only) for pooled/neutral;
#'   \item `masked_scores` — OOS-masked score tables for pooled/neutral;
#'   \item `perf_tables` — performance tables (Full/Pre/Post);
#'   \item `turnover_neutral` — turnover series for the sector-neutral run;
#'   \item `cost_sweep_neutral` — performance under gross/net across `cost_bps`.
#' }
#'
#' @seealso [ml_backtest_multi()], [scores_oos_only()], [PortfolioTesteR::perf_metrics()]
#' @family Chapter3-helpers
#' @examples
#' \donttest{
#' if (requireNamespace("PortfolioTesteR", quietly = TRUE)) {
#'   library(PortfolioTesteR)
#'   data(sample_prices_weekly, package = "PortfolioTesteR")
#'   # Simple feature
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)
#'   feats <- list(mom12 = mom12)
#'   fit_first     <- function(X, y, ...) list()
#'   predict_first <- function(model, Xnew, ...) as.numeric(Xnew[[1]])
#'   sch <- list(is = 52L, oos = 26L, step = 26L)
#'   syms <- setdiff(names(sample_prices_weekly), "Date")
#'   gmap <- data.frame(Symbol = syms,
#'                      Group  = rep(c("G1","G2"), length.out = length(syms)))
#'   bt_pooled  <- ml_backtest_multi(feats, sample_prices_weekly, c(4L),
#'                                   fit_first, predict_first, sch,
#'                                   selection = list(top_k = 5L),
#'                                   weighting = list(method = "softmax", temperature = 12),
#'                                   caps = list(max_per_symbol = 0.10),
#'                                   group_mode = "pooled")
#'   bt_neutral <- ml_backtest_multi(feats, sample_prices_weekly, c(4L),
#'                                   fit_first, predict_first, sch,
#'                                   selection = list(top_k = 5L),
#'                                   weighting = list(method = "softmax", temperature = 12),
#'                                   caps = list(max_per_symbol = 0.10),
#'                                   group_mode = "per_group",
#'                                   group_map = gmap)
#'   out <- pt_collect_results(
#'     bt_pooled_list  = bt_pooled,
#'     bt_neutral_list = bt_neutral,
#'     prices          = sample_prices_weekly,
#'     horizons        = c(4L),
#'     split_date      = as.Date("2019-01-01"),
#'     cost_bps        = c(5, 15),
#'     freq            = 52,
#'     ic_roll_window  = 13L
#'   )
#'   names(out)
#'   str(out[["H4w"]]$perf_tables)
#' }}
#'
#' @export
#' @importFrom data.table as.data.table setorder shift melt fifelse rbindlist
#' @importFrom stats setNames
pt_collect_results <- function(bt_pooled_list,
                               bt_neutral_list,
                               weekly_prices,              # kept for backwards compat
                               horizons,
                               split_date,
                               cost_bps = c(5,10,15,20,25,50),
                               freq = 52,
                               # --- new, all default to current behaviour ---
                               prices = NULL,              # alias for weekly_prices
                               ic_roll_window = 26L,
                               mask_scores_to_decision_dates = TRUE,
                               cost_model = function(turnover, bps) (bps/10000) * turnover) {

  prices_in <- if (!is.null(prices)) prices else weekly_prices

  # ---- internal helpers ----
  scores_oos_only_ <- function(scores_dt, weights_wide) {
    if (!mask_scores_to_decision_dates) return(data.table::as.data.table(scores_dt))
    scores_oos_only(scores_dt, weights_wide)  # uses the exported helper above
  }

  turnover_series_ <- function(weights_wide) {
    W <- data.table::as.data.table(weights_wide)
    if (!nrow(W)) return(data.table::data.table(Date = as.Date(character()), turnover = numeric()))
    if (is.numeric(W$Date)) W[, Date := as.Date(Date, origin = "1970-01-01")]
    sy <- setdiff(names(W), "Date")
    if (!length(sy)) return(data.table::data.table(Date = as.Date(character()), turnover = numeric()))
    Wn <- W[rowSums(as.matrix(W[, ..sy]) != 0) > 0]
    if (!nrow(Wn)) return(data.table::data.table(Date = as.Date(character()), turnover = numeric()))
    L  <- data.table::melt(Wn, id.vars = "Date", variable.name = "Symbol", value.name = "w")
    data.table::setorder(L, Date, Symbol)
    L[, w_lag := data.table::shift(w, 1L), by = Symbol]
    L[, .(turnover = 0.5 * sum(abs(w - data.table::fifelse(is.na(w_lag), 0, w_lag)))), by = Date][!is.na(turnover)]
  }

  as_dt_metrics_ <- function(M) {
    if (data.table::is.data.table(M)) return(M)
    if (is.data.frame(M))            return(data.table::as.data.table(M))
    if (is.list(M))                  return(data.table::as.data.table(M))
    data.table::as.data.table(as.list(M))
  }

  perf_splits_ <- function(bt, split_date, freq) {
    b <- bt$backtest
    dates <- as.Date(b$dates, origin = "1970-01-01")
    rvec  <- if (!is.null(b$returns)) as.numeric(b$returns) else as.numeric(b$active_returns)
    if (!length(rvec)) return(data.table::data.table(Slice = character()))
    R <- data.table::data.table(Date = dates[seq_along(rvec)], Ret = rvec)

    summarize <- function(D) {
      M <- PortfolioTesteR::perf_metrics(data.frame(Date = D$Date, ret = D$Ret), freq = freq)
      as_dt_metrics_(M)
    }
    full <- summarize(R);                      full[,  Slice := "Full"]
    pre  <- summarize(R[Date <  split_date]);  pre[,   Slice := "Pre"]
    post <- summarize(R[Date >= split_date]);  post[,  Slice := "Post"]
    data.table::rbindlist(list(full, pre, post), use.names = TRUE, fill = TRUE)
  }

  cost_sweep_ <- function(returns_dt, turnover_dt, cost_bps, freq) {
    if (!nrow(returns_dt)) return(data.table::data.table())
    D <- base::merge(data.table::copy(returns_dt),
                     data.table::copy(turnover_dt),
                     by = "Date", all.x = TRUE)
    D[is.na(turnover), turnover := 0]

    out <- lapply(cost_bps, function(bps) {
      D[, Ret_net := Ret - cost_model(turnover, bps)]
      g <- PortfolioTesteR::perf_metrics(data.frame(Date = D$Date, ret = D$Ret),      freq = freq)
      n <- PortfolioTesteR::perf_metrics(data.frame(Date = D$Date, ret = D$Ret_net),  freq = freq)
      gdt <- as_dt_metrics_(g)[, `:=`(Cost_bps = bps, Net = FALSE)]
      ndt <- as_dt_metrics_(n)[, `:=`(Cost_bps = bps, Net = TRUE)]
      list(gdt, ndt)
    })
    data.table::rbindlist(unlist(out, recursive = FALSE), use.names = TRUE, fill = TRUE)
  }

  # ---- sanity ----
  need_names <- paste0("H", horizons, "w")
  stopifnot(all(need_names %in% names(bt_pooled_list)),
            all(need_names %in% names(bt_neutral_list)))

  out <- stats::setNames(vector("list", length(horizons)), need_names)

  # ---- loop ----
  for (h in horizons) {
    key <- paste0("H", h, "w")
    btP <- bt_pooled_list[[key]]
    btN <- bt_neutral_list[[key]]

    Yh <- PortfolioTesteR::make_labels(prices = prices_in, horizon = h, type = "log")

    covP <- if (nrow(btP$scores)) PortfolioTesteR::coverage_by_date(btP$scores) else data.table::data.table()
    covN <- if (nrow(btN$scores)) PortfolioTesteR::coverage_by_date(btN$scores) else data.table::data.table()
    icP  <- if (nrow(btP$scores)) PortfolioTesteR::ml_ic_series_on_scores(btP$scores, Yh) else data.table::data.table()
    icN  <- if (nrow(btN$scores)) PortfolioTesteR::ml_ic_series_on_scores(btN$scores, Yh) else data.table::data.table()

    S_P_oos <- scores_oos_only_(btP$scores, btP$weights)
    S_N_oos <- scores_oos_only_(btN$scores, btN$weights)

    icrollP_oos <- if (nrow(S_P_oos)) PortfolioTesteR::ml_plot_ic_roll(S_P_oos, Yh, window = ic_roll_window) else data.table::data.table()
    icrollN_oos <- if (nrow(S_N_oos)) PortfolioTesteR::ml_plot_ic_roll(S_N_oos, Yh, window = ic_roll_window) else data.table::data.table()

    perfP <- perf_splits_(btP, split_date = split_date, freq = freq)
    perfN <- perf_splits_(btN, split_date = split_date, freq = freq)

    turnsN <- turnover_series_(btN$weights)

    bN <- btN$backtest
    datesN <- as.Date(bN$dates, origin = "1970-01-01")
    rvecN  <- if (!is.null(bN$returns)) as.numeric(bN$returns) else as.numeric(bN$active_returns)
    retN   <- if (length(rvecN)) data.table::data.table(Date = datesN[seq_along(rvecN)], Ret = rvecN) else data.table::data.table()

    costsN <- cost_sweep_(retN, turnsN, cost_bps = cost_bps, freq = freq)

    out[[key]] <- list(
      bt_pooled            = btP,
      bt_neutral           = btN,
      coverage             = list(pooled = covP, neutral = covN),
      ic_series            = list(pooled = icP, neutral = icN),
      icroll_oos_26w       = list(pooled = icrollP_oos, neutral = icrollN_oos),
      masked_scores        = list(pooled = S_P_oos, neutral = S_N_oos),
      perf_tables          = list(pooled = perfP, neutral = perfN),
      turnover_neutral     = turnsN,
      cost_sweep_neutral   = costsN
    )
  }

  attr(out, "config") <- list(
    horizons           = horizons,
    split_date         = split_date,
    cost_bps           = cost_bps,
    freq               = freq,
    ic_roll_window     = ic_roll_window,
    masked_to_decision = mask_scores_to_decision_dates
  )
  out
}



#' Mask score tables to out-of-sample decision dates
#'
#' Utility used in the chapter’s diagnostics: keep scores only on dates when a
#' portfolio decision was actually made (non-zero realised weights); set other
#' dates to `NA`. Inputs are wide by symbol with a `Date` column.
#'
#' @param scores_dt Wide table of model scores with columns `Date`, `SYM1`,
#'   `SYM2`, … (one column per symbol).
#' @param weights_wide Wide table of realised portfolio weights with columns
#'   `Date`, `SYM1`, `SYM2`, … (one column per symbol). Decision dates are
#'   inferred as rows where any symbol weight is non-zero.
#'
#' @return
#' A copy of `scores_dt` where rows not matching decision dates are set to `NA`
#' (except the `Date` column). If either input is empty, returns `scores_dt[0]`.
#'
#' @examples
#' # Toy example
#' dates <- as.Date("2020-01-01") + 7*(0:5)
#' scores <- data.frame(
#'   Date = dates,
#'   AAA = seq(0.1, 0.6, length.out = 6),
#'   BBB = rev(seq(0.1, 0.6, length.out = 6))
#' )
#' weights <- data.frame(
#'   Date = dates,
#'   AAA  = c(0, 0.1, 0, 0.2, 0, 0.15),
#'   BBB  = c(0, 0,   0, 0,   0, 0   )
#' )
#' scores_oos_only(scores, weights)
#'
#' @seealso [pt_collect_results()]
#' @family Chapter3-helpers
#' @export
#' @importFrom data.table as.data.table setorder shift melt fifelse
scores_oos_only <- function(scores_dt, weights_wide) {
  S <- data.table::as.data.table(scores_dt)
  W <- data.table::as.data.table(weights_wide)
  if (!nrow(S) || !nrow(W)) return(S[0])
  if (is.numeric(S$Date)) S[, Date := as.Date(Date, origin = "1970-01-01")]
  if (is.numeric(W$Date)) W[, Date := as.Date(Date, origin = "1970-01-01")]
  sy <- setdiff(names(W), "Date")
  if (!length(sy)) return(S[0])
  ddates <- W$Date[rowSums(as.matrix(W[, ..sy]) != 0) > 0]
  if (!length(ddates)) return(S[0])
  S[!(Date %in% ddates), (setdiff(names(S), "Date")) := NA][]
}

