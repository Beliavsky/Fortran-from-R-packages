
# ---- Data & leakage control -------------------------------------------------
utils::globalVariables(c("Date"))
# Internal: rough memory estimate for seq designs (to avoid freezing RStudio)
.seq_mem_warn <- function(n_samples, steps, p, max_bytes = 1e8) {
  est_bytes <- as.numeric(n_samples) * as.numeric(steps) * as.numeric(p) * 8  # 8 bytes per double
  if (!is.finite(est_bytes)) return(invisible(est_bytes))
  if (est_bytes > max_bytes) {
    warning(sprintf(
      "Sequence design may be large: ~%.1f MB (samples=%s, steps=%s, features=%s).",
      est_bytes / 1e6, n_samples, steps, p
    ))
  }
  invisible(est_bytes)
}

#' Join multiple panels on intersecting dates (unique symbol names)
#'
#' Align and join a list of wide panels on their common dates. Each input
#' panel must have one \code{Date} column and a disjoint set of symbol columns.
#'
#' @details
#' All panels are first aligned to the intersection of their \code{Date} values.
#' Symbol names must be unique across panels; if the same symbol appears in
#' multiple inputs, an error is raised.
#'
#' @param panels List of wide panels (\code{data.frame} or \code{data.table}), each with columns \code{Date} + symbols.
#'
#' @return A \code{data.table} with \code{Date} plus the union of all symbol columns, restricted to common dates.
#' @export
join_panels <- function(panels) {
  stopifnot(is.list(panels), length(panels) >= 1)
  as_dt <- function(x) { if (data.table::is.data.table(x)) x else data.table::as.data.table(x) }
  clean_one <- function(dt) {
    dt <- as_dt(dt)
    if (!"Date" %in% names(dt)) stop("Each panel must have a 'Date' column")
    dups <- setdiff(names(dt)[duplicated(names(dt))], "Date")
    if (length(dups)) dt <- dt[, !dups, with = FALSE]
    dt
  }
  L <- lapply(panels, clean_one)
  out <- Reduce(function(a,b) merge(a, b, by = "Date", all = FALSE), L)
  data.table::setorderv(out, "Date")
  out
}

#' Lag each symbol column by k steps
#'
#' Given a wide panel with a `Date` column followed by symbol columns, returns
#' the same shape with each symbol column lagged by `k` periods. The `Date`
#' column is preserved; leading values introduced by the lag are `NA_real_`.
#'
#' @param df data.frame or data.table with columns `Date` then symbols.
#' @param k Integer lag (>= 1).
#'
#' @return A `data.table` with the same columns as `df`, lagged by `k`.
#' @export
#' @examples
#' x <- data.frame(Date = as.Date("2020-01-01") + 0:2, A = 1:3, B = 11:13)
#' panel_lag(x, 1L)
panel_lag <- function(df, k = 1L) {
  dt <- if (data.table::is.data.table(df)) data.table::copy(df) else data.table::as.data.table(df)
  if (!"Date" %in% names(dt)) stop("'df' must have a Date column")
  syms <- setdiff(names(dt), "Date")
  if (k <= 0L) return(dt)
  n <- nrow(dt)
  for (nm in syms) {
    x <- dt[[nm]]
    if (is.factor(x)) x <- as.numeric(x)
    dt[[nm]] <- c(rep(NA_real_, k), head(as.numeric(x), n - k))
  }
  dt
}


#' Make future-return labels aligned to the decision date
#'
#' Compute forward returns over a fixed horizon and align them to the decision date.
#'
#' @details
#' For each date \eqn{t}, the label is computed from prices at \eqn{t} and \eqn{t + h}.
#' \itemize{
#'   \item \code{type = "simple"}: \eqn{p_{t+h}/p_t - 1}
#'   \item \code{type = "log"}: \eqn{\log(p_{t+h}) - \log(p_t)}
#'   \item \code{type = "sign"}: \code{sign(simple return)}
#' }
#' Trailing dates that do not have \code{horizon} steps ahead are set to \code{NA}.
#'
#' @param prices Wide price panel (\code{Date} + symbols).
#' @param horizon Integer \code{>= 1}, number of forward steps for the label.
#' @param type Character, one of \code{"log"}, \code{"simple"}, \code{"sign"}.
#'
#' @return A \code{data.table} with \code{Date} + symbols containing the labels.
#' @export
make_labels <- function(prices, horizon = 4L, type = c("log","simple","sign")) {
  type <- match.arg(type)
  dt <- if (data.table::is.data.table(prices)) data.table::copy(prices) else data.table::as.data.table(prices)
  if (!"Date" %in% names(dt)) stop("'prices' must have a Date column")
  syms <- setdiff(names(dt), "Date")
  n <- nrow(dt)
  lead_h <- function(x, h) c(x[(h+1):length(x)], rep(NA_real_, h))
  out <- data.table::data.table(Date = dt$Date)
  for (nm in syms) {
    p <- as.numeric(dt[[nm]])
    fwd <- lead_h(p, horizon)
    lab <- switch(type,
                  log    = log(fwd) - log(p),
                  simple = fwd / p - 1,
                  sign   = as.numeric((fwd / p - 1) > 0)
    )
    out[[nm]] <- lab
  }
  out
}


#' Quick leakage guard: date alignment & NA expectations
#'
#' @param features Wide feature panel with a `Date` column.
#' @param labels   Wide label panel (same `Date` index / symbols as features).
#' @param verbose  If TRUE, prints diagnostic messages.
#'
#' @return TRUE/FALSE (invisible), with messages if `verbose = TRUE`.
#' @export
validate_no_leakage <- function(features, labels, verbose = TRUE) {
  fx <- if (data.table::is.data.table(features)) features else data.table::as.data.table(features)
  ly <- if (data.table::is.data.table(labels))   labels   else data.table::as.data.table(labels)
  if (!identical(fx$Date, ly$Date)) {
    if (verbose) message("Date misalignment between features and labels")
    return(FALSE)
  }
  syms <- intersect(setdiff(names(fx), "Date"), setdiff(names(ly), "Date"))
  ok <- TRUE; msgs <- character()
  # Tail of labels should have NAs (future unknown)
  tail_rows <- max(1L, nrow(ly) - 3L):nrow(ly)
  tail_all_na <- all(colSums(is.na(ly[tail_rows, ..syms])) > 0)
  if (!tail_all_na) { ok <- FALSE; msgs <- c(msgs, "Labels tail not NA (horizon?)") }
  # Features should have at least 1-period NA head if you lagged
  head_rows <- 1L:min(5L, nrow(fx))
  head_has_na <- any(colSums(is.na(fx[head_rows, ..syms])) > 0)
  if (!head_has_na) { ok <- FALSE; msgs <- c(msgs, "Features head has no NA (did you panel_lag()?)") }
  if (verbose && length(msgs)) message(paste(msgs, collapse = " | "))
  ok
}

.roll_align_common <- function(panels) {
  stopifnot(is.list(panels), length(panels) >= 1)

  coerce_dt <- function(x) {
    dx <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    nm <- names(dx)
    if (!("Date" %in% nm)) {
      cand <- which(tolower(nm) == "date")
      if (!length(cand)) cand <- 1L
      data.table::setnames(dx, cand[1L], "Date")
    }
    if (!inherits(dx[["Date"]], "Date")) {
      if (is.numeric(dx[["Date"]])) dx[, Date := as.Date(Date, origin = "1970-01-01")] else dx[, Date := as.Date(Date)]
    }
    data.table::setorderv(dx, "Date")
    dx
  }

  L <- lapply(panels, coerce_dt)

  symsets <- lapply(L, function(d) setdiff(names(d), "Date"))
  syms <- Reduce(intersect, symsets)
  if (!length(syms)) stop("No common symbol columns across panels.")

  dates <- Reduce(intersect, lapply(L, function(d) d[["Date"]]))
  if (!length(dates)) stop("No common dates across panels.")

  L2 <- lapply(L, function(z) {
    z <- z[z[["Date"]] %in% dates]
    z <- z[, c("Date", syms), with = FALSE]
    data.table::setorderv(z, "Date")
    z
  })

  list(panels = L2, dates = dates, symbols = syms)
}


#' Rolling fit/predict for tabular features (pooled / per-symbol / per-group)
#'
#' @param features_list list of wide panels (each: `Date` + symbols).
#' @param label wide panel of future returns (same symbol set as features).
#' @param fit_fn function `(X, y) -> model` trained on in-sample stacked rows.
#' @param predict_fn function `(model, Xnew) -> numeric` scores.
#' @param is_periods,oos_periods,step ints; in-sample length, out-of-sample length,
#'   and step size for the rolling window.
#' @param group one of `"pooled"`, `"per_symbol"`, `"per_group"`.
#' @param group_map optional `data.frame(Symbol, Group)` if `group = "per_group"`.
#' @param na_action `"omit"` or `"zero"` for feature NA handling.
#' @return wide panel of scores: `Date + symbols`.
#' @export
#' @examples
#' \donttest{
#' data(sample_prices_weekly); data(sample_prices_daily)
#' mom <- panel_lag(calc_momentum(sample_prices_weekly, 12), 1)
#' vol <- panel_lag(align_to_timeframe(
#'   calc_rolling_volatility(sample_prices_daily, 20),
#'   sample_prices_weekly$Date, "forward_fill"), 1)
#' Y <- make_labels(sample_prices_weekly, horizon = 4, type = "log")
#' fit_lm  <- function(X,y){ Xc <- cbind(1,X); list(coef=stats::lm.fit(Xc,y)$coefficients) }
#' pred_lm <- function(m,X){ as.numeric(cbind(1,X) %*% m$coef) }
#' S <- roll_fit_predict(list(mom=mom, vol=vol), Y, fit_lm, pred_lm, 52, 4, 4)
#' head(S)
#' }
roll_fit_predict <- function(features_list, label, fit_fn, predict_fn,
                             is_periods = 156L, oos_periods = 4L, step = 4L,
                             group = c("pooled","per_symbol","per_group"),
                             group_map = NULL,
                             na_action = c("omit","zero")) {
  if (data.table::is.data.table(features_list)) features_list <- list(features_list)
  group <- match.arg(group)
  na_action <- match.arg(na_action)
  stopifnot(is.list(features_list), length(features_list) >= 1)

  aligned <- .roll_align_common(c(features_list, list(label)))
  dates   <- aligned$dates
  syms    <- aligned$symbols
  p       <- length(features_list)

  L <- aligned$panels
  lab <- L[[length(L)]]
  feats <- L[seq_len(length(L) - 1L)]

  if (group == "per_group") {
    if (is.null(group_map)) stop("group_map is required when group='per_group'")
    group_map <- data.table::as.data.table(group_map)
    if (!all(c("Symbol","Group") %in% names(group_map))) stop("group_map must have Symbol, Group")
    sym2grp <- setNames(group_map$Group, group_map$Symbol)
  }

  n <- length(dates)
  S <- data.table::data.table(Date = dates)
  for (nm in syms) S[[nm]] <- NA_real_

  get_x <- function(t, nm) {
    v <- numeric(p)
    for (j in seq_len(p)) v[j] <- as.numeric(feats[[j]][[nm]][t])
    if (any(!is.finite(v))) {
      if (na_action == "omit") return(NULL)
      if (na_action == "zero") v[!is.finite(v)] <- 0
    }
    v
  }

  i <- 1L
  while (TRUE) {
    is_start <- i
    is_end   <- i + is_periods - 1L
    oos_end  <- is_end + oos_periods
    if (oos_end > n) break

    if (group == "pooled") {
      Xis <- NULL; yis <- NULL
      for (t in is_start:is_end) {
        for (nm in syms) {
          v <- get_x(t, nm)
          y <- as.numeric(lab[[nm]][t])
          if (!is.null(v) && is.finite(y)) {
            Xis <- if (is.null(Xis)) matrix(v, nrow = 1L) else rbind(Xis, v)
            yis <- c(yis, y)
          }
        }
      }
      model <- if (!is.null(Xis) && length(yis)) fit_fn(Xis, yis) else NULL

    } else if (group == "per_symbol") {
      model <- vector("list", length(syms)); names(model) <- syms
      for (nm in syms) {
        Xis <- NULL; yis <- NULL
        for (t in is_start:is_end) {
          v <- get_x(t, nm); y <- as.numeric(lab[[nm]][t])
          if (!is.null(v) && is.finite(y)) {
            Xis <- if (is.null(Xis)) matrix(v, nrow = 1L) else rbind(Xis, v)
            yis <- c(yis, y)
          }
        }
        if (!is.null(Xis) && length(yis)) model[[nm]] <- fit_fn(Xis, yis)
      }

    } else { # per_group
      groups <- unique(group_map$Group)
      model <- vector("list", length(groups)); names(model) <- groups
      for (g in groups) {
        sy_g <- intersect(syms, group_map$Symbol[group_map$Group == g])
        if (!length(sy_g)) next
        Xis <- NULL; yis <- NULL
        for (t in is_start:is_end) {
          for (nm in sy_g) {
            v <- get_x(t, nm); y <- as.numeric(lab[[nm]][t])
            if (!is.null(v) && is.finite(y)) {
              Xis <- if (is.null(Xis)) matrix(v, nrow = 1L) else rbind(Xis, v)
              yis <- c(yis, y)
            }
          }
        }
        if (!is.null(Xis) && length(yis)) model[[g]] <- fit_fn(Xis, yis)
      }
    }

    for (t in (is_end + 1L):oos_end) {
      for (nm in syms) {
        v <- get_x(t, nm); if (is.null(v)) next
        sc <- NA_real_
        if (group == "pooled") {
          if (!is.null(model)) sc <- as.numeric(predict_fn(model, matrix(v, nrow = 1L)))
        } else if (group == "per_symbol") {
          m <- model[[nm]]; if (!is.null(m)) sc <- as.numeric(predict_fn(m, matrix(v, nrow = 1L)))
        } else {
          g <- sym2grp[[nm]]; m <- model[[g]]
          if (!is.null(m)) sc <- as.numeric(predict_fn(m, matrix(v, nrow = 1L)))
        }
        S[[nm]][t] <- sc
      }
    }
    i <- i + step
  }
  S
}



# ---- Score transforms -------------------------------------------------------
#' Per-date score transform (z-score or rank)
#' @param scores wide panel `Date + symbols`.
#' @param method `"zscore"` or `"rank"`.
#' @param per_date logical; currently must be TRUE.
#' @param robust logical; robust z-score via median/MAD.
#' @return panel of transformed scores.
#' @export
transform_scores <- function(scores, method = c("zscore","rank"), per_date = TRUE, robust = FALSE) {
  method <- match.arg(method)
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)
  syms <- setdiff(names(S), "Date")
  stopifnot(per_date)
  for (r in seq_len(nrow(S))) {
    v <- as.numeric(S[r, ..syms]); ok <- is.finite(v)
    if (!any(ok)) next
    if (method == "zscore") {
      if (isTRUE(robust)) {
        m <- stats::median(v[ok]); s <- stats::mad(v[ok], constant = 1.4826)
      } else {
        m <- mean(v[ok]); s <- stats::sd(v[ok])
      }
      if (!is.finite(s) || s == 0) { v[ok] <- 0 } else { v[ok] <- (v[ok] - m)/s }
    } else {
      rk <- rank(v[ok], ties.method = "average")
      v[ok] <- (rk - min(rk))/(max(rk) - min(rk))
    }
    S[r, (syms) := as.list(v)]
  }
  S
}

# Explicit ranker
rank_scores <- function(scores, ties = c("average","first","min","max"), normalize = TRUE) {
  ties <- match.arg(ties)
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)
  syms <- setdiff(names(S), "Date")
  for (r in seq_len(nrow(S))) {
    v <- as.numeric(S[r, ..syms]); ok <- is.finite(v)
    if (!any(ok)) next
    rk <- rank(v[ok], ties.method = ties)
    if (normalize) rk <- (rk - min(rk))/(max(rk) - min(rk))
    v[ok] <- rk; S[r, (syms) := as.list(v)]
  }
  S
}


#' Select top-K scores per date
#' @param scores score panel.
#' @param k integer; how many to keep per date.
#' @param ties Ties method passed to base::rank (e.g., 'first','average').
#' @return logical mask panel (`Date + symbols`) marking selected names.
#' @export
select_top_k_scores <- function(scores, k, ties = "first") {
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)
  syms <- setdiff(names(S), "Date")
  out <- data.table::data.table(Date = S$Date)
  for (nm in syms) out[[nm]] <- FALSE
  for (r in seq_len(nrow(S))) {
    v <- as.numeric(S[r, ..syms]); ok <- which(is.finite(v))
    if (!length(ok)) next
    ord <- order(v[ok], decreasing = TRUE, method = "radix")
    keep_idx <- ok[ head(ord, k) ]
    for (j in keep_idx) out[[ syms[j] ]][r] <- TRUE
  }
  out
}

#' Map scores to portfolio weights
#' @param method `"softmax"`, `"rank"`, `"linear"`, `"equal"`.
#' @param temperature softmax temperature (higher=flatter).
#' @param floor non-negative number added before normalization.
#' @param scores Wide score panel (Date + symbols).
#' @return weight panel (`Date + symbols`) with non-negative rows summing to 1
#'   over active dates.
#' @export
weight_from_scores <- function(scores, method = c("softmax","rank","linear","equal"),
                               temperature = 10, floor = 0) {
  method <- match.arg(method)
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)
  syms <- setdiff(names(S), "Date")
  W <- data.table::data.table(Date = S$Date)
  for (nm in syms) W[[nm]] <- 0
  for (r in seq_len(nrow(S))) {
    v <- as.numeric(S[r, ..syms]); sel <- is.finite(v)
    if (!any(sel)) next
    if (method == "equal") {
      w <- rep(1, sum(sel))
    } else if (method == "softmax") {
      x <- v[sel]
      x <- x / ifelse(temperature > 0, temperature, 1)
      x <- x - max(x)
      ex <- exp(x)
      w <- ex
    } else if (method == "rank") {
      rk <- rank(v[sel], ties.method = "average")
      w <- rk
    } else { # linear
      x <- pmax(v[sel], 0)
      if (!any(x > 0)) { w <- rep(1, length(x)) } else { w <- x }
    }
    w <- w + floor
    w <- w / sum(w)
    idx <- which(sel)
    for (j in seq_along(idx)) W[[ syms[idx[j]] ]][r] <- w[j]
  }
  W
}

# Post-weight exposure caps (per symbol and/or per group)
`%||%` <- function(a,b) if (!is.null(a)) a else b

#' Apply post-weight exposure caps
#' @param caps list with `max_per_symbol` and optional `max_per_group`.
#' @param group_map optional `data.frame(Symbol, Group)`.
#' @param renorm logical; renormalize each active row to 1.
#' @param weights Wide weight panel (Date + symbols) before caps.
#' @export
cap_exposure <- function(weights,
                         caps = list(max_per_symbol = 0.10, max_per_group = 0.25),
                         group_map = NULL, renorm = TRUE) {
  W <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(W), "Date")
  max_sym <- caps$max_per_symbol %||% Inf
  # symbol cap
  for (nm in syms) W[[nm]] <- pmin(W[[nm]], max_sym)
  # group cap (optional)
  if (!is.null(group_map)) {
    stopifnot(all(c("Symbol","Group") %in% names(group_map)))
    gmap <- data.table::as.data.table(group_map)
    for (r in seq_len(nrow(W))) {
      row <- as.numeric(W[r, ..syms]); names(row) <- syms
      agg <- tapply(row, gmap$Group[match(names(row), gmap$Symbol)], sum, na.rm = TRUE)
      over <- names(agg)[which(agg > (caps$max_per_group %||% 1))]
      if (length(over)) {
        for (grp in over) {
          in_grp <- gmap$Symbol[gmap$Group == grp]
          idx <- which(syms %in% in_grp)
          if (length(idx)) {
            scale <- (caps$max_per_group %||% 1) / sum(row[idx])
            row[idx] <- row[idx] * scale
          }
        }
      }
      if (isTRUE(renorm) && sum(row, na.rm = TRUE) > 0) row <- row / sum(row, na.rm = TRUE)
      W[r, (syms) := as.list(row)]
    }
  } else if (isTRUE(renorm)) {
    for (r in seq_len(nrow(W))) {
      row <- as.numeric(W[r, ..syms])
      s <- sum(row, na.rm = TRUE)
      if (s > 0) W[r, (syms) := as.list(row / s)]
    }
  }
  W
}

# ---- Filters (optional minimal) --------------------------------------------

# Generic metric gate; returns a logical mask panel
filter_by_metric <- function(metric_df, op = c("below","above","between"),
                             threshold, threshold2 = NULL, min_qualified = 0L) {
  op <- match.arg(op)
  M <- if (data.table::is.data.table(metric_df)) data.table::copy(metric_df) else data.table::as.data.table(metric_df)
  syms <- setdiff(names(M), "Date")
  out <- data.table::data.table(Date = M$Date)
  for (nm in syms) out[[nm]] <- FALSE
  for (r in seq_len(nrow(M))) {
    v <- as.numeric(M[r, ..syms]); ok <- is.finite(v)
    if (!any(ok)) next
    keep <- switch(op,
                   below   = ok & (v < threshold),
                   above   = ok & (v > threshold),
                   between = ok & (pmin(v, threshold2) == v & pmax(v, threshold) == v)
    )
    if (sum(keep, na.rm = TRUE) < min_qualified && any(ok)) {
      ord <- order(abs(v[ok] - (if (op == "below") -Inf else if (op == "above") Inf else (threshold + threshold2)/2)), decreasing = FALSE)
      choose <- which(ok)[ head(ord, min_qualified) ]
      keep[choose] <- TRUE
    }
    out[r, (syms) := as.list(keep)]
  }
  out
}




#' Combine multiple score panels (mean / weighted / rank-average / trimmed)
#'
#' Combine several wide score panels (`Date` + symbols) into a single panel
#' by applying one of several aggregation methods.
#'
#' @details
#' \itemize{
#'   \item \code{method = "mean"}: simple column-wise mean across panels.
#'   \item \code{method = "weighted"}: weighted mean; see \code{weights}.
#'   \item \code{method = "rank_avg"}: average of within-date normalized ranks.
#'   \item \code{method = "trimmed_mean"}: mean with \code{trim} fraction removed at both tails.
#' }
#'
#' @param scores_list List of wide score panels to combine (each has columns \code{Date} + symbols).
#' @param method Character, one of \code{"mean"}, \code{"weighted"}, \code{"rank_avg"}, \code{"trimmed_mean"}.
#' @param weights Optional numeric vector of length equal to \code{length(scores_list)}; used only when \code{method = "weighted"}.
#' @param trim Numeric in \code{[0, 0.5)}; fraction to trim from each tail for \code{method = "trimmed_mean"}.
#'
#' @return A \code{data.table} with columns \code{Date} + symbols, containing the combined scores.
#' @export
combine_scores <- function(scores_list,
                           method = c("mean","weighted","rank_avg","trimmed_mean"),
                           weights = NULL, trim = 0.1) {
  method <- match.arg(method)
  stopifnot(is.list(scores_list), length(scores_list) >= 1)

  aligned <- .roll_align_common(scores_list)
  L <- aligned$panels
  dates <- aligned$dates
  syms  <- aligned$symbols

  out <- data.table::data.table(Date = dates)
  for (nm in syms) out[[nm]] <- NA_real_

  K <- length(L)
  if (method == "weighted") {
    if (is.null(weights)) stop("weights required for method='weighted'")
    if (length(weights) != K) stop("weights length must equal number of score panels")
    weights <- weights / sum(weights)
  }

  for (r in seq_along(dates)) {
    mats <- lapply(L, function(z) as.numeric(z[r, ..syms]))
    M <- do.call(rbind, mats)   # K x Nsym

    if (method == "mean") {
      v <- colMeans(M, na.rm = TRUE)

    } else if (method == "weighted") {
      v <- as.numeric(c(weights %*% M))

    } else if (method == "rank_avg") {
      R <- matrix(NA_real_, nrow = K, ncol = length(syms))
      for (k in seq_len(K)) {
        x  <- M[k, ]
        ok <- is.finite(x)
        if (!any(ok)) next
        rk <- rank(x[ok], ties.method = "average")
        z  <- rep(NA_real_, length(x))
        z[ok] <- (rk - min(rk)) / (max(rk) - min(rk))
        R[k, ] <- z
      }
      v <- colMeans(R, na.rm = TRUE)

    } else { # trimmed_mean
      v <- apply(M, 2L, function(col) {
        x <- col[is.finite(col)]
        if (!length(x)) return(NA_real_)
        base::mean(x, trim = trim, na.rm = TRUE)
      })
    }

    out[r, (syms) := as.list(v)]
  }

  out
}

#


# ==== TURNOVER CAP ==========================================================
#' Cap turnover sequentially across dates
#' @param weights desired weight panel.
#' @param max_turnover maximum per-rebalance turnover (0..1).
#' @return executed weights after turnover capping.
#' @export
cap_turnover <- function(weights, max_turnover = 0.2) {
  W <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(W), "Date")
  W[,(syms) := lapply(.SD, function(x) {x[!is.finite(x)] <- 0; x}), .SDcols = syms]

  prev <- rep(0, length(syms))
  for (r in seq_len(nrow(W))) {
    row <- as.numeric(W[r, ..syms])
    s <- sum(row, na.rm = TRUE)
    if (s > 0) row <- row / s
    tcur <- 0.5 * sum(abs(row - prev))
    if (tcur > max_turnover && tcur > 0) {
      alpha <- max_turnover / tcur
      row <- prev + alpha * (row - prev)
    }
    s <- sum(row, na.rm = TRUE)
    if (s > 0) row <- row / s
    W[r, (syms) := as.list(row)]
    prev <- row
  }
  W
}



#' Validate a symbol-to-group mapping
#'
#' Normalizes and checks a symbol → group mapping for a given set of symbols.
#' Accepts either a data.frame/data.table with columns `Symbol` and `Group`,
#' or a named character vector `c(symbol = "group", ...)`.
#' Errors if any requested symbol is missing or mapped more than once.
#'
#' @param symbols Character vector of symbols to validate/keep.
#' @param group_map Data frame/data.table with columns `Symbol`,`Group`,
#'   or a named character vector mapping `symbol -> group`.
#'
#' @return A two-column `data.frame` with columns `Symbol` and `Group`
#'   (one row per symbol), sorted by `Symbol`.
#' @export
#' @examples
#' validate_group_map(
#'   c("A","B"),
#'   data.frame(Symbol = c("A","B"), Group = c("G1","G1"))
#' )
validate_group_map <- function(symbols, group_map) {
  if (is.null(group_map)) stop("group_map is NULL")
  GM <- if (data.table::is.data.table(group_map)) data.table::copy(group_map) else data.table::as.data.table(group_map)
  if (!all(c("Symbol","Group") %in% names(GM))) stop("group_map must have columns: Symbol, Group")
  GM <- GM[Symbol %in% symbols]
  if (anyDuplicated(GM$Symbol)) stop("group_map must map each symbol to exactly one group")
  missing <- setdiff(symbols, GM$Symbol)
  if (length(missing)) stop("group_map missing symbols: ", paste(missing, collapse = ", "))
  data.table::setorder(GM, Symbol)
  GM
}



#' Demo sector (group) map for examples/tests
#'
#' @param symbols character vector of tickers.
#' @param n_groups integer number of groups to assign (cycled).
#' @return A data.table with columns `Symbol` and `Group` (title-case), for demo use.
#'
#' @examples
#' # Minimal usage
#' syms <- c("AAPL","MSFT","AMZN","XOM","JPM")
#' gdf  <- demo_sector_map(syms, n_groups = 3L)  # columns: Symbol, Group
#' print(gdf)
#'
#' # Use with cap_exposure(): convert to a named vector (names = symbols)
#' gmap <- stats::setNames(gdf$Group, gdf$Symbol)
#'
#' data(sample_prices_weekly)
#' mom12 <- calc_momentum(sample_prices_weekly, 12)
#' sel10 <- filter_top_n(mom12, 10)
#' w_eq  <- weight_equally(sel10)
#'
#' w_cap <- cap_exposure(
#'   weights          = w_eq,
#'   max_per_symbol   = 0.10,
#'   group_map        = gmap,     # <- named vector, OK for current cap_exposure()
#'   max_per_group    = 0.45,
#'   renormalize_down = TRUE
#' )
#' head(w_cap)
#'
#' @export
demo_sector_map <- function(symbols, n_groups = 2L) {
  groups <- paste0("G", seq_len(n_groups))
  data.table::data.table(Symbol = symbols, Group = groups[(seq_along(symbols) - 1L) %% n_groups + 1L])
}



scores_to_weights <- function(scores,
                              top_k = NULL,
                              group_map = NULL,
                              max_per_group = NULL,
                              transform = c("none","zscore","rank"),
                              weighting = c("softmax","equal","rank","linear"),
                              temperature = 10, floor = 0,
                              caps = list(max_per_symbol = 0.10)) {
  transform  <- match.arg(transform)
  weighting  <- match.arg(weighting)
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)

  if (transform == "zscore") S <- transform_scores(S, "zscore")
  if (transform == "rank")   S <- transform_scores(S, "rank")

  syms <- setdiff(names(S), "Date")

  if (!is.null(top_k)) {
    if (!is.null(group_map) && !is.null(max_per_group)) {
      group_map <- validate_group_map(syms, group_map)
      mask <- select_top_k_scores_by_group(S, k = top_k, group_map = group_map, max_per_group = max_per_group)
    } else {
      mask <- select_top_k_scores(S, k = top_k)
    }
    Ssel <- data.table::copy(S)
    for (nm in syms) Ssel[[nm]][!mask[[nm]]] <- NA_real_
  } else {
    Ssel <- S
  }

  W0 <- weight_from_scores(Ssel, method = weighting, temperature = temperature, floor = floor)
  cap_exposure(W0, caps = caps)
}

#' One-call backtest wrapper (tabular features)
#'
#' @inheritParams roll_fit_predict
#' @param schedule list with elements `is`, `oos`, `step`.
#' @param transform `"none"`, `"zscore"`, `"rank"` for per-date score transform.
#' @param selection list: `top_k`, `max_per_group` (optional).
#' @param weighting list: `method`, `temperature`, `floor`.
#' @param caps list: `max_per_symbol`, optionally `max_per_group`.
#' @param prices price panel used by the backtester (`Date` + symbols).
#' @param initial_capital starting capital.
#' @param name string for the backtest result.
#' @param labels Wide label panel (Date + symbols).
#' @param turnover Optional turnover cap settings (currently advisory/unused).
#' @return list: `scores`, `mask`, `weights`, `backtest`.
#' @export
#' @examples
#' \donttest{
#' data(sample_prices_weekly); data(sample_prices_daily)
#' mom <- panel_lag(calc_momentum(sample_prices_weekly, 12), 1)
#' vol <- panel_lag(align_to_timeframe(
#'   calc_rolling_volatility(sample_prices_daily, 20),
#'   sample_prices_weekly$Date, "forward_fill"), 1)
#' Y <- make_labels(sample_prices_weekly, 4, "log")
#' fit_lm  <- function(X,y){ Xc <- cbind(1,X); list(coef=stats::lm.fit(Xc,y)$coefficients) }
#' pred_lm <- function(m,X){ as.numeric(cbind(1,X) %*% m$coef) }
#' res <- ml_backtest(list(mom=mom, vol=vol), Y, fit_lm, pred_lm,
#'                    schedule = list(is=52,oos=4,step=4),
#'                    transform = "zscore",
#'                    selection = list(top_k=10),
#'                    weighting = list(method="softmax", temperature=12),
#'                    caps = list(max_per_symbol=0.10),
#'                    prices = sample_prices_weekly, initial_capital = 1e5)
#' res$backtest
#' }
ml_backtest <- function(features_list,
                        labels,
                        fit_fn, predict_fn,
                        schedule        = list(is = 156L, oos = 4L, step = 4L),
                        group           = c("pooled","per_symbol","per_group"),
                        transform       = c("none","zscore","rank"),
                        selection       = list(top_k = 15L, max_per_group = NULL),
                        group_map       = NULL,
                        weighting       = list(method = "softmax", temperature = 12, floor = 0),
                        caps            = list(max_per_symbol = 0.10, max_per_group = NULL),
                        turnover        = NULL,
                        prices,
                        initial_capital = 1e5,
                        name            = "ML backtest") {

  `%||%` <- function(a,b) if (!is.null(a)) a else b
  group     <- match.arg(group)
  transform <- match.arg(transform)

  as_dt_sorted <- function(x) {
    dx <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!"Date" %in% names(dx)) stop("Panel missing 'Date' column")
    if (!inherits(dx$Date, "Date")) dx[, Date := as.Date(Date)]
    data.table::setorderv(dx, "Date")
    dx
  }
  F_list <- lapply(features_list, as_dt_sorted)
  Lbl    <- as_dt_sorted(labels)
  Px     <- as_dt_sorted(prices)

  # Align across features + labels
  aligned <- .roll_align_common(c(F_list, list(Lbl)))
  L       <- aligned$panels
  dates   <- aligned$dates
  syms    <- aligned$symbols
  lab     <- L[[length(L)]]
  feats   <- L[seq_len(length(L) - 1L)]

  # Rolling fit/predict ONCE with all features
  S <- roll_fit_predict(
    features_list = feats,
    label         = lab,
    fit_fn        = fit_fn,
    predict_fn    = predict_fn,
    is_periods    = schedule$is,
    oos_periods   = schedule$oos,
    step          = schedule$step,
    group         = group,
    group_map     = group_map
  )

  # Optional transform
  if (transform != "none") S <- transform_scores(S, method = transform)

  # Selection (with optional per-group cap)
  mask <- if (!is.null(selection$max_per_group) && !is.null(group_map)) {
    select_top_k_scores_by_group(S,
                                 k             = selection$top_k %||% 15L,
                                 group_map     = group_map,
                                 max_per_group = selection$max_per_group)
  } else {
    select_top_k_scores(S, k = selection$top_k %||% 15L)
  }

  # Mask -> weights -> caps
  Sc_sel <- data.table::copy(S)
  for (nm in syms) Sc_sel[[nm]][!mask[[nm]]] <- NA_real_

  W0 <- weight_from_scores(
    Sc_sel,
    method      = weighting$method %||% "softmax",
    temperature = weighting$temperature %||% 12,
    floor       = weighting$floor %||% 0
  )
  W <- cap_exposure(
    W0,
    caps      = list(max_per_symbol = caps$max_per_symbol %||% 0.10,
                     max_per_group  = caps$max_per_group  %||% NULL),
    group_map = group_map,
    renorm    = TRUE
  )

  # Align prices and backtest
  Px2 <- Px[Px[["Date"]] %in% dates, c("Date", syms), with = FALSE]
  BT  <- run_backtest(prices = Px2, weights = W, initial_capital = initial_capital, name = name)

  list(scores = S, mask = mask, weights = W, backtest = BT)
}

#' Select top-k symbols per group by score
#'
#' @description
#' For each date, choose the top `k` symbols **within each group** based on a
#' score panel. Returns a logical selection panel aligned to the input.
#'
#' @details
#' \itemize{
#'   \item Group membership comes from `group_map` (symbol -> group).
#'   \item Selection is computed independently by group on each date.
#'   \item Ties follow the ordering implied by `order(..., method = "radix")`.
#' }
#'
#' @param scores Wide score panel (`Date` + symbols).
#' @param k Positive integer: number of symbols to select per group.
#' @param group_map Named character vector or 2-column data.frame
#'   (`symbol`, `group`) mapping symbols to groups.
#' @param max_per_group Integer cap per group (default `3L`).
#'
#' @return Logical selection panel (`Date` + symbols) where `TRUE` marks
#'         selected symbols.
#' @examples
#' set.seed(42)
#' scores <- data.frame(
#'   Date = as.Date("2020-01-01") + 0:1,
#'   A = runif(2), B = runif(2), C = runif(2), D = runif(2), E = runif(2), F = runif(2)
#' )
#' gmap <- data.frame(Symbol = c("A","B","C","D","E","F"),
#'                    Group  = c("G1","G1","G2","G2","G3","G3"))
#' sel <- select_top_k_scores_by_group(scores, k = 4, group_map = gmap, max_per_group = 2)
#' head(sel)
#' @export
select_top_k_scores_by_group <- function(scores, k, group_map, max_per_group = 3L) {
  S <- if (data.table::is.data.table(scores)) data.table::copy(scores) else data.table::as.data.table(scores)
  syms <- setdiff(names(S), "Date")

  # Normalize map: returns data.frame with columns Symbol, Group
  GM <- validate_group_map(syms, group_map)
  sym2grp <- setNames(GM$Group, GM$Symbol)
  groups  <- unique(GM$Group)

  out <- data.table::data.table(Date = S$Date)
  for (nm in syms) out[[nm]] <- FALSE

  for (r in seq_len(nrow(S))) {
    v  <- as.numeric(S[r, ..syms])
    ok <- which(is.finite(v))
    if (!length(ok)) next

    ord <- ok[order(v[ok], decreasing = TRUE, method = "radix")]
    picked <- character(0)
    grp_ct <- setNames(integer(length(groups)), groups)

    for (j in ord) {
      s <- syms[j]; g <- sym2grp[[s]]
      if (is.null(g) || is.na(g)) next
      if (grp_ct[g] < max_per_group) {
        picked <- c(picked, s)
        grp_ct[g] <- grp_ct[g] + 1L
        if (length(picked) >= k) break
      }
    }
    if (length(picked)) for (s in picked) out[[s]][r] <- TRUE
  }
  out
}



# ==== Score diagnostics (robust alignment) ===================================
#' Evaluate scores vs labels (IC and hit-rate)
#' @param top_frac fraction in the top bucket for hit-rate.
#' @param method `"spearman"` or `"pearson"`.
#' @param scores Wide score panel.
#' @param labels Wide label panel aligned to scores.
#' @return data.table with `Date, IC, hit_rate`; `ICIR` on `attr(result,"ICIR")`.
#' @export
evaluate_scores <- function(scores, labels, top_frac = 0.2, method = c("spearman","pearson")) {
  method <- match.arg(method)
  stopifnot(is.numeric(top_frac), top_frac > 0, top_frac <= 1)

  # Align by common dates & symbols
  aligned <- .roll_align_common(list(scores, labels))
  S <- aligned$panels[[1L]]
  L <- aligned$panels[[2L]]
  dates <- aligned$dates
  syms  <- aligned$symbols

  res <- data.table::data.table(Date = dates, IC = NA_real_, hit_rate = NA_real_)

  for (r in seq_along(dates)) {
    x <- as.numeric(S[r, ..syms])
    y <- as.numeric(L[r, ..syms])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3L) {
      # Information Coefficient
      ic <- suppressWarnings(stats::cor(x[ok], y[ok], method = method))
      res$IC[r] <- if (is.finite(ic)) ic else NA_real_

      # Hit-rate of top bucket by scores
      k <- max(1L, floor(sum(ok) * top_frac))
      ord <- order(x[ok], decreasing = TRUE)
      top_idx <- which(ok)[ head(ord, k) ]
      hr <- mean(y[top_idx] > 0, na.rm = TRUE)
      res$hit_rate[r] <- if (is.finite(hr)) hr else NA_real_
    }
  }

  ic_vals <- res$IC[is.finite(res$IC)]
  ICIR <- if (length(ic_vals) >= 2L && stats::sd(ic_vals) > 0)
    base::mean(ic_vals) / stats::sd(ic_vals) else NA_real_

  attr(res, "ICIR") <- ICIR
  res
}



# ==== Research diagnostics ===================================================
#' Count finite entries per date
#' @param panel wide panel `Date + symbols`.
#' @return data.table with `Date, n_finite`.
#' @export
coverage_by_date <- function(panel) {
  DT <- if (data.table::is.data.table(panel)) data.table::copy(panel) else data.table::as.data.table(panel)
  syms <- setdiff(names(DT), "Date")
  cov <- apply(as.matrix(DT[, ..syms]), 1L, function(x) sum(is.finite(x)))
  data.table::data.table(Date = DT$Date, n_finite = cov)
}



#' Information Coefficient time series
#'
#' @param scores Wide score panel (Date + symbols).
#' @param labels Wide label panel aligned to scores.
#' @param method IC method ('spearman' or 'pearson').
#'
#' @return data.table with `Date` and `IC`.
#' @export
ic_series <- function(scores, labels, method = c("spearman","pearson")) {
  method <- match.arg(method)
  aligned <- .roll_align_common(list(scores, labels))
  S <- aligned$panels[[1L]]; L <- aligned$panels[[2L]]
  syms <- aligned$symbols
  IC <- rep(NA_real_, nrow(S))
  for (r in seq_len(nrow(S))) {
    x <- as.numeric(S[r, ..syms]); y <- as.numeric(L[r, ..syms])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3L) IC[r] <- suppressWarnings(stats::cor(x[ok], y[ok], method = method))
  }
  data.table::data.table(Date = S$Date, IC = IC)
}


#' Bucketed label analysis by score rank
#' @param n_buckets number of equal-count buckets.
#' @param label_type `"log"` or `"simple"` for interpretation of labels.
#' @param scores Wide score panel (Date + symbols).
#' @param labels Wide label panel aligned to scores (Date + symbols).
#' @param n_buckets Number of equal-count buckets.
#' @param label_type Use 'log' or 'simple' label arithmetic.
#' @return data.table of per-date bucket returns; cumulative series attached
#'   as `attr(result, "cum")`.
#' @export
bucket_returns <- function(scores, labels, n_buckets = 10L, label_type = c("log","simple")) {
  label_type <- match.arg(label_type)
  aligned <- .roll_align_common(list(scores, labels))
  S <- aligned$panels[[1L]]; L <- aligned$panels[[2L]]
  syms <- aligned$symbols

  out <- data.table::data.table(Date = S$Date)
  for (b in seq_len(n_buckets)) out[[paste0("B", b)]] <- NA_real_

  for (r in seq_len(nrow(S))) {
    sc <- as.numeric(S[r, ..syms])
    lb <- as.numeric(L[r, ..syms])
    ok <- is.finite(sc) & is.finite(lb)
    if (sum(ok) < n_buckets) next
    # equal-count buckets by rank
    rk <- rank(sc[ok], ties.method = "average")
    qs <- cut(rk, breaks = quantile(rk, probs = seq(0,1,length.out = n_buckets+1), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
    for (b in seq_len(n_buckets)) {
      idx <- which(ok)[qs == b]
      if (!length(idx)) next
      out[[paste0("B", b)]][r] <- mean(lb[idx], na.rm = TRUE)
    }
  }

  # attach cumulative series (log->exp(cumsum), simple->cumprod(1+r))
  cum <- data.table::data.table(Date = out$Date)
  for (b in seq_len(n_buckets)) {
    col <- paste0("B", b)
    rts <- out[[col]]
    if (label_type == "log") {
      cum[[col]] <- exp(cumsum(replace(rts, !is.finite(rts), 0)))
    } else {
      rts[!is.finite(rts)] <- 0
      cum[[col]] <- cumprod(1 + rts)
    }
  }
  attr(out, "cum") <- cum
  out
}

#' Rolling IC mean, standard deviation, and ICIR
#'
#' Compute rolling information coefficient (IC) statistics from a per-date IC series.
#'
#' @details
#' For each rolling window, compute the mean IC, the standard deviation of IC,
#' and the information coefficient ratio (ICIR = mean / sd). Windows with fewer
#' than two finite IC values yield \code{NA} for ICIR.
#'
#' @param ic_dt Data frame/data.table produced by \code{\link{ic_series}} with columns \code{Date} and \code{IC}.
#' @param window Integer window length for the rolling statistics.
#'
#' @return A \code{data.table} with columns \code{Date}, \code{IC_mean}, \code{IC_sd}, and \code{ICIR}.
#' @export
roll_ic_stats <- function(ic_dt, window = 26L) {
  stopifnot(all(c("Date","IC") %in% names(ic_dt)))
  x <- ic_dt$IC
  rm <- zoo::rollapply(x, width = window, FUN = function(z) mean(z, na.rm = TRUE), fill = NA, align = "right")
  rs <- zoo::rollapply(x, width = window, FUN = function(z) stats::sd(z, na.rm = TRUE),   fill = NA, align = "right")
  data.table::data.table(Date = ic_dt$Date, roll_mean = rm, roll_sd = rs, roll_ICIR = rm / rs)
}



# ==== Group-aware selection ===================================================



#' Rebalance calendar (rows with non-zero allocation)
#' @return data.table with `Date, row` indices where a rebalance occurred.
#' @param weights Wide weight panel (Date + symbols).
#' @export
rebalance_calendar <- function(weights) {
  DT <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(DT), "Date")
  idx  <- which(rowSums(as.matrix(DT[, ..syms])) > 0)
  data.table::data.table(Date = DT$Date[idx], row = idx)
}


#' Carry-forward weights between rebalances (validation helper)
#' @return weight panel with rows filled forward and each active row sum=1.
#' @param weights Wide weight panel (Date + symbols) with weights only on rebalance rows.
#' @export
carry_forward_weights <- function(weights) {
  DT <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(DT), "Date")
  last <- rep(0, length(syms)); names(last) <- syms
  for (r in seq_len(nrow(DT))) {
    w <- as.numeric(DT[r, ..syms])
    if (sum(w, na.rm = TRUE) > 0) {
      last <- w
    } else if (sum(last) > 0) {
      DT[r, (syms) := as.list(last)]
    }
  }
  # normalize each active row to 1
  for (r in seq_len(nrow(DT))) {
    w <- as.numeric(DT[r, ..syms]); s <- sum(w, na.rm = TRUE)
    if (s > 0) DT[r, (syms) := as.list(w / s)]
  }
  DT
}




# ==== RETURNS & PERFORMANCE ===================================================
#' Panel simple returns from prices
#'
#' Converts a wide price panel (Date + symbols) into arithmetic simple returns
#' at the same cadence, dropping the first row per symbol.
#'
#' @param prices A data frame or data.table with columns \code{Date} and one column
#'   per symbol containing adjusted prices at a common frequency (daily, weekly, monthly).
#'
#' @return A data frame with \code{Date} and one column per symbol containing simple returns
#'   \eqn{R_{t} = P_{t}/P_{t-1} - 1}.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   rets <- panel_returns_simple(sample_prices_weekly)
#'   head(rets)
#' }
#'
#' @export
panel_returns_simple <- function(prices) {
  DT <- if (data.table::is.data.table(prices)) data.table::copy(prices) else data.table::as.data.table(prices)
  stopifnot("Date" %in% names(DT))
  syms <- setdiff(names(DT), "Date")
  R <- data.table::data.table(Date = DT$Date)
  for (nm in syms) {
    p <- as.numeric(DT[[nm]])
    R[[nm]] <- c(NA_real_, p[-1L] / head(p, -1L) - 1)
  }
  R
}
#' Portfolio returns from weights and prices (CASH-aware)
#'
#' Computes the portfolio simple return series by applying (lagged) portfolio
#' weights to next-period asset returns, optionally net of proportional costs.
#'
#' **CASH support:** if `weights` contains a column named `"CASH"` (case-insensitive)
#' but `prices` has no matching column, a synthetic flat price series is added
#' internally (price = 1 \eqn{\Rightarrow}{=>} return = 0). In that case the function does **not**
#' re-normalise the non-CASH weights; the row is treated as a complete budget
#' (symbols + CASH = 1).
#'
#' @param weights A data.frame/data.table of portfolio weights on rebalance dates:
#'   first column `Date`, remaining columns one per symbol (numeric weights).
#'   Weights decided at \eqn{t-1} are applied to returns over \eqn{t}.
#' @param prices A data.frame/data.table of adjusted prices at the same cadence:
#'   first column `Date`, remaining columns one per symbol.
#' @param cost_bps One-way proportional cost per side in basis points (e.g., `10`
#'   for 10 bps). Default `0`. If `> 0` and your package exposes a turnover
#'   helper, it will be used; otherwise costs are ignored with a warning.
#'
#' @return A `data.table` with columns `Date` and `ret` (portfolio simple return).
#'
#' @details
#' The function carries forward the latest available weights to each return row
#' via the usual one-period decision lag. Transaction cost handling is conservative:
#' if a turnover helper is not available, costs are skipped.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)
#'   sel10 <- PortfolioTesteR::filter_top_n(mom12, 10)
#'   w_eq  <- PortfolioTesteR::weight_equally(sel10)
#'
#'   pr <- portfolio_returns(w_eq, sample_prices_weekly, cost_bps = 0)
#'   head(pr)
#' }
#'
#' @seealso PortfolioTesteR::panel_returns_simple
#' @export
portfolio_returns <- function(weights, prices, cost_bps = 0) {
  stopifnot(is.data.frame(weights), is.data.frame(prices))
  stopifnot("Date" %in% names(weights), "Date" %in% names(prices))

  W <- data.table::as.data.table(data.table::copy(weights))
  P <- data.table::as.data.table(data.table::copy(prices))
  data.table::setkey(W, Date); data.table::setkey(P, Date)

  w_syms <- setdiff(names(W), "Date")
  p_syms <- setdiff(names(P), "Date")

  # Detect CASH (case-insensitive) in weights
  cash_col <- w_syms[tolower(w_syms) == "cash"]
  if (length(cash_col) > 1L) {
    stop("portfolio_returns(): multiple CASH-like columns in weights.")
  }
  has_cash_in_w <- length(cash_col) == 1L

  # If CASH present in weights but missing in prices, synthesize a flat price
  if (has_cash_in_w && !(cash_col %in% p_syms)) {
    P[, (cash_col) := 1]  # flat price \eqn{\Rightarrow}{=>} 0 simple return
    p_syms <- c(p_syms, cash_col)
  }

  # Intersect symbols actually available in both
  syms <- intersect(w_syms, p_syms)
  if (!length(syms)) stop("portfolio_returns(): no overlapping symbols between weights and prices.")

  # Asset simple returns (Date + syms)
  R <- PortfolioTesteR::panel_returns_simple(P[, c("Date", syms), with = FALSE])

  # Decision lag: weights at t-1 apply to returns at t
  R2   <- R[-1]
  Wlag <- W[-.N]
  stopifnot(nrow(R2) == nrow(Wlag))

  # Compute portfolio returns row-wise
  Xw <- as.matrix(Wlag[, ..syms])
  Xr <- as.matrix(R2[, ..syms])
  pr <- rowSums(Xw * Xr, na.rm = TRUE)

  # Optional transaction costs (best-effort; ignored if helper missing)
  if (is.numeric(cost_bps) && isTRUE(cost_bps > 0)) {
    add_costs <- FALSE
    # Try to use a turnover helper if available (best-effort)
    if (exists("turnover_by_date", where = asNamespace("PortfolioTesteR"), inherits = FALSE)) {
      # Attempt a safe call; if it fails, skip costs with a warning
      to_dt <- try(PortfolioTesteR::turnover_by_date(W, P), silent = TRUE)
      if (!inherits(to_dt, "try-error") && is.data.frame(to_dt) && "turnover" %in% names(to_dt)) {
        # One-way cost on gross turnover (common convention)
        c_perc <- as.numeric(cost_bps) / 1e4
        cost   <- c_perc * to_dt$turnover
        cost   <- cost[seq_along(pr)]  # align length defensively
        cost[!is.finite(cost)] <- 0
        pr <- pr - cost
        add_costs <- TRUE
      }
    }
    if (!add_costs) {
      warning("portfolio_returns(): cost_bps > 0 but no compatible turnover helper found; costs ignored.")
    }
  }

  data.table::data.table(Date = R2$Date, ret = pr)[]
}

#' Portfolio performance metrics
#'
#' @param ret_dt data.table/data.frame with columns `Date` and `ret`.
#' @param freq Integer periods-per-year for annualization (e.g., 52 or 252).
#' @return list with `total_return`, `ann_return`, `ann_vol`, `sharpe`, `max_drawdown`.
#' @export
perf_metrics <- function(ret_dt, freq = 52) {
  stopifnot(all(c("Date","ret") %in% names(ret_dt)))
  r <- ret_dt$ret[is.finite(ret_dt$ret)]
  if (!length(r)) {
    return(list(total_return = NA_real_, ann_return = NA_real_,
                ann_vol = NA_real_, sharpe = NA_real_, max_drawdown = NA_real_))
  }
  n  <- length(r)
  TR <- prod(1 + r) - 1
  ann_ret <- (1 + TR)^(freq / n) - 1
  ann_vol <- stats::sd(r) * sqrt(freq)
  sharpe  <- if (is.finite(ann_vol) && ann_vol > 0) ann_ret / ann_vol else NA_real_

  wealth <- cumprod(replace(1 + ret_dt$ret, !is.finite(ret_dt$ret), 0))
  dd     <- wealth / cummax(wealth) - 1
  max_dd <- min(dd, na.rm = TRUE)

  list(total_return = TR, ann_return = ann_ret, ann_vol = ann_vol,
       sharpe = sharpe, max_drawdown = max_dd)
}


# ==== SIMPLE PARAMETER GRID TUNER ============================================

#' Quick grid tuning for tabular pipeline
#' @param grid list of vectors: `top_k`, `temperature`, `method`, `transform`.
#' @param cost_bps optional one-way cost in basis points for net performance.
#' @param freq re-annualization frequency (e.g., 52).
#' @param features_list List of feature panels.
#' @param labels Label panel.
#' @param prices Price panel used for backtests (Date + symbols).
#' @param fit_fn,predict_fn Model fit and predict functions.
#' @param schedule List with elements is, oos, step.
#' @param group Grouping mode for roll_fit_predict ('pooled'/'per_symbol'/'per_group').
#' @param selection_defaults Default selection settings (e.g., top_k).
#' @param weighting_defaults Default weighting settings (e.g., method, temperature).
#' @param caps Exposure caps (e.g., max_per_symbol/max_per_group).
#' @param group_map Optional Symbol->Group mapping.
#' @return data.table with metrics per grid row.
#' @export
tune_ml_backtest <- function(features_list, labels, prices,
                             fit_fn, predict_fn,
                             schedule = list(is = 104L, oos = 4L, step = 4L),
                             grid = list(
                               top_k = c(10L, 15L),
                               temperature = c(8, 12),
                               method = c("softmax","rank"),
                               transform = c("zscore")
                             ),
                             group = "pooled",
                             selection_defaults = list(top_k = 15L, max_per_group = NULL),
                             weighting_defaults = list(method = "softmax", temperature = 12, floor = 0),
                             caps = list(max_per_symbol = 0.08),
                             group_map = NULL,
                             cost_bps = 0,
                             freq = 52) {

  # build all parameter combinations
  combos <- do.call(base::expand.grid, c(grid, stringsAsFactors = FALSE))
  rows <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    params <- as.list(combos[i, ])

    sel <- selection_defaults
    if (!is.null(params$top_k))          sel$top_k <- as.integer(params$top_k)
    if (!is.null(params$max_per_group))  sel$max_per_group <- as.integer(params$max_per_group)

    wt <- weighting_defaults
    if (!is.null(params$method))       wt$method <- params$method
    if (!is.null(params$temperature))  wt$temperature <- as.numeric(params$temperature)
    tr <- if (!is.null(params$transform)) params$transform else "none"

    res <- ml_backtest(
      features_list   = features_list,
      labels          = labels,
      fit_fn          = fit_fn,
      predict_fn      = predict_fn,
      schedule        = schedule,
      group           = group,
      transform       = tr,
      selection       = sel,
      group_map       = group_map,
      weighting       = wt,
      caps            = caps,
      prices          = prices,
      initial_capital = 1e5,
      name            = "grid"
    )

    # cost-aware metrics (independent of run_backtest internals)
    pr  <- portfolio_returns(res$weights, prices, cost_bps = cost_bps)
    met <- perf_metrics(pr, freq = freq)

    rows[[i]] <- c(params, list(
      total_return = met$total_return,
      ann_return   = met$ann_return,
      ann_vol      = met$ann_vol,
      sharpe       = met$sharpe,
      max_dd       = met$max_drawdown
    ))
  }

  data.table::rbindlist(lapply(rows, data.table::as.data.table), fill = TRUE)[order(-sharpe)]
}



# ==== SEQUENCE BUILDING & RAM GUARDS =========================================

.seq_mem_bytes <- function(n_samples, steps, p) {
  as.double(n_samples) * as.double(steps) * as.double(p) * 8.0  # doubles
}
.human_bytes <- function(b) {
  units <- c("B","KB","MB","GB","TB"); i <- 1L
  while (b >= 1024 && i < length(units)) { b <- b/1024; i <- i+1L }
  sprintf("%.2f %s", b, units[i])
}

# Align a list of feature panels + a label panel on common symbols/dates
.align_features_labels <- function(features_list, labels) {
  as_dt_sorted <- function(x) {
    dx <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!"Date" %in% names(dx)) stop("Panel missing 'Date' column")
    if (!inherits(dx$Date, "Date")) dx[, Date := as.Date(Date)]
    data.table::setorderv(dx, "Date")
    dx
  }
  F_list <- lapply(features_list, as_dt_sorted)
  Lbl    <- as_dt_sorted(labels)

  feat_syms <- lapply(F_list, function(d) setdiff(names(d), "Date"))
  lab_syms  <- setdiff(names(Lbl), "Date")
  syms      <- Reduce(intersect, c(feat_syms, list(lab_syms)))
  if (!length(syms)) stop("No common symbol columns across features and labels.")

  dates <- Reduce(intersect, c(lapply(F_list, `[[`, "Date"), list(Lbl$Date)))
  if (!length(dates)) stop("No common dates across features and labels.")

  F_list <- lapply(F_list, function(d) d[d[["Date"]] %in% dates, c("Date", syms), with = FALSE])
  Lbl    <- Lbl[Lbl[["Date"]] %in% dates, c("Date", syms), with = FALSE]
  list(features = F_list, labels = Lbl, dates = dates, symbols = syms)
}

# Build (samples x steps x p) sequences + y safely, with RAM guard
# NOTE: for research hygiene, DO NOT normalize here (to avoid leakage).
as_sequences <- function(features_list, labels, steps = 26L, horizon = 4L,
                         symbols = NULL,
                         min_dates_per_symbol = NULL,
                         max_samples = NULL,
                         sample_mode = c("tail","head","random"),
                         mem_budget_gb = getOption("PortfolioTesteR.ml.max_mem_gb", 1),
                         return = c("array")) {
  sample_mode <- match.arg(sample_mode)
  return      <- match.arg(return)

  ALN <- .align_features_labels(features_list, labels)
  F_list <- ALN$features; L <- ALN$labels
  dates  <- ALN$dates;    syms <- if (is.null(symbols)) ALN$symbols else intersect(ALN$symbols, symbols)
  p <- length(F_list); n <- length(dates)

  if (is.null(min_dates_per_symbol)) min_dates_per_symbol <- steps + horizon
  keep_syms <- syms[vapply(syms, function(s) {
    all(vapply(F_list, function(FX) sum(is.finite(FX[[s]])) >= min_dates_per_symbol, FUN.VALUE = TRUE))
  }, FUN.VALUE = TRUE)]
  if (!length(keep_syms)) stop("No symbols meet 'min_dates_per_symbol' requirement.")
  syms <- keep_syms

  # possible decision indices t (sequence ends at t; label at t)
  t_start <- steps
  t_end   <- n - horizon
  if (t_end < t_start) stop("Not enough dates for the requested steps/horizon.")
  t_idx   <- t_start:t_end

  # pre-count samples (pooled across symbols)
  n_samples <- length(t_idx) * length(syms)
  est_bytes <- .seq_mem_bytes(n_samples, steps, p)
  if (est_bytes > (mem_budget_gb * 1024^3) && is.null(max_samples)) {
    stop(sprintf(
      "Sequence build would allocate ~%s (samples~%d, steps=%d, p=%d) > budget=%.2f GB.\nSet max_samples=... or reduce steps/features.",
      .human_bytes(est_bytes), n_samples, steps, p, mem_budget_gb
    ))
  }
  if (!is.null(max_samples) && max_samples < n_samples) {
    # subselect (t, symbol) pairs
    pairs <- expand.grid(t = t_idx, s = syms, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    pick_ix <- switch(sample_mode,
                      tail   = tail(seq_len(nrow(pairs)), max_samples),
                      head   = seq_len(max_samples),
                      random = sample.int(nrow(pairs), max_samples)
    )
    pairs <- pairs[pick_ix, , drop = FALSE]
  } else {
    pairs <- expand.grid(t = t_idx, s = syms, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }

  S <- nrow(pairs)
  X <- array(NA_real_, dim = c(S, steps, p))
  y <- rep(NA_real_, S)

  # fill
  for (i in seq_len(S)) {
    t <- pairs$t[i]
    s <- pairs$s[i]
    # build [steps x p]
    for (j in seq_len(p)) {
      v <- F_list[[j]][[s]]
      seg <- v[(t - steps + 1):t]
      X[i, , j] <- as.numeric(seg)
    }
    y[i] <- as.numeric(L[[s]][t])
  }

  list(
    X = X,                 # [S x steps x p]
    y = y,                 # [S]
    pairs = pairs,         # data.frame with (t, s)
    dates = dates,
    symbols = syms,
    steps = steps,
    p = p,
    horizon = horizon,
    meta = list(
      est_bytes = est_bytes,
      est_human = .human_bytes(est_bytes),
      built_bytes = .seq_mem_bytes(S, steps, p),
      built_human = .human_bytes(.seq_mem_bytes(S, steps, p))
    )
  )
}

# Flatten [S x steps x p] to [S x (steps*p)] (for linear / tree models)
.flatten_seq3d_to_2d <- function(X3) {
  d <- dim(X3); stopifnot(length(d) == 3L)
  matrix(aperm(X3, c(1,2,3)), nrow = d[1], ncol = d[2]*d[3])
}

# Simple IS-only scaler (mean/sd per column) to avoid leakage
.compute_scaler_2d <- function(X2) {
  m <- colMeans(X2, na.rm = TRUE)
  s <- apply(X2, 2L, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  list(mean = m, sd = s)
}
.apply_scaler_2d <- function(X2, scaler) {
  sweep(sweep(X2, 2L, scaler$mean, "-"), 2L, scaler$sd, "/")
}

#' Rolling fit/predict for sequence models (flattened steps-by-p features)
#'
#' @param features_list list of panels to be stacked over `steps` history.
#' @param labels future-return panel aligned to the features.
#' @param steps int; lookback length (e.g., 26).
#' @param horizon int; label horizon (e.g., 4).
#' @param normalize `"none"`, `"zscore"`, or `"minmax"` applied using IS data only.
#' @param min_train_samples Optional minimum IS samples required to fit; if not met, skip fit.
#' @inheritParams roll_fit_predict
#' @return wide panel of scores.
#' @export
#' @examples
#' \donttest{
#' data(sample_prices_weekly); data(sample_prices_daily)
#' mom <- panel_lag(calc_momentum(sample_prices_weekly, 12), 1)
#' vol <- panel_lag(align_to_timeframe(
#'   calc_rolling_volatility(sample_prices_daily, 20),
#'   sample_prices_weekly$Date, "forward_fill"), 1)
#' Y <- make_labels(sample_prices_weekly, horizon = 4, type = "log")
#' fit_lm  <- function(X,y){ Xc <- cbind(1,X); list(coef=stats::lm.fit(Xc,y)$coefficients) }
#' pred_lm <- function(m,X){ as.numeric(cbind(1,X) %*% m$coef) }
#' S <- roll_fit_predict_seq(list(mom=mom, vol=vol), Y,
#'                           steps = 26, horizon = 4,
#'                           fit_fn = fit_lm, predict_fn = pred_lm,
#'                           is_periods = 104, oos_periods = 4, step = 4)
#' head(S)
#' }
roll_fit_predict_seq <- function(features_list, labels,
                                 steps = 26L, horizon = 4L,
                                 fit_fn, predict_fn,
                                 is_periods = 156L, oos_periods = 4L, step = 4L,
                                 group = c("pooled","per_symbol","per_group"),
                                 group_map = NULL,
                                 normalize = c("none","zscore","minmax"),
                                 min_train_samples = 50L,
                                 na_action = c("omit","zero")) {

  group     <- match.arg(group)
  normalize <- match.arg(normalize)
  na_action <- match.arg(na_action)

  # --- Coerce/sort & align common symbols/dates
  as_dt_sorted <- function(x) {
    dx <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!"Date" %in% names(dx)) stop("Panel missing 'Date' column")
    if (!inherits(dx$Date, "Date")) dx[, Date := as.Date(Date)]
    data.table::setorderv(dx, "Date")
    dx
  }
  F_list <- lapply(features_list, as_dt_sorted)
  Lbl    <- as_dt_sorted(labels)

  feat_syms    <- lapply(F_list, function(d) setdiff(names(d), "Date"))
  lab_syms     <- setdiff(names(Lbl), "Date")
  syms         <- Reduce(intersect, c(feat_syms, list(lab_syms)))
  if (!length(syms)) stop("No common symbol columns across features and labels.")
  dates_common <- Reduce(intersect, c(lapply(F_list, `[[`, "Date"), list(Lbl$Date)))
  if (!length(dates_common)) stop("No common dates across features and labels.")

  F_list <- lapply(F_list, function(d) d[d[["Date"]] %in% dates_common, c("Date", syms), with = FALSE])
  Lbl    <- Lbl[Lbl[["Date"]] %in% dates_common, c("Date", syms), with = FALSE]

  dates <- Lbl$Date
  p     <- length(F_list)               # number of feature panels
  n     <- length(dates)

  # memory heads-up
  .seq_mem_warn(n_samples = max(0L, (is_periods - steps + 1L)) * length(syms),
                steps = steps, p = p, max_bytes = 1e8)

  # helpers ----------------------------------------------------------
  seq_vec_at <- function(t_idx, sym) {
    if (t_idx - steps + 1L < 1L) return(NULL)
    out <- numeric(steps * p)
    pos <- 1L
    for (j in seq_len(p)) {
      x <- as.numeric(F_list[[j]][[sym]])
      seg <- x[(t_idx - steps + 1L):t_idx]
      out[pos:(pos + steps - 1L)] <- seg
      pos <- pos + steps
    }
    if (any(!is.finite(out))) {
      if (na_action == "omit") return(NULL)
      if (na_action == "zero") out[!is.finite(out)] <- 0
    }
    out
  }

  normalize_fit <- function(Xis) {
    if (normalize == "none") return(list(transform = function(M) M))
    if (normalize == "zscore") {
      mu <- colMeans(Xis, na.rm = TRUE)
      sdv <- apply(Xis, 2L, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
      list(transform = function(M) sweep(sweep(M, 2L, mu, "-"), 2L, sdv, "/"))
    } else { # minmax
      mn <- apply(Xis, 2L, min, na.rm = TRUE)
      mx <- apply(Xis, 2L, max, na.rm = TRUE)
      rng <- mx - mn; rng[!is.finite(rng) | rng == 0] <- 1
      list(transform = function(M) sweep(sweep(M, 2L, mn, "-"), 2L, rng, "/"))
    }
  }

  if (group == "per_group") {
    if (is.null(group_map)) stop("group_map is required when group='per_group'")
    GM <- data.table::as.data.table(group_map)
    if (!all(c("Symbol","Group") %in% names(GM))) stop("group_map must have Symbol, Group")
    GM <- GM[Symbol %in% syms]
    sym2grp <- setNames(GM$Group, GM$Symbol)
    groups  <- unique(GM$Group)
  }

  # output scores
  S <- data.table::data.table(Date = dates)
  for (nm in syms) S[[nm]] <- NA_real_

  i <- max(steps, 1L)
  while (TRUE) {
    is_start <- i
    is_end   <- i + is_periods - 1L
    oos_end  <- is_end + oos_periods
    if (oos_end > n) break

    build_is_design <- function(sym_subset = syms) {
      Xis <- NULL; yis <- NULL
      for (t in is_start:is_end) for (nm in sym_subset) {
        xv <- seq_vec_at(t, nm); yv <- as.numeric(Lbl[[nm]][t])
        if (!is.null(xv) && is.finite(yv)) {
          Xis <- if (is.null(Xis)) matrix(xv, nrow = 1L) else rbind(Xis, xv)
          yis <- c(yis, yv)
        }
      }
      list(Xis = Xis, yis = yis)
    }

    if (group == "pooled") {
      D <- build_is_design(syms)
      if (!is.null(D$Xis) && length(D$yis) >= min_train_samples) {
        scal  <- normalize_fit(D$Xis)
        model <- list(obj = fit_fn(scal$transform(D$Xis), D$yis), scal = scal)
      } else model <- NULL

    } else if (group == "per_symbol") {
      model <- vector("list", length(syms)); names(model) <- syms
      for (nm in syms) {
        D <- build_is_design(sym_subset = nm)
        if (!is.null(D$Xis) && length(D$yis) >= min_train_samples) {
          scal <- normalize_fit(D$Xis)
          model[[nm]] <- list(obj = fit_fn(scal$transform(D$Xis), D$yis), scal = scal)
        }
      }

    } else { # per_group
      model <- vector("list", length(groups)); names(model) <- groups
      for (g in groups) {
        sym_g <- intersect(syms, GM$Symbol[GM$Group == g])
        if (!length(sym_g)) next
        D <- build_is_design(sym_subset = sym_g)
        if (!is.null(D$Xis) && length(D$yis) >= min_train_samples) {
          scal <- normalize_fit(D$Xis)
          model[[g]] <- list(obj = fit_fn(scal$transform(D$Xis), D$yis), scal = scal)
        }
      }
    }

    # OOS scoring
    for (t in (is_end + 1L):oos_end) {
      for (nm in syms) {
        xv <- seq_vec_at(t, nm); if (is.null(xv)) next
        sc <- NA_real_
        if (group == "pooled") {
          if (!is.null(model)) sc <- as.numeric(predict_fn(model$obj, model$scal$transform(matrix(xv, nrow = 1L))))
        } else if (group == "per_symbol") {
          mdl <- model[[nm]]
          if (!is.null(mdl)) sc <- as.numeric(predict_fn(mdl$obj, mdl$scal$transform(matrix(xv, nrow = 1L))))
        } else {
          g <- sym2grp[[nm]]; mdl <- model[[g]]
          if (!is.null(mdl)) sc <- as.numeric(predict_fn(mdl$obj, mdl$scal$transform(matrix(xv, nrow = 1L))))
        }
        S[[nm]][t] <- sc
      }
    }

    i <- i + step
  }

  S
}



#' One-call backtest wrapper (sequence features)
#' @inheritParams roll_fit_predict_seq
#' @inheritParams ml_backtest
#' @return list: `scores`, `mask`, `weights`, `backtest`.
#' @export
#' @examples
#' \donttest{
#' # as above, but with steps/horizon and normalize
#' }
ml_backtest_seq <- function(features_list,
                            labels,
                            steps = 26L, horizon = 4L,
                            fit_fn, predict_fn,
                            schedule        = list(is = 156L, oos = 4L, step = 4L),
                            group           = c("pooled","per_symbol","per_group"),
                            group_map       = NULL,
                            normalize       = c("none","zscore","minmax"),
                            selection       = list(top_k = 15L, max_per_group = NULL),
                            weighting       = list(method = "softmax", temperature = 12, floor = 0),
                            caps            = list(max_per_symbol = 0.10, max_per_group = NULL),
                            prices,
                            initial_capital = 1e5,
                            name            = "SEQ backtest") {
  group     <- match.arg(group)
  normalize <- match.arg(normalize)
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # roll sequence scores
  S <- roll_fit_predict_seq(features_list, labels,
                            steps = steps, horizon = horizon,
                            fit_fn = fit_fn, predict_fn = predict_fn,
                            is_periods = schedule$is, oos_periods = schedule$oos, step = schedule$step,
                            group = group, group_map = group_map,
                            normalize = normalize)

  # optional transform
  Sc <- transform_scores(S, method = "zscore")

  # selection (with optional group cap)
  mask <- if (!is.null(selection$max_per_group) && !is.null(group_map)) {
    select_top_k_scores_by_group(Sc, k = selection$top_k %||% 15L,
                                 group_map = group_map, max_per_group = selection$max_per_group)
  } else {
    select_top_k_scores(Sc, k = selection$top_k %||% 15L)
  }
  Sc_sel <- data.table::copy(Sc)
  syms <- setdiff(names(Sc_sel), "Date")
  for (nm in syms) Sc_sel[[nm]][!mask[[nm]]] <- NA_real_

  # weights
  W0 <- weight_from_scores(Sc_sel,
                           method      = weighting$method %||% "softmax",
                           temperature = weighting$temperature %||% 12,
                           floor       = weighting$floor %||% 0)
  W  <- cap_exposure(W0,
                     caps      = list(max_per_symbol = caps$max_per_symbol %||% 0.10,
                                      max_per_group  = caps$max_per_group %||% NULL),
                     group_map = group_map,
                     renorm    = TRUE)

  # align prices to S$Date & backtest
  Px <- if (data.table::is.data.table(prices)) data.table::copy(prices) else data.table::as.data.table(prices)
  if (!inherits(Px$Date, "Date")) Px[, Date := as.Date(Date)]
  data.table::setorderv(Px, "Date")
  Px_aligned <- Px[Px$Date %in% S$Date, c("Date", syms), with = FALSE]

  BT <- run_backtest(Px_aligned, W, initial_capital, name)
  list(scores = S, mask = mask, weights = W, backtest = BT)
}


# ==== LIGHTWEIGHT DIAGNOSTICS FOR SEQUENCES ===================================


#' Membership stability across dates
#' @param mask logical mask panel from selection.
#' @return data.table with `Date, stability` (Jaccard vs previous date).
#' @export
membership_stability <- function(mask) {
  DT <- if (data.table::is.data.table(mask)) data.table::copy(mask) else data.table::as.data.table(mask)
  syms <- setdiff(names(DT), "Date")
  if (nrow(DT) < 2) return(DT[, .(Date, stability = NA_real_)])
  stab <- rep(NA_real_, nrow(DT))
  prev <- as.logical(DT[1, ..syms])
  for (r in 2:nrow(DT)) {
    cur <- as.logical(DT[r, ..syms])
    both <- sum(prev | cur)
    kept <- sum(prev & cur)
    stab[r] <- if (both > 0) kept / both else NA_real_
    prev <- cur
  }
  data.table::data.table(Date = DT$Date, stability = stab)
}



#' Turnover by date
#' @param weights weight panel.
#' @return data.table with `Date, turnover` (0.5 * L1 change).
#' @export
turnover_by_date <- function(weights) {
  DT <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(DT), "Date")
  to <- rep(NA_real_, nrow(DT))
  prev <- rep(0, length(syms))
  for (r in seq_len(nrow(DT))) {
    row <- as.numeric(DT[r, ..syms]); row[!is.finite(row)] <- 0
    s <- sum(row); if (s > 0) row <- row / s
    to[r] <- 0.5 * sum(abs(row - prev))
    prev <- row
  }
  data.table::data.table(Date = DT$Date, turnover = to)
}

group_exposure <- function(weights, group_map) {
  DT <- if (data.table::is.data.table(weights)) data.table::copy(weights) else data.table::as.data.table(weights)
  syms <- setdiff(names(DT), "Date")
  GM <- validate_group_map(syms, group_map)
  groups <- unique(GM$Group)
  out <- data.table::data.table(Date = DT$Date)
  for (g in groups) out[[g]] <- NA_real_
  for (r in seq_len(nrow(DT))) {
    w <- as.numeric(DT[r, ..syms]); names(w) <- syms
    for (g in groups) {
      in_g <- GM$Symbol[GM$Group == g]
      out[[g]][r] <- sum(w[names(w) %in% in_g], na.rm = TRUE)
    }
  }
  out
}




# ==== PURGED/EMBARGOED CV FOR SEQUENCES ======================================

# Contiguous K-folds on decision indices with purge/embargo around validation.
purged_folds <- function(t_idx, k = 3L, purge = 0L, embargo = 0L) {
  t_idx <- sort(unique(as.integer(t_idx)))
  n <- length(t_idx)
  if (k < 2L || n < k) stop("Not enough samples for k folds")
  cuts <- floor(seq(0, n, length.out = k + 1L))
  folds <- vector("list", k)
  for (i in seq_len(k)) {
    vr  <- seq.int(cuts[i] + 1L, cuts[i + 1L])
    val <- t_idx[vr]
    lb  <- min(val) - purge
    ub  <- max(val) + embargo
    train <- t_idx[t_idx < lb | t_idx > ub]
    folds[[i]] <- list(train = train, val = val)
  }
  folds
}

# Internal: build sequences for an arbitrary set of decision indices (pooled).
.seq_build_from_indices <- function(F_list, Lbl, t_idx, syms, steps,
                                    na_action = c("omit","zero")) {
  na_action <- match.arg(na_action)
  if (!length(t_idx)) return(NULL)
  p <- length(F_list)
  pairs <- expand.grid(t = t_idx, s = syms, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  S <- nrow(pairs)
  X3 <- array(NA_real_, dim = c(S, steps, p))
  y  <- rep(NA_real_, S)
  keep <- rep(TRUE, S)
  for (i in seq_len(S)) {
    t <- pairs$t[i]; s <- pairs$s[i]
    if ((t - steps + 1L) < 1L) { keep[i] <- FALSE; next }
    ok_sample <- TRUE
    for (j in seq_len(p)) {
      v <- F_list[[j]][[s]][(t - steps + 1L):t]
      if (any(!is.finite(v))) {
        if (na_action == "omit") { ok_sample <- FALSE; break }
        v[!is.finite(v)] <- 0
      }
      X3[i, , j] <- as.numeric(v)
    }
    yy <- as.numeric(Lbl[[s]][t])
    if (!is.finite(yy) && na_action == "omit") ok_sample <- FALSE
    if (!ok_sample) { keep[i] <- FALSE } else { y[i] <- if (is.finite(yy)) yy else 0 }
  }
  if (!any(keep)) return(NULL)
  list(X3 = X3[keep, , , drop = FALSE],
       y  = y [keep],
       pairs = pairs[keep, , drop = FALSE])
}



#' Purged/embargoed K-fold CV for sequence models (inside IS window)
#'
#' @param features_list List of feature panels (Date + symbols) for sequences.
#' @param labels        Label panel aligned to features (Date + symbols).
#' @param is_start,is_end Inclusive IS window indices (on aligned dates).
#' @param steps_grid    Integer vector of candidate sequence lengths.
#' @param horizon       Forecast horizon (periods ahead).
#' @param fit_fn        Function (X, y) -> model for sequences.
#' @param predict_fn    Function (model, Xnew) -> numeric scores.
#' @param k             Number of CV folds.
#' @param purge         Gap (obs) removed between train/val within folds (default uses steps).
#' @param embargo       Embargo (obs) after validation to avoid bleed (default uses horizon).
#' @param group         'pooled', 'per_symbol', or 'per_group'.
#' @param max_train_samples Optional cap on IS samples per fold.
#' @param max_val_samples   Optional cap on validation samples per fold.
#' @param na_action     How to handle NA features ('omit' or 'zero').
#' @param metric        IC metric ('spearman' or 'pearson').
#'
#' @return data.table with columns like `steps`, `folds`, and CV `score`.
#' @export
cv_tune_seq <- function(features_list, labels,
                        is_start, is_end,
                        steps_grid = c(12L, 26L, 52L),
                        horizon = 4L,
                        fit_fn, predict_fn,
                        k = 3L, purge = NULL, embargo = NULL,
                        group = c("pooled"),
                        max_train_samples = NULL, max_val_samples = NULL,
                        na_action = c("omit","zero"),
                        metric = c("spearman","rmse")) {
  group    <- match.arg(group)
  stopifnot(group == "pooled")  # per-symbol CV omitted to keep it light
  na_action<- match.arg(na_action)
  metric   <- match.arg(metric)

  ALN <- .align_features_labels(features_list, labels)
  F_list <- ALN$features; Lbl <- ALN$labels
  dates <- ALN$dates;     syms <- ALN$symbols
  n <- length(dates)

  out <- data.table::data.table(steps = steps_grid, folds = NA_integer_,
                                score = NA_real_, metric = metric)

  for (g in seq_along(steps_grid)) {
    steps <- as.integer(steps_grid[g])
    # earliest usable t must be >= steps
    is_lo <- max(is_start, steps)
    is_hi <- min(is_end, n - horizon)
    if (is_hi < is_lo) { next }
    t_idx <- seq.int(is_lo, is_hi)

    # default purge/embargo
    pur <- if (is.null(purge)) steps else purge
    emb <- if (is.null(embargo)) horizon else embargo

    FOLDS <- purged_folds(t_idx, k = k, purge = pur, embargo = emb)
    fold_scores <- c()

    for (fd in FOLDS) {
      TR <- .seq_build_from_indices(F_list, Lbl, fd$train, syms, steps, na_action)
      VL <- .seq_build_from_indices(F_list, Lbl, fd$val,   syms, steps, na_action)
      if (is.null(TR) || is.null(VL)) next

      # optional downsampling for RAM/time control
      if (!is.null(max_train_samples) && nrow(TR$pairs) > max_train_samples) {
        keep <- sample.int(nrow(TR$pairs), max_train_samples)
        TR$X3 <- TR$X3[keep, , , drop = FALSE]; TR$y <- TR$y[keep]
      }
      if (!is.null(max_val_samples) && nrow(VL$pairs) > max_val_samples) {
        keep <- sample.int(nrow(VL$pairs), max_val_samples)
        VL$X3 <- VL$X3[keep, , , drop = FALSE];  VL$y <- VL$y[keep]
      }

      Xtr <- .flatten_seq3d_to_2d(TR$X3)
      Xvl <- .flatten_seq3d_to_2d(VL$X3)

      scaler <- .compute_scaler_2d(Xtr)
      Xtr_s  <- .apply_scaler_2d(Xtr, scaler)
      Xvl_s  <- .apply_scaler_2d(Xvl, scaler)

      model <- fit_fn(Xtr_s, TR$y)
      pred  <- as.numeric(predict_fn(model, Xvl_s))

      if (metric == "spearman") {
        sc <- suppressWarnings(stats::cor(pred, VL$y, method = "spearman", use = "complete.obs"))
      } else { # rmse
        sc <- sqrt(mean((pred - VL$y)^2, na.rm = TRUE))
        sc <- -sc  # higher is better
      }
      if (is.finite(sc)) fold_scores <- c(fold_scores, sc)
    }

    if (length(fold_scores)) {
      out$folds[g] <- length(fold_scores)
      out$score[g] <- mean(fold_scores)
    }
  }

  out[order(-score)]
}



# ==== WALK-FORWARD ENUMERATION & SWEEP (TABULAR) ==============================

.wf_enumerate_windows <- function(dates, is, oos, step) {
  n <- length(dates)
  out <- data.table::data.table(
    is_start = integer(), is_end = integer(),
    oos_start = integer(), oos_end = integer()
  )
  i <- 1L
  while (TRUE) {
    is_start <- i
    is_end   <- i + is - 1L
    oos_end  <- is_end + oos
    if (oos_end > n) break
    out <- rbind(out, data.table::data.table(
      is_start = is_start, is_end = is_end,
      oos_start = is_end + 1L, oos_end = oos_end
    ))
    i <- i + step
  }
  out
}

.expand_grid_list <- function(grid) {
  stopifnot(is.list(grid), length(grid) >= 1)
  keys <- names(grid)
  comb <- do.call(expand.grid, c(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  data.table::as.data.table(comb)
}

# Compute per-window OOS metrics from Scores/Weights/Labels/Prices slices
.wf_window_metrics <- function(scores, labels, weights, prices, oos_dates, freq = 52,
                               ic_method = "spearman") {
  # Slice to OOS date range
  in_rng <- scores$Date %in% oos_dates
  S <- scores[in_rng]
  L <- labels[labels$Date %in% oos_dates]
  W <- weights[weights$Date %in% oos_dates]
  P <- prices[prices$Date %in% oos_dates]

  # IC over OOS dates
  ics <- ic_series(S, L, method = ic_method)
  ic_mean <- mean(ics$IC, na.rm = TRUE)

  # Portfolio returns (optionally cost-aware; call site can pre-apply)
  rets <- portfolio_returns(W, P)
  r <- rets$ret
  mu <- mean(r, na.rm = TRUE) * freq
  sd <- stats::sd(r, na.rm = TRUE) * sqrt(freq)
  sharpe <- if (is.finite(sd) && sd > 0) mu / sd else NA_real_

  data.table::data.table(ic = ic_mean, ann_return = mu, sharpe = sharpe)
}


#' Walk-forward sweep of tabular configs (window-wise distribution of metrics)
#' @param max_windows optional limit for speed.
#' @param features_list List of feature panels.
#' @param labels Label panel.
#' @param prices Price panel for backtests.
#' @param fit_fn,predict_fn Model fit and predict functions.
#' @param schedule List with elements is, oos, step.
#' @param grid Named list of parameter vectors to sweep (e.g., top_k, temperature, method, transform).
#' @param caps Exposure caps (e.g., max_per_symbol/max_per_group).
#' @param group_map Optional Symbol->Group mapping for group caps/selection.
#' @param freq Compounding frequency for annualization (e.g., 52 for weekly).
#' @param cost_bps One-way trading cost in basis points (applied on rebalance).
#' @param ic_method IC method ('spearman' or 'pearson').
#' @return data.table with medians/means across OOS windows.
#' @export
wf_sweep_tabular <- function(features_list, labels, prices,
                             fit_fn, predict_fn,
                             schedule = list(is = 104L, oos = 4L, step = 4L),
                             grid = list(top_k = c(10L,15L),
                                         temperature = c(8,12),
                                         method = c("softmax","rank"),
                                         transform = c("zscore")),
                             caps = list(max_per_symbol = 0.08, max_per_group = NULL),
                             group_map = NULL,
                             freq = 52,           # periods/year for annualization
                             cost_bps = 0,        # cost applied when computing returns, if > 0
                             max_windows = NULL,  # throttle for speed
                             ic_method = "spearman") {

  # Align once
  ALN <- .align_features_labels(features_list, labels)
  F_list <- ALN$features
  Lbl    <- ALN$labels
  Px     <- if (data.table::is.data.table(prices)) data.table::copy(prices) else data.table::as.data.table(prices)
  Px     <- Px[Px$Date %in% ALN$dates, c("Date", ALN$symbols), with = FALSE]

  # Enumerate WF windows
  WIN <- .wf_enumerate_windows(ALN$dates, schedule$is, schedule$oos, schedule$step)
  if (!is.null(max_windows) && nrow(WIN) > max_windows) {
    WIN <- WIN[seq_len(max_windows)]
  }

  # Expand grid
  CFG <- .expand_grid_list(grid)
  out <- data.table::data.table()
  if (!nrow(CFG) || !nrow(WIN)) return(out)

  # Precompute raw scores (pooled) per feature -> then combine (mean)
  #  wrap each feature panel in a list()
  scores_per_feat <- lapply(F_list, function(FX) {
    roll_fit_predict(list(FX), Lbl, fit_fn, predict_fn,
                     is_periods = schedule$is, oos_periods = schedule$oos, step = schedule$step)
  })

  # Combine scores: simple mean across features
  S <- data.table::data.table(Date = ALN$dates)
  for (sym in ALN$symbols) {
    mat <- do.call(cbind, lapply(scores_per_feat, function(sc) as.numeric(sc[[sym]])))
    S[[sym]] <- rowMeans(mat, na.rm = TRUE)
  }

  # For each config, transform/select/weight once; then slice by windows for metrics
  for (i in seq_len(nrow(CFG))) {
    cfg <- CFG[i]
    Sx <- data.table::copy(S)

    # transform
    trf <- as.character(cfg$transform)
    if (!is.na(trf) && trf != "none") Sx <- transform_scores(Sx, trf)

    # selection
    top_k <- as.integer(cfg$top_k)
    if (!is.null(caps$max_per_group) && !is.null(group_map)) {
      mask <- select_top_k_scores_by_group(Sx, k = top_k, group_map = group_map,
                                           max_per_group = caps$max_per_group)
    } else {
      mask <- select_top_k_scores(Sx, k = top_k)
    }
    Ssel <- data.table::copy(Sx)
    syms <- ALN$symbols
    for (nm in syms) Ssel[[nm]][!mask[[nm]]] <- NA_real_

    # weighting
    method <- as.character(cfg$method)
    temp   <- if ("temperature" %in% names(cfg)) as.numeric(cfg$temperature) else 12
    W0 <- weight_from_scores(Ssel, method = method, temperature = temp)
    W  <- cap_exposure(W0, caps = caps, group_map = group_map, renorm = TRUE)

    # Optional costs: fold into a per-date net return when computing metrics
    # We'll pass raw weights here; .wf_window_metrics will call portfolio_returns(W, Px)
    # If you want costs, compute returns here and overwrite W with an equivalent "net ret" path,
    # but simpler: keep cost outside this sweep or use cost-aware backtests separately.

    # Per-window metrics
    rows <- vector("list", nrow(WIN))
    for (w in seq_len(nrow(WIN))) {
      rng <- WIN[w]
      oos_dates <- ALN$dates[rng$oos_start:rng$oos_end]
      # if costs requested, create a temporary cost-adjusted return series
      if (cost_bps > 0) {
        rets <- portfolio_returns(W[W$Date %in% oos_dates], Px[Px$Date %in% oos_dates], cost_bps = cost_bps)
        # compute window metrics from cost-aware returns
        mu <- mean(rets$ret, na.rm = TRUE) * freq
        sd <- stats::sd(rets$ret, na.rm = TRUE) * sqrt(freq)
        sharpe <- if (is.finite(sd) && sd > 0) mu / sd else NA_real_
        ics <- ic_series(S[S$Date %in% oos_dates], Lbl[Lbl$Date %in% oos_dates], method = ic_method)
        ic_mean <- mean(ics$IC, na.rm = TRUE)
        rows[[w]] <- data.table::data.table(ic = ic_mean, ann_return = mu, sharpe = sharpe)
      } else {
        rows[[w]] <- .wf_window_metrics(S, Lbl, W, Px, oos_dates, freq = freq, ic_method = ic_method)
      }
    }
    M <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)

    out <- rbind(out, data.table::data.table(
      top_k       = top_k,
      temperature = temp,
      method      = method,
      transform   = trf,
      n_windows   = nrow(M),
      ic_median   = stats::median(M$ic, na.rm = TRUE),
      ic_mean     = base::mean(M$ic, na.rm = TRUE),
      annret_med  = stats::median(M$ann_return, na.rm = TRUE),
      annret_mean = base::mean(M$ann_return, na.rm = TRUE),
      sharpe_med  = stats::median(M$sharpe, na.rm = TRUE),
      sharpe_mean = base::mean(M$sharpe, na.rm = TRUE)
    ), fill = TRUE)
  }

  data.table::setorder(out, -sharpe_med, -annret_med)
  out
}





