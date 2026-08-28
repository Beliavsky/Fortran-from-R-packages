plot_residuals <- function(x, R = NULL, cutoff = 0.10, main = NULL) {

  # --- Input validation ---
  if (!inherits(x, "GPArotation"))
    stop("'x' must be a GPArotation object.", call. = FALSE)

  # --- Correlation matrix ---
  if (is.null(R)) {
    if (!is.null(x$correlation)) {
      R <- x$correlation
    } else {
      stop("'R' is required. Supply the observed correlation matrix, ",
           "or rotate a factanal object so R is stored automatically.",
           call. = FALSE)
    }
  }
  if (!is.matrix(R))
    R <- as.matrix(R)
  R <- cov2cor(R)

  # --- Loadings in R row order ---
  L_res   <- unclass(x$loadings)
  Phi_res <- if (is.null(x$Phi)) diag(ncol(L_res)) else x$Phi
  L_res   <- L_res[rownames(R), , drop = FALSE]

  # --- Residual matrix ---
  Resid       <- residuals(x)
  nms         <- rownames(R)
  n           <- length(nms)
  srmr        <- sqrt(mean(Resid[lower.tri(Resid)]^2))

  old_par <- par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))
  on.exit(par(old_par))

  # --- Panel 1: residual heatmap via rect() ---
  max_abs <- max(abs(Resid))
  pal     <- colorRampPalette(c("steelblue", "white", "tomato"))(100)
  breaks  <- seq(-max_abs, max_abs, length.out = 101)

  main_title <- if (is.null(main))
    paste0("Residual Correlation Matrix\nSRMR = ", round(srmr, 3))
  else
    main

  plot(NULL,
       xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       main = main_title)
       
  for (i in 1:n) {
    for (j in 1:n) {
      col_idx <- findInterval(Resid[i, j], breaks)
      col_idx <- max(1, min(col_idx, 100))
      rect(j - 0.5, n - i + 0.5,
           j + 0.5, n - i + 1.5,
           col = pal[col_idx], border = NA)
    }
  }

  # grid lines
  abline(h = 0.5:(n + 0.5), col = "white", lwd = 0.5)
  abline(v = 0.5:(n + 0.5), col = "white", lwd = 0.5)
  box()

  # x-axis: CCAI8 at left, same order as R
  axis(1, at = 1:n, labels = nms,      las = 2, cex.axis = 0.7)
  # y-axis: CCAI8 at top (row 1 is drawn at top via n - i + 1)
  axis(2, at = 1:n, labels = rev(nms), las = 1, cex.axis = 0.7)

  # --- Panel 2: sorted absolute residuals ---
  res_vals <- sort(abs(Resid[lower.tri(Resid)]))
  bar_cols <- ifelse(res_vals > cutoff, "tomato",
              ifelse(res_vals > 0.05,   "darkgreen",
                                        "lightgreen"))

  barplot(res_vals,
          col    = bar_cols,
          border = NA,
          main   = "Sorted Absolute Residuals",
          ylab   = "Absolute residual",
          xlab   = "Item pair (sorted)")
  abline(h = cutoff, lty = 2, col = "gray40", lwd = 1.5)
  abline(h = 0.05,   lty = 3, col = "gray60", lwd = 1.0)
  legend("topleft", bty = "n", cex = 0.8,
         legend = c(paste0("> ", cutoff,  " (concern)"),
                    paste0("0.05 - ", cutoff, " (watch)"),
                    "< 0.05 (acceptable)"),
         fill   = c("tomato", "darkgreen", "lightgreen"),
         border = NA)

  invisible(list(
    Resid  = Resid,
    SRMR   = srmr,
    cutoff = cutoff
  ))
}


audit_residuals <- function(x, R = NULL, cutoff = 0.10, allDiff = FALSE) {

  # --- Input validation ---
  if (!inherits(x, "GPArotation"))
    stop("'x' must be a GPArotation object.", call. = FALSE)

  # --- Correlation matrix ---
  if (is.null(R)) {
    if (!is.null(x$correlation)) {
      R <- x$correlation
    } else {
      stop("'R' is required. Supply the observed correlation matrix, ",
           "or rotate a factanal object so R is stored automatically.",
           call. = FALSE)
    }
  }
  if (!is.matrix(R))
    R <- as.matrix(R)
  R <- cov2cor(R)

  # --- Loadings in R row order ---
  L_res   <- unclass(x$loadings)
  Phi_res <- if (is.null(x$Phi)) diag(ncol(L_res)) else x$Phi
  L_res   <- L_res[rownames(R), , drop = FALSE]

  # --- Residual matrix ---
  Resid <- residuals(x)

  # --- Summary statistics ---
  srmr       <- sqrt(mean(Resid[lower.tri(Resid)]^2))
  row_strain <- rowMeans(abs(Resid))

  # --- Flagged pairs ---
  idx   <- which(lower.tri(Resid), arr.ind = TRUE)
  pairs <- data.frame(
    Item_A   = rownames(L_res)[idx[, 1]],
    Item_B   = rownames(L_res)[idx[, 2]],
    Residual = round(Resid[lower.tri(Resid)], 4)
  )
  pairs <- pairs[abs(pairs$Residual) > cutoff, ]
  pairs <- pairs[order(abs(pairs$Residual), decreasing = TRUE), ]

  # --- Output ---
  cat(sprintf("SRMR: %.4f\n\n", srmr))
  cat("Mean absolute residual per item:\n")
  print(round(sort(row_strain, decreasing = TRUE), 4))

  if (nrow(pairs) == 0) {
    cat("\nNo item pairs with absolute residual exceeding", cutoff, "\n")
  } else {
    cat("\nItem pairs with absolute residual exceeding", cutoff, ":\n")
    print(pairs, row.names = FALSE)
  }

  if (allDiff) {
    cat("\nFull residual matrix (R - R_hat):\n")
    print(round(Resid, 4))
  }

  invisible(list(
    Resid      = Resid,
    SRMR       = srmr,
    row_strain = row_strain,
    pairs      = pairs,
    cutoff     = cutoff
  ))
}

