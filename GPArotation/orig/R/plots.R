#' Plot Method for GPArotation Objects
#'
#' @param x An object of class "GPArotation"
#' @param type Character string indicating plot type: "salient", "pairs", or "profile"
#' @param hide_hyperplane Logical. If TRUE, blanks out values < cutoff in the salient plot
#'                                 or hides the legend in the profile plot.
#' @param salient Numeric. Values above this are considered salient.
#' @param cutoff Numeric. Values below this are considered in the hyperplane.
#' @param ... Additional arguments passed to base plot functions.
#'
#' @export

plot.GPArotation <- function(x,
            type = c("salient", "pairs", "profile", "vector", "target", 
                     "heatmap","trajectory", "diagnostics", "residuals"),
            hide_hyperplane = FALSE, salient = 0.40, cutoff = .1, R = NULL,
            factors = c(1,2), items = NULL, main = NULL, ...){
    
  type <- match.arg(type)
  L   <- unclass(x$loadings)
  Phi <- x$Phi
 
  k <- ncol(L)
  v <- nrow(L)
 
  fmt_num <- function(val) {
    txt <- sprintf("%.2f", val)
    txt <- sub("^-0.", "-.", txt)
    txt <- sub("^0.", ".", txt)
    return(txt)
  }
 
  # ---------------------------------------------------------
  # TYPE 1: THE APA SALIENT DASHBOARD
  # ---------------------------------------------------------
  if (type == "salient") {

    # Sort for display only
    x_sorted <- .sortGPALoadings(x)
    L   <- unclass(x_sorted$loadings)
    Phi <- x_sorted$Phi

    old_par <- par(mar = c(2, 5, 4, 1), oma = c(2, 0, 2, 0))
    on.exit(par(old_par))
   
    if (!is.null(Phi)) {
      S <- L %*% Phi
      rownames(S) <- rownames(L)
      colnames(S) <- colnames(L)
     
      lay_mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
      layout(lay_mat, heights = c(v + 2, max(k * 3, 8)), widths = c(1, 1))
      # layout(lay_mat, heights = c(v + 2, k + 2), widths = c(1, 1))
    } else {
      layout(matrix(1, nrow = 1))
    }
   
    breaks <- c(0, cutoff, salient, max(abs(L), 1) + 0.1)
    cols   <- c("white", "grey95", "grey80")
   
    draw_matrix <- function(mat, title, exempt_diag = FALSE) {
      nr <- nrow(mat)
      nc <- ncol(mat)
      mat_rev <- mat[nr:1, , drop = FALSE]
     
      text_cex <- min(1.2, 25 / max(nr, nc))
      axis_cex <- min(1.0, text_cex * 1.1)

      # For shading: cap diagonal at 0 so it stays white
      mat_shade <- abs(mat_rev)
      if (exempt_diag) {
        for (d in 1:min(nr, nc))
          mat_shade[nr - d + 1, d] <- 0
      }
     
      image(x = 1:nc, y = 1:nr, z = t(mat_shade),
            col = cols, breaks = breaks,
            axes = FALSE, xlab = "", ylab = "")
     
	  if (nchar(title) > 0)
         mtext(title, side = 3, line = 2.5, font = 2, cex = min(1.0, axis_cex * 1.1))

      axis(3, at = 1:nc, labels = colnames(mat), tick = FALSE,
           line = 0.5, cex.axis = axis_cex)
      axis(2, at = 1:nr, labels = rownames(mat_rev), tick = FALSE,
           las = 2, line = -0.5, cex.axis = axis_cex)
     
      box()
      if (nc > 1) abline(v = (1:(nc-1)) + 0.5, col = "grey90")
      if (nr > 1) abline(h = (1:(nr-1)) + 0.5, col = "grey90")
     
      for (i in 1:nr) {
        for (j in 1:nc) {
          val <- mat_rev[i, j]
          if (hide_hyperplane && abs(val) < cutoff) next
          # Exempt diagonal from bolding
          is_diag    <- exempt_diag && (i == (nr - j + 1))
          font_style <- ifelse(!is_diag && abs(val) >= salient, 2, 1)
          text(j, i, labels = fmt_num(val), font = font_style, cex = text_cex)
        }
      }
    }
   
    # 1. Pattern Matrix (Top Left)
    if (is.null(Phi)) 
      draw_matrix(L, "")
   
    if (!is.null(Phi)) {
      # 1. Pattern Matrix (Top Left)
      draw_matrix(L, "Pattern Matrix")

      # 2. Structure Matrix (Top Right)
      draw_matrix(S, "Structure Matrix")
     
      # 3. Phi Matrix (Bottom Left) - exempt diagonal
      dimnames(Phi) <- list(colnames(L), colnames(L))
      draw_matrix(Phi, "", exempt_diag = TRUE)
     
      # 4. Fit stats (Bottom Right)
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      box(col = "grey80")

      stats_list <- c()
      if (!is.null(x$n.obs)) stats_list <- c(stats_list, paste("N =", x$n.obs))
      if (!is.null(x$dof))   stats_list <- c(stats_list, paste("df =", x$dof))
      if (!is.null(x$RMSEA)) stats_list <- c(stats_list, paste("RMSEA =", fmt_num(x$RMSEA)))
      if (!is.null(x$SRMR))  stats_list <- c(stats_list, paste("SRMR =",  fmt_num(x$SRMR)))
     
      if (length(stats_list) > 0) {
        stats_text <- paste(stats_list, collapse = "\n")
        text(0.5, 0.5, labels = stats_text, cex = 1.3, font = 2)
      }
    }
  }

  # ---------------------------------------------------------
  # TYPE 2: PAIRS PLOT
  # ---------------------------------------------------------
  else if (type == "pairs") {
    old_par <- par(mfrow = c(k, k), mar = c(2, 2, 1, 1), oma = c(4, 4, 4, 1))
    on.exit(par(old_par))

    # Consistent y-axis across all histograms
    max_count <- max(sapply(1:k, function(i)
      max(hist(L[, i], plot = FALSE, breaks = 15)$counts)))
    hist_ylim <- c(0, max_count * 1.1)

    for (row_idx in 1:k) {
      for (col_idx in 1:k) {
        if (row_idx == col_idx) {
          hist(L[, row_idx], main = "", xlab = "", ylab = "",
               col = "lightblue", border = "white", breaks = 15,
               xlim = c(-1, 1), ylim = hist_ylim)
          legend("topleft", bty = "n", legend = colnames(L)[row_idx], text.font = 2)
          box(col = "grey80")
          } else if (row_idx < col_idx) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            if (!is.null(Phi)) {
              text(0,  0.15, labels = paste0("r = ",  fmt_num(Phi[row_idx, col_idx])),
                   cex = 1.5, font = 2)
              text(0, -0.15, labels = paste0("r\u00B2 = ", fmt_num(Phi[row_idx, col_idx]^2)),
                   cex = 1.2, font = 1)
              # TODO: additional factor relationship summaries
             } else {
               text(0, 0, labels = "Orthogonal", col = "grey", cex = 1.2)
             }
            box(col = "grey80")
          }
          else {
          plot(L[, col_idx], L[, row_idx], pch = 20, col = "black",
               xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
          abline(h = 0, v = 0, col = "grey50", lty = 2)
          box(col = "grey80")
        }
      }
    }
  }
  
  # ---------------------------------------------------------
  # TYPE 3: PROFILE BAR CHARTS
  # ---------------------------------------------------------
  else if (type == "profile") {
    old_par <- par(mfrow = c(k, 1), mar = c(1, 4, 1, 1), oma = c(6, 0, 2, 0))
    on.exit(par(old_par))
    for (i in 1:k) {
      bar_cols <- ifelse(abs(L[, i]) >= salient, "steelblue", "grey80")
      bp <- barplot(L[, i], names.arg = rep("", v), ylim = c(-1, 1),
                    col = bar_cols, border = NA, ylab = colnames(L)[i], las = 2)
      abline(h = 0, col = "black")
      abline(h = c(-cutoff, cutoff), col = "red", lty = 2)
      if (!hide_hyperplane)
        legend("topright", bty = "o", cex = 0.8, bg = "white",
               legend = c(paste("hyperplane:", sum(abs(L[, i]) < cutoff)),
                          paste("salient:   ", sum(abs(L[, i]) >= salient))))
      if (i == k) {
        axis(1, at = bp, labels = rownames(L), las = 2, cex.axis = 0.8)
        mtext("Item", side = 1, line = 5, cex = 0.8)
      }
    }
  }
  
  # ---------------------------------------------------------
  # TYPE 4: TARGET PLOTS
  # ---------------------------------------------------------
  else if (type == "target") {

    Target <- list(...)$Target
    if (is.null(Target) && !is.null(x$Target)) Target <- x$Target

    if (is.null(Target)) {
      warning("No Target matrix found. Please supply one ",
              "(e.g., plot(x, type='target', Target=my_matrix)).")
    } else {

      A_unrot <- attr(x, "A_unrotated")
      L_rot   <- L

      f1 <- if (!is.null(list(...)$f1)) list(...)$f1 else 1
      f2 <- if (!is.null(list(...)$f2)) list(...)$f2 else 2

      # Pre-rotation positions
      Tinit <- if (!is.null(x$Tinit)) x$Tinit else diag(ncol(A_unrot))
      L_pre <- if (x$orthogonal) {
        A_unrot %*% Tinit
      } else {
        A_unrot %*% t(solve(Tinit))
      }

      old_par <- par(mfrow = c(1, 2), mar = c(4, 5, 4, 1))
      on.exit(par(old_par), add = TRUE)

      # ----- Panel 1: Arrow plot -----
      f1 <- if (!is.null(list(...)$f1)) list(...)$f1 else 1
      f2 <- if (!is.null(list(...)$f2)) list(...)$f2 else 2

      Tinit <- if (!is.null(x$Tinit)) x$Tinit else diag(ncol(A_unrot))
      L_pre <- if (x$orthogonal) {
        A_unrot %*% Tinit
      } else {
        A_unrot %*% t(solve(Tinit))
      }

      target_complete <- !apply(Target[, c(f1, f2), drop = FALSE], 1, anyNA)
      n_specified     <- sum(target_complete)

      if (n_specified < 2) {
        # PST with few fully specified rows --- show message instead
        plot.new()
        text(0.5, 0.5,
             paste0("Target partially specified\n(",
                    n_specified, " fully specified rows in factors ",
                    f1, " and ", f2, ").\nSee heatmap for discrepancies."),
             cex = 0.9, col = "grey40")
      } else {
        # Full or mostly specified target --- draw arrow plot
        all_coords <- rbind(L_pre[, c(f1, f2)],
                            L_rot[, c(f1, f2)],
                            Target[target_complete, c(f1, f2)])
        lims <- range(all_coords, na.rm = TRUE) * 1.15
        plot(NULL,
             xlim = lims, ylim = lims,
             xlab = paste("Factor", f1),
             ylab = paste("Factor", f2),
             main = if (is.null(main)) "Target Rotation: Item Trajectories" else main,
             asp  = 1)
        abline(h = 0, v = 0, col = "grey70", lty = 2)

        arrows(x0 = L_pre[, f1], y0 = L_pre[, f2],
               x1 = L_rot[, f1], y1 = L_rot[, f2],
               length = 0.06, col = "grey60", lwd = 0.8)

        points(Target[target_complete, f1],
               Target[target_complete, f2],
               pch = 19, col = "steelblue", cex = 0.8)

        points(L_pre[, f1], L_pre[, f2],
               pch = 17, col = "tomato", cex = 0.8)

        points(L_rot[, f1], L_rot[, f2],
               pch = 15, col = "gold3", cex = 0.8)

        legend("topright", bty = "n", cex = 0.75,
               legend = c("Target", "Pre-rotation", "Post-rotation"),
               col    = c("steelblue", "tomato", "gold3"),
               pch    = c(19, 17, 15))
      }

      # ----- Panel 2: Discrepancy heatmap -----
      Discrepancy <- L_rot - Target
      nr    <- nrow(Discrepancy)
      nc    <- ncol(Discrepancy)
      D_rev <- Discrepancy[nr:1, , drop = FALSE]

      error_cols  <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
      breaks_disc <- seq(-1, 1, length.out = 101)

      image(x = 1:nc, y = 1:nr, z = t(D_rev),
            col = error_cols, breaks = breaks_disc,
            axes = FALSE, xlab = "", ylab = "",
            main = "Rotated minus Target")

      text_cex <- min(1.2, 25 / max(nr, nc))
      axis(3, at = 1:nc, labels = colnames(Discrepancy), tick = FALSE,
           line = -0.5, cex.axis = min(1.0, text_cex * 1.1))
      axis(2, at = 1:nr, labels = rownames(D_rev), tick = FALSE,
           las = 2, line = -0.5, cex.axis = min(1.0, text_cex * 1.1))
      box()
      if (nc > 1) abline(v = (1:(nc-1)) + 0.5, col = "grey90")
      if (nr > 1) abline(h = (1:(nr-1)) + 0.5, col = "grey90")

      for (i in 1:nr) {
        for (j in 1:nc) {
          val <- D_rev[i, j]
          if (is.na(val)) {
            text(j, i, labels = "-", cex = text_cex, col = "grey70")
          } else {
            text_col <- ifelse(abs(val) > 0.4, "white", "black")
            text(j, i, labels = fmt_num(val), cex = text_cex, col = text_col)
          }
        }
      }
    }
  }
     
  # ---------------------------------------------------------
  # TYPE 5: VECTOR PLOTS
  # ---------------------------------------------------------
  else if (type == "vector") {

    old_par <- par(mar = c(3, 3, 3, 1), oma = c(2, 2, 2, 1))
    on.exit(par(old_par))

    # (k-1) x (k-1) grid - lower triangle pairs only
    par(mfrow = c(k - 1, k - 1))

    for (row_idx in 2:k) {
      for (col_idx in 1:(k - 1)) {

        if (col_idx >= row_idx) {
          # Empty panel for upper triangle
          plot.new()
        } else {
          # r value as panel title
          r_label <- if (!is.null(Phi))
            paste0("r = ", fmt_num(Phi[row_idx, col_idx]),
                   "  r\u00B2 = ", fmt_num(Phi[row_idx, col_idx]^2))
          else
            "Orthogonal"

          # Vector plot
          plot(0, 0, type = "n", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
               axes = FALSE, xlab = "", ylab = "", main = r_label)
          abline(h = 0, v = 0, col = "grey50", lty = 2)
          symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "grey80")
          arrows(x0 = 0, y0 = 0,
                 x1 = L[, col_idx], y1 = L[, row_idx],
                 length = 0.05, col = "steelblue", lwd = 1.2)
          text(L[, col_idx], L[, row_idx],
               labels = rownames(L),
               pos = 3, cex = 0.6, col = "grey30")
          box(col = "grey80")

          # Factor names as axis labels
          mtext(colnames(L)[col_idx], side = 1, line = 1.5, cex = 0.9, font = 2)
          mtext(colnames(L)[row_idx], side = 2, line = 1.5, cex = 0.9, font = 2)
        }
      }
    }
  }
  # ---------------------------------------------------------
  # TYPE 6: HEATMAP
  # ---------------------------------------------------------
  else if (type == "heatmap") {
    layout(matrix(1, nrow = 1))  # reset layout

    # Sort for display only
    x_sorted <- .sortGPALoadings(x)
    L   <- unclass(x_sorted$loadings)
    Phi <- x_sorted$Phi

    old_par <- par(mar = c(2, 5, 4, 1), oma = c(2, 0, 2, 0))
    on.exit(par(old_par))

    if (!is.null(Phi)) {
      S <- L %*% Phi
      rownames(S) <- rownames(L)
      colnames(S) <- colnames(L)

      lay_mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
      layout(lay_mat, heights = c(v + 2, max(k * 3, 8)), widths = c(1, 1))
    } else {
      layout(matrix(1, nrow = 1))
    }

    # Diverging palette: blue -> white -> red
    pal    <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
    breaks <- seq(-1, 1, length.out = 101)

    draw_heatmap <- function(mat, title, show_values = TRUE, exempt_diag = FALSE) {
      nr <- nrow(mat)
      nc <- ncol(mat)
      mat_rev <- mat[nr:1, , drop = FALSE]

      text_cex <- min(1.2, 25 / max(nr, nc))
      axis_cex <- min(1.0, text_cex * 1.1)

      # Clip values to [-1, 1] for color mapping
      mat_clipped <- pmax(pmin(mat_rev, 1), -1)

      # Exempt diagonal from coloring - set to 0 (white)
      if (exempt_diag) {
        for (d in 1:min(nr, nc))
          mat_clipped[nr - d + 1, d] <- 0
      }

      image(x = 1:nc, y = 1:nr, z = t(mat_clipped),
            col = pal, breaks = breaks,
            axes = FALSE, xlab = "", ylab = "")

      if (nchar(title) > 0)
        mtext(title, side = 3, line = 2.5, font = 2,
              cex = min(1.0, axis_cex * 1.1))

      axis(3, at = 1:nc, labels = colnames(mat), tick = FALSE,
           line = 0.5, cex.axis = axis_cex)
      axis(2, at = 1:nr, labels = rownames(mat_rev), tick = FALSE,
           las = 2, line = -0.5, cex.axis = axis_cex)

      box()
      if (nc > 1) abline(v = (1:(nc-1)) + 0.5, col = "grey70")
      if (nr > 1) abline(h = (1:(nr-1)) + 0.5, col = "grey70")

      if (show_values) {
        for (i in 1:nr) {
          for (j in 1:nc) {
            val        <- mat_rev[i, j]
            is_diag    <- exempt_diag && (i == (nr - j + 1))
            txt_col    <- ifelse(!is_diag && abs(val) > 0.6, "white", "black")
            font_style <- ifelse(!is_diag && abs(val) >= salient, 2, 1)
            text(j, i, labels = fmt_num(val),
                 col = txt_col, font = font_style, cex = text_cex)
          }
        }
      }
    }

    # 1. Pattern Matrix (Top Left)
    if (is.null(Phi))
      draw_heatmap(L, "")

    if (!is.null(Phi)) {
      # 1. Pattern Matrix (Top Left)
      draw_heatmap(L, "Pattern Matrix")

      # 2. Structure Matrix (Top Right)
      draw_heatmap(S, "Structure Matrix")

      # 3. Phi Matrix (Bottom Left)
      dimnames(Phi) <- list(colnames(L), colnames(L))
      draw_heatmap(Phi, "", show_values = TRUE, exempt_diag = TRUE)

      # 4. Fit stats (Bottom Right)
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      box(col = "grey80")

      stats_list <- c()
      if (!is.null(x$n.obs)) stats_list <- c(stats_list, paste("N =", x$n.obs))
      if (!is.null(x$dof))   stats_list <- c(stats_list, paste("df =", x$dof))
      if (!is.null(x$RMSEA)) stats_list <- c(stats_list, paste("RMSEA =", fmt_num(x$RMSEA)))
      if (!is.null(x$SRMR))  stats_list <- c(stats_list, paste("SRMR =",  fmt_num(x$SRMR)))

      if (length(stats_list) > 0) {
        stats_text <- paste(stats_list, collapse = "\n")
        text(0.5, 0.5, labels = stats_text, cex = 1.3, font = 2)
      }
    }
  }  
  else if (type == "trajectory") {
    plot_trajectory(x, factors = factors, items = items, main = main)
  }
  else if (type == "diagnostics"){
    plot_gpa_diagnostics(x, main = main)
  }
  else if (type == "residuals") {
    plot_residuals(x, R = R, cutoff = cutoff, main = main)
  }
  

  invisible(x)
}





