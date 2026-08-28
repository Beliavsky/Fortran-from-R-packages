
plot_trajectory <- function(res, factors = c(1, 2), items = NULL, main = NULL) {

  trajectory_data <- .reconstruct_trajectory(res)

  A_unrotated <- attr(res, "A_unrotated")
  p <- nrow(A_unrotated)
  if (is.null(items)) items <- 1:p

  # Build main title
  if (is.null(main)) {
    main_title <- paste("Factor Migration trajectory:", res$method)
  } else {
    main_title <- main  # "" suppresses title, any string overrides
  }

  plot(NULL, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), 
       xlab = paste("Factor", factors[1]),
       ylab = paste("Factor", factors[2]),
       main = main_title)

  # Draw manifold (Unit Circle)
  theta <- seq(0, 2*pi, length.out = 100)
  lines(cos(theta), sin(theta), col = "gray90")
  abline(h = 0, v = 0, col = "gray94", lty = 2)

  cols <- rainbow(length(items))

  for (i in seq_along(items)) {
    var_idx <- items[i]

    trajectory_coords <- t(sapply(trajectory_data, function(Ti) {
      if (res$orthogonal) {
        L_step <- A_unrotated %*% Ti
      } else {
        L_step <- A_unrotated %*% t(solve(Ti))
      }
      return(L_step[var_idx, factors])
    }))

##    points(trajectory_coords[1, 1], trajectory_coords[1, 2],
##           col = cols[i], pch = 4, cex = 0.5)
##    lines(trajectory_coords[, 1], trajectory_coords[, 2],
##          col = cols[i], lwd = 1)
##    points(trajectory_coords[nrow(trajectory_coords), 1],
##           trajectory_coords[nrow(trajectory_coords), 2],
##           col = cols[i], pch = 16, cex = 0.8)
##    text(trajectory_coords[nrow(trajectory_coords), 1],
##         trajectory_coords[nrow(trajectory_coords), 2],
##         labels = var_idx, pos = 3, cex = 0.7, col = cols[i])

    points(trajectory_coords[1, 1], trajectory_coords[1, 2],
           col = "black", pch = 4, cex = 0.5)
    lines(trajectory_coords[, 1], trajectory_coords[, 2],
          col = cols[i], lwd = 1)
    points(trajectory_coords[nrow(trajectory_coords), 1],
           trajectory_coords[nrow(trajectory_coords), 2],
           col = "black", pch = 16, cex = 0.8)
    text(trajectory_coords[nrow(trajectory_coords), 1],
         trajectory_coords[nrow(trajectory_coords), 2],
         labels = var_idx, pos = 3, cex = 0.7, col = "black")
  }

  stats_text <- paste0(
    "Algorithm: ", res$algorithm, "\n",
    "fwindow: ", res$fwindow, "\n",
    "Iterations: ", nrow(res$Table) - 1, "\n",
    "Final f: ", signif(res$Table[nrow(res$Table), 2], 3), "\n",
    "Final s: ", signif(res$Table[nrow(res$Table), 3], 3)
  )
  legend("topleft", legend = stats_text, bty = "n", cex = 0.8,
         adj = c(0, 0.5), inset = c(0, -0.05))
}

#############################################
## diagnostics plot
#############################################

plot_gpa_diagnostics <- function(res, trajectory_history = .reconstruct_trajectory(res), main = NULL) {
    # --- Input validation ---
    if (missing(res) || is.null(res))
        stop("'res' must be a GPArotation object.", call. = FALSE)
    if (!inherits(res, "GPArotation"))
        stop("'res' must be a GPArotation object, got: ", class(res), call. = FALSE)
    if (is.null(res$Table))
        stop("'res$Table' is missing - was this object created by GPForth or GPFoblq?",
             call. = FALSE)
    if (is.null(attr(res, "A_unrotated")))
        stop("'A_unrotated' attribute missing - object may have been created by a legacy function.",
             call. = FALSE)

    tab         <- data.frame(res$Table)
    A_unrotated <- attr(res, "A_unrotated")

    if (res$orthogonal) {
        mid_metric <- sapply(trajectory_history, function(Ti) {
            L <- A_unrotated %*% Ti
            var(as.vector(L^2))
        })
        mid_label <- "Var(Loadings^2)"
        mid_title <- "Simple Structure Emergence"
    } else {
        mid_metric <- sapply(trajectory_history, function(Ti) {
            Phi <- t(Ti) %*% Ti
            max(abs(Phi[lower.tri(Phi)]))
        })
        mid_label <- "Max |phi_ij|"
        mid_title <- "Factor Redundancy Trace"
    }

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    layout(matrix(c(1, 2, 3), ncol = 1), heights = c(2, 1, 1))

    # 1. Top panel - criterion gap
    f_gap <- pmax(tab$f - min(tab$f), .Machine$double.eps^0.5)
    par(mar = c(0, 4, 3, 2))
    plot(tab$iter, log10(f_gap), type = "l", col = "blue", lwd = 2,
         xaxt = "n", xlab = "", ylab = "log10(f-f[min])")
    main_title <- if (is.null(main))
      paste0(res$method, " - ",
           ifelse(res$convergence, "Converged", "Did not converge"))
    else
      main

    title(main = main_title)

    legend("topright", legend = "Criterion Gap", bty = "n", cex = 0.85)
    grid(lty = 3, col = "gray90")

    # 2. Middle panel - structural metric
    par(mar = c(0, 4, 2, 2))
    plot(0:(length(trajectory_history) - 1), mid_metric,
         type = "l", col = "darkgreen", lwd = 1.5,
         xaxt = "n", xlab = "", ylab = mid_label)
    if (mid_metric[length(mid_metric)] > mid_metric[1]) {
        legend("bottomright", legend = mid_title, bty = "n", cex = 0.85)
    } else {
        legend("topright",    legend = mid_title, bty = "n", cex = 0.85)
    }
    grid(lty = 3, col = "gray90")

    # 3. Bottom panel - gradient norm
    par(mar = c(4, 4, 2, 2))
    plot(tab$iter, tab$"log10.s.", type = "l", col = "darkred",
         xlab = "Iteration", ylab = "log10(s)")
    legend("topright", legend = "Gradient Norm", bty = "n", cex = 0.85)
    abline(h = log10(1e-5), lty = 2, col = "gray")
    grid(lty = 3, col = "gray90")
}



##########################################################
## landscape trajectory plot for 2 factor rotations only
##########################################################
     

plot_algorithm_comparison <- function(A, method = "varimax",
                                      factors = c(1, 2),
                                      maxit = 500, n = 20000,
                                      items = NULL, ...) {

  # Run all three algorithms
  rot_fn  <- get(method, envir = asNamespace("GPArotation"))
  res.leg <- rot_fn(A, algorithm = "legacy",  fwindow =  1,
                    maxit = maxit, ...)
  res.bb  <- rot_fn(A, algorithm = "bb",      fwindow = 10,
                    maxit = maxit, ...)
  res.cay <- rot_fn(A, algorithm = "cayley",  fwindow = 10,
                    maxit = maxit, ...)

  # --- Checks ---
  if (!res.leg$orthogonal)
    stop("plot_algorithm_comparison requires an orthogonal rotation method.\n",
         "Use one of: varimax, quartimax, bentlerT, geominT, parsimax, ",
         "simplimaxT, tandemI, tandemII, targetT, entropy, oblimax.",
         call. = FALSE)

  if (ncol(A) != 2)
    stop("plot_algorithm_comparison requires a 2-factor loading matrix.\n",
         "The landscape plot is only defined for 2-factor orthogonal rotations.",
         call. = FALSE)

  old_par <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old_par))

  # --- Top row: landscape trajectories ---
  suppressWarnings(
    plot_landscape_trajectory(res.leg, n = n, show_main = FALSE, show_sub = FALSE))
#  title(main = paste0("Legacy  |  Iters: ", nrow(res.leg$Table) - 1), cex.main = 0.9)

  suppressWarnings(
    plot_landscape_trajectory(res.bb, n = n, show_main = FALSE, show_sub = FALSE))
#  title(main = paste0("BB fw=10  |  Iters: ", nrow(res.bb$Table) - 1), cex.main = 0.9)

  suppressWarnings(
    plot_landscape_trajectory(res.cay, n = n, show_main = FALSE, show_sub = TRUE))
#  title(main = paste0("Cayley fw=10  |  Iters: ", nrow(res.cay$Table) - 1), cex.main = 0.9)

  # --- Overall title ---
  mtext(paste0("Algorithm comparison: ", method),
        outer = TRUE, cex = 1.1, font = 2, line = 0.5)

  # --- Bottom row: factor migration trajectories ---
  plot(res.leg, "trajectory", factors = factors, items = items)
  plot(res.bb,  "trajectory", factors = factors, items = items)
  plot(res.cay, "trajectory", factors = factors, items = items)

  invisible(list(legacy = res.leg, bb = res.bb, cayley = res.cay))
}


##########################################################################
## 1. Internal Utility: Trajectory Reconstruction Engine
##########################################################################
.reconstruct_trajectory <- function(res) {
  res <- .sortGPALoadings(res)
  
  A_unrotated <- attr(res, "A_unrotated")
  Tmat        <- res$Tinit
  method      <- sub("vgQ.", "", res$methodName)
  vgQfun      <- paste("vgQ", method, sep = ".")
  max_iter    <- nrow(res$Table) - 1
  fwindow     <- res$fwindow
  algorithm   <- res$algorithm
  methodArgs  <- res$methodArgs
  
  trajectory_history <- list()
  trajectory_history[[1]] <- Tmat
  
  f_history <- numeric(max_iter + 1)
  
  L      <- if(res$orthogonal) A_unrotated %*% Tmat else A_unrotated %*% t(solve(Tmat))
  VgQ    <- do.call(vgQfun, append(list(L), methodArgs))
  f      <- VgQ$f
  f_history[1] <- f
  G      <- if(res$orthogonal) crossprod(A_unrotated, VgQ$Gq) else -t(t(L) %*% VgQ$Gq %*% solve(Tmat))
  
  Tmat_prev <- NULL
  G_prev    <- NULL
  alpha     <- 1

  for (iter in 0:max_iter) {
    if (res$orthogonal) {
      M <- crossprod(Tmat, G); S <- (M + t(M)) / 2
      Gp <- G - Tmat %*% S
    } else {
      Gp <- G - Tmat %*% diag(colSums(Tmat * G))
    }
    s <- sqrt(sum(Gp^2))
    
    if (iter >= max_iter) break

    if (algorithm %in% c("bb", "cayley") && !is.null(Tmat_prev)) {
      dT <- Tmat - Tmat_prev; dGp <- Gp - G_prev
      bb_denom <- sum(dGp^2)
      if (bb_denom > 0) {
        alpha <- sum(dT^2) / abs(sum(dT * dGp))
        alpha <- max(1e-10, min(alpha, 10))
      }
    } else {
      alpha <- 2 * alpha
    }

    f_win_idx <- max(1, (iter + 1) - fwindow + 1):(iter + 1)
    target_f  <- max(f_history[f_win_idx], na.rm = TRUE)

    for (i in 0:10) {
      if (res$orthogonal) {
        if (algorithm == "cayley") {
          W <- Gp %*% t(Tmat) - Tmat %*% t(Gp); I <- diag(ncol(Tmat))
          Tmatt <- solve(I + (alpha/2) * W) %*% (I - (alpha/2) * W) %*% Tmat
        } else {
          X <- Tmat - alpha * Gp; UDV <- svd(X); Tmatt <- UDV$u %*% t(UDV$v)
        }
        L_test <- A_unrotated %*% Tmatt
      } else {
        X <- Tmat - alpha * Gp; v <- 1 / sqrt(colSums(X^2))
        Tmatt <- X %*% diag(v); L_test <- A_unrotated %*% t(solve(Tmatt))
      }
      
      VgQt <- do.call(vgQfun, append(list(L_test), methodArgs))
      if ((target_f - VgQt$f) > 0.5 * s^2 * alpha) break
      alpha <- alpha / 2
    }

    Tmat_prev <- Tmat; G_prev <- Gp
    Tmat <- Tmatt; f <- VgQt$f
    G <- if(res$orthogonal) crossprod(A_unrotated, VgQt$Gq) else -t(t(L_test) %*% VgQt$Gq %*% solve(Tmat))
    
    f_history[iter + 2] <- f
    trajectory_history[[iter + 2]] <- Tmat
  }
  
  return(trajectory_history)
}

##########################################################
## 2. Global Surface Visualizer Engine
##########################################################
plotRotationLandscape <- function(A, method = "quartimax", n = 20000,
                                  main = NULL, lwd = 2, col = "black", 
                                  xlim = NULL, ylim = NULL, ...) {
  if (ncol(A) != 2)
    stop("plotRotationLandscape only works for 2-factor solutions.")

  vgQfun_fn <- get(paste("vgQ", method, sep = "."), envir = asNamespace("GPArotation"))

  if (is.null(main))
    main <- paste("Rotation landscape:", method)

  theta  <- seq(0, 2 * pi, length.out = n)
  f_vals <- numeric(n)

  for (i in seq_along(theta)) {
    Tmat      <- matrix(c( cos(theta[i]), sin(theta[i]),
                          -sin(theta[i]), cos(theta[i])), 2, 2)
    L         <- A %*% Tmat
    VgQ       <- do.call(vgQfun_fn, append(list(L), list(...)))
    f_vals[i] <- VgQ$f
  }

  min_idx <- which.min(f_vals)
  is_inset <- !is.null(xlim)

  plot(theta, f_vals,
       type  = "l", lwd   = lwd, col   = col, main  = main,
       xlab  = if(is_inset) "" else expression(theta ~ "(radians)"),
       ylab  = if(is_inset) "" else "f",
       xaxt  = "n", yaxt  = if(is_inset) "n" else "s",
       xlim  = xlim, ylim  = ylim,
       panel.first = if(is_inset) rect(-1000, -1000, 1000, 1000, col="white", border=NA) else NULL)

  if (!is_inset) {
    axis(1, at     = c(0, pi/2, pi, 3*pi/2, 2*pi),
            labels = c("0", expression(pi/2), expression(pi),
                       expression(3*pi/2), expression(2*pi)))
                       
    text(theta[min_idx], f_vals[min_idx],
         labels = paste0("min f = ", round(f_vals[min_idx], 4)),
         pos    = 4, cex = 0.8)
  }

  abline(h = f_vals[min_idx], col = "grey80", lty = 2)
  points(theta[min_idx], f_vals[min_idx], col = "darkorange", pch = 17, cex = 1.5)

  invisible(data.frame(theta = theta, f = f_vals))
}

##########################################################
## 3. Macro/Micro Optimization Trajectory Path Overlay
##########################################################

plot_landscape_trajectory <- function(x, n = 20000, show_main = TRUE, 
                                      show_legend = TRUE, xlim = NULL, ylim = NULL, ...) {

  if (!inherits(x, "GPArotation"))
    stop("'x' must be a GPArotation object.", call. = FALSE)
  if (ncol(x$loadings) != 2)
    stop("plot_landscape_trajectory only works for 2-factor solutions.", call. = FALSE)
  if (!x$orthogonal)
    stop("plot_landscape_trajectory only works for orthogonal rotations.", call. = FALSE)
    
  # Silently reconstruct the path vectors using our engine
  traj       <- .reconstruct_trajectory(x)
  theta_path <- sapply(traj, function(Th) atan2(Th[2, 1], Th[1, 1]) %% (2 * pi))
  f_path     <- x$Table[1:length(traj), "f"]
  
  A_unrot    <- attr(x, "A_unrotated")
  if (is.null(A_unrot)) {
    A_unrot  <- x$loadings %*% solve(t(traj[[length(traj)]]))
  }
  
  method     <- sub("vgQ.", "", x$methodName)
  methodArgs <- x$methodArgs
  main_title <- if (show_main) NULL else ""
  
  do.call(plotRotationLandscape,
            c(list(A = A_unrot, method = method, n = n, 
            lwd = 0.8, col = "grey70", main = main_title,
            xlim = xlim, ylim = ylim), methodArgs))
            
  dx   <- diff(theta_path)
  dy   <- diff(f_path)
  dist <- sqrt(dx^2 + dy^2)
  
  # Set threshold to 1e-3 to suppress zero-length sub-pixel arrow warnings
  keep <- if (!is.null(xlim)) dist > 1e-3 else dist > 1e-2
  
  if (any(keep)) {
    suppressWarnings(
      arrows(theta_path[-length(theta_path)][keep],
             f_path[-length(f_path)][keep],
             theta_path[-1][keep],
             f_path[-1][keep],
             length = 0.05, col = "steelblue", lwd = 1.2)
      )
  }
  
  points(theta_path, f_path, col = "steelblue", pch = 19, cex = 0.6)
  points(theta_path[1], f_path[1], col = "darkgreen", pch = 19, cex = 1.3)
  points(theta_path[length(theta_path)], f_path[length(f_path)], col = "tomato", pch = 19, cex = 1.3)
    
  # Only draw the structural metadata text if this is the large macro landscape
#  if (is.null(xlim)) {
#    title(sub = paste0(x$algorithm, "  |  fwindow: ", x$fwindow,
#              "  |  Iters: ", nrow(x$Table) - 1), cex.sub = 0.8)
#  }
    
  if (show_legend) {
    legend("topright", bty = "n", cex = 0.75,
           legend  = c("Start", "Path", "End", "Global min"),
           col     = c("darkgreen", "steelblue", "tomato", "darkorange"),
           pch     = c(19, NA, 19, 17), 
           lty     = c(NA, 1, NA, NA),  
           lwd     = c(NA, 1.5, NA, NA))
  }
        
  invisible(list(theta = theta_path, f = f_path))
}

##########################################################################
## 4. Main Exported Function: Multi-Engine Benchmark Interface
##########################################################################
plot2fOrthComparison <- function(A, method = "quartimax", items = NULL,
                                 nang = 6000, maxit = 500, magnify = FALSE, 
                                 inset_pos = "auto", randomStart = FALSE, ...) {

  n <- nang
  if (!inherits(A, c("matrix", "array", "loadings", "factanal")))
    stop("'A' must be a matrix, loadings, or factanal object.", call. = FALSE)

  A_check <- if (inherits(A, "factanal")) unclass(A$loadings) else unclass(A)
  if (ncol(A_check) != 2)
    stop("This visualization strictly requires a 2-factor loading matrix.\n",
         "The landscape plot is only defined for 2-factor orthogonal rotations.", call. = FALSE)

  if (randomStart){
    if (inherits(A, "factanal")){
       A$loadings <- A$loadings %*% Random.Start(2)
    }
    else
    A <- A %*% Random.Start(2)
  }

  rot_fn  <- get(method, envir = asNamespace("GPArotation"))
  res.leg <- rot_fn(A, algorithm = "legacy", fwindow =  1, maxit = maxit, ...)
  res.bb  <- rot_fn(A, algorithm = "bb",     fwindow = 10, maxit = maxit, ...)
  res.cay <- rot_fn(A, algorithm = "cayley", fwindow = 10, maxit = maxit, ...)

  if (!res.leg$orthogonal)
    stop("This visualization requires an orthogonal rotation method.\n")

  # Establish a solid 2x3 column-by-column viewport matrix mapping
  layout_matrix <- matrix(1:6, nrow = 2, ncol = 3, byrow = FALSE)
  old_layout    <- layout(layout_matrix, heights = c(1.35, 1))
  
  old_par <- par(mar = c(3.8, 4.2, 2.5, 1.2), oma = c(1, 0, 3, 0), cex.axis = 0.9)
  on.exit({
    par(old_par)
    layout(1)
  })
  
  clean_method <- paste0(toupper(substring(method, 1, 1)), substring(method, 2))
  
  # Helper: Build bounds centered around the endgame destination coordinate
  calc_zoom <- function(traj) {
    final_t <- traj$theta[length(traj$theta)]
    final_f <- traj$f[length(traj$f)]
    n_tail  <- min(8, length(traj$f)) 
    n       <- length(traj$theta)
    
    x_span  <- min(max(abs(traj$theta[(n - n_tail + 1):n] - final_t)), 0.08)
    zoom_x  <- c(final_t - max(x_span * 1.5, 0.02), final_t + max(x_span * 1.5, 0.02))
    
    y_span  <- min(max(abs(traj$f[(n - n_tail + 1):n] - final_f)), 0.005)
    y_span  <- min(y_span, 0.005) 
    zoom_y  <- c(final_f - max(y_span/5, 1e-4), final_f + max(y_span * 1.5, 1e-3))
    
    list(x = zoom_x, y = zoom_y)
  }
  
  # Helper: Draw adaptive collision-avoiding overlay
  draw_inset <- function(res_obj, traj_obj) {
    withCallingHandlers({
      old_plt <- par("plt")
      u <- par("usr")
      w_usr <- u[2] - u[1]
      h_usr <- u[4] - u[3]
      
      # AUTOMATIC HEURISTIC: Kick box to upper left if trajectory crosses bottom-left
      current_pos <- inset_pos
      if (current_pos == "auto") {
        in_left_zone   <- (traj_obj$theta < (u[1] + w_usr * 0.25))
        in_bottom_zone <- (traj_obj$f < (u[3] + h_usr * 0.35))
        if (sum(in_left_zone & in_bottom_zone, na.rm = TRUE) > 1) {
          current_pos <- "topleft"
        } else {
          current_pos <- "bottomleft"
        }
      }
      
      # Map positions cleanly inside our panel slot frame
      inset_x1 <- old_plt[1] + (old_plt[2] - old_plt[1]) * 0.03
      inset_x2 <- old_plt[1] + (old_plt[2] - old_plt[1]) * 0.33
      
      if (current_pos == "topleft") {
        inset_y1 <- old_plt[3] + (old_plt[4] - old_plt[3]) * 0.55
        inset_y2 <- old_plt[3] + (old_plt[4] - old_plt[3]) * 0.90
      } else {
        inset_y1 <- old_plt[3] + (old_plt[4] - old_plt[3]) * 0.05
        inset_y2 <- old_plt[3] + (old_plt[4] - old_plt[3]) * 0.40
      }
      
      par(plt = c(inset_x1, inset_x2, inset_y1, inset_y2), new = TRUE)
      z_box <- calc_zoom(traj_obj)
      plot_landscape_trajectory(res_obj, n = n, show_main = FALSE,
                                show_legend = FALSE,
                                xlim = z_box$x, ylim = z_box$y)
      box(col = "gray40", lwd = 1.5)
      mtext("Endgame", side = 3, line = 0.1, cex = 0.7, font = 2)
      par(plt = old_plt)
    }, warning = function(w) {
      if (grepl("zero-length arrow", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  }

  # -------------------------------------------------------------------------
  # COLUMN 1: LEGACY
  # -------------------------------------------------------------------------
  traj_leg <- suppressWarnings(plot_landscape_trajectory(res.leg, n = n, show_main = FALSE))
  title(main = "Legacy, fwindow: 1", line = 0.8, cex.main = 1.0)
  if (magnify) suppressWarnings(draw_inset(res.leg, traj_leg))

  plot(res.leg, "trajectory", items = items, main = "")
  #mtext("Factor 1", side = 1, adj = 0.5, line = 1.9, cex = 0.65, outer = FALSE)
  # title(main = paste0("2a: Factor Migration Path"), line = 0.5, cex.main = 0.9, font.main = 3)

  # -------------------------------------------------------------------------
  # COLUMN 2: BARZILAI-BORWEIN (BB)
  # -------------------------------------------------------------------------
  traj_bb <- suppressWarnings(plot_landscape_trajectory(res.bb, n = n, show_main = FALSE))
  title(main = "Barzilai-Borwein, fwindow = 10", line = 0.8, cex.main = 1.0)
  if (magnify) suppressWarnings(draw_inset(res.bb,  traj_bb))
  
  plot(res.bb, "trajectory", items = items, main = "")
  #  mtext("Factor 1", side = 1, adj = 0.5, line = 1.9, cex =0.65,outer = FALSE)

  # title(main = paste0("2b: Factor Migration Path"), line = 0.5, cex.main = 0.9, font.main = 3)
  
  # -------------------------------------------------------------------------
  # COLUMN 3: CAYLEY
  # -------------------------------------------------------------------------
  traj_cay <- suppressWarnings(plot_landscape_trajectory(res.cay, n = n, show_main = FALSE))
  title(main = "Cayley, fwindow = 10", line = 0.8, cex.main = 1.0)
  if (magnify) suppressWarnings(draw_inset(res.cay, traj_cay))
  
  plot(res.cay, "trajectory", items = items, main = "")
  #  mtext("Factor 1", side = 1, adj = 0.5, line = 1.9, cex=0.65, outer = FALSE)

  # title(main = paste0("2c: Factor Migration Path"), line = 0.5, cex.main = 0.9, font.main = 3)
  
  # -------------------------------------------------------------------------
  # CANVAS HEADER BANNER
  # -------------------------------------------------------------------------
  #mtext(text = paste0("Comparison of Numerical Optimization Landscapes (", clean_method, " Rotation)"),
  #      side = 3, outer = TRUE, line = 0.8, cex = 1.2, font = 2)
  mtext(text = paste0("Rotation Landscapes and Migration Trajectories (", clean_method, " Rotation)"),
        side = 3, outer = TRUE, line = 0.8, cex = 1.2, font = 2)
  
  invisible(list(legacy = res.leg, bb = res.bb, cayley = res.cay))
}