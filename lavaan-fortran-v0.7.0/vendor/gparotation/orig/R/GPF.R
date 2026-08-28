GPForth <- function(A,
                    Tmat       = diag(ncol(A)),
                    normalize  = FALSE,
                    eps        = 1e-05,
                    maxit      = 2000,
                    method     = "varimax",
                    methodArgs = NULL,
                    algorithm  = "bb",
                    fwindow    = 10) {
  # Gradient Projection Algorithm for orthogonal rotation (Bernaards & Jennrich, 2005).
  #
  # Args:
  #   A          : Unrotated factor loading matrix (p x k)
  #   Tmat       : Initial rotation matrix (default: identity)
  #   normalize  : Kaiser normalization - logical or weight matrix
  #   eps        : Convergence tolerance on gradient norm (default 1e-5)
  #   maxit      : Maximum iterations (default 2000)
  #   method     : Rotation criterion, e.g. "varimax", "quartimin" (default "varimax")
  #   methodArgs : Additional arguments passed to the vgQ criterion function
  #   algorithm  : Algorithm used for step size alpha
  #   fwindow    : width of the non-monotone line search
  #
  # Returns:
  #   A GPArotation object with:
  #     loadings    - rotated loading matrix
  #     Th          - final rotation matrix
  #     Table       - iteration history (iter, f, log10(gradient norm), alpha)
  #     method      - method label from criterion function
  #     orthogonal  - TRUE (orthogonal rotation)
  #     convergence - TRUE if converged within maxit
  #     Gq          - final gradient matrix

  # --- Accept factanal object ---
  extracted   <- .extract_A(A)
  A           <- extracted$A
  n.obs       <- extracted$n.obs
  correlation <- extracted$correlation

  # --- Kaiser normalization ---
  if ((!is.logical(normalize)) || normalize) {
    W         <- NormalizingWeight(A, normalize = normalize)
    normalize <- TRUE
    A         <- A / W
  }

  if (ncol(A) < 2)
    stop("Rotation does not make sense for single-factor models.")

  A_unrotated <- A

  # --- validation ---
  algorithm <- tolower(algorithm)
  if (!(algorithm %in% c("legacy", "bb", "cayley"))){
    # fallback to package default without crashing
    algorithm <- "bb"
    }
  # ensure fwindow is an integer even if the user passed NULL directly
  if (is.null(fwindow)) fwindow <- ifelse(algorithm == "legacy", 1, 10)

  # --- Initialization ---
  Tinit  <- Tmat
  vgQfun <- paste("vgQ", method, sep = ".")
  alpha  <- 1
  Tmat_prev <- NULL
  Gp_prev    <- NULL
  L      <- A %*% Tmat
  VgQ    <- do.call(vgQfun, append(list(L), methodArgs))
  f      <- VgQ$f
  G      <- crossprod(A, VgQ$Gq)
  VgQt   <- do.call(vgQfun, append(list(L), methodArgs))
  Table  <- matrix(NA_real_, nrow = maxit + 1, ncol = 4,
                   dimnames = list(NULL, c("iter", "f", "log10(s)", "alpha")))

  # --- Main loop ---
  for (iter in 0:maxit) {

    M  <- crossprod(Tmat, G)
    S  <- (M + t(M)) / 2
    Gp <- G - Tmat %*% S
    s  <- sqrt(sum(Gp^2))

    Table[iter + 1, ] <- c(iter, f, log10(s), alpha)

    if (s < eps) break

    # --- Barzilai-Borwein step size ---
    if (algorithm %in% c("bb", "cayley") && !is.null(Tmat_prev)) {
      dT       <- Tmat - Tmat_prev
      dGp      <- Gp    - Gp_prev
      bb_denom <- sum(dGp^2)
      if (bb_denom > 0){
        # alpha <- abs(sum(dT * dGp) / bb_denom) # alternate BB
        alpha <- sum(dT^2) / abs(sum(dT * dGp))  # standard BB
        alpha <- max(1e-10, min(alpha, 20))}
        # BB sets alpha directly - enter line search without doubling
    } else {
      alpha <- 2 * alpha  # original behaviour
    }

    fwindow_range <- max(1, (iter + 1) - fwindow + 1):(iter + 1)
    target_f <- max(Table[fwindow_range ,2], na.rm = TRUE)

    # Line search
    for (i in 0:10) {
      if (algorithm == "cayley") {
        # --- CAYLEY (The Curve) ---
        W <- Gp %*% t(Tmat) - Tmat %*% t(Gp)
        I <- diag(ncol(Tmat))
        X_inv <- solve(I + (alpha / 2) * W)
        Tmatt <- X_inv %*% (I - (alpha / 2) * W) %*% Tmat
      } else {
        # --- SVD (The Snap-Back) ---
        X     <- Tmat - alpha * Gp
        UDV   <- svd(X)
        Tmatt <- UDV$u %*% t(UDV$v)
      }
      L     <- A %*% Tmatt
      VgQt  <- do.call(vgQfun, append(list(L), methodArgs))
      if ((target_f - VgQt$f) > 0.5 * s^2 * alpha) 
        break
      alpha <- alpha / 2
    }

    Tmat_prev <- Tmat
    Gp_prev    <- Gp
    Tmat <- Tmatt
    f    <- VgQt$f
    G    <- crossprod(A, VgQt$Gq)
  }
  
  # --- Convergence check ---
  convergence <- (s < eps)
  if (iter == maxit && !convergence) {
    conv_gap <- paste0("Algorithm stopped at s = ", signif(s, 4),
                       " (target eps = ", eps, ").")
    suggestion <- if (algorithm == "legacy")
      "Consider: increase maxit or try algorithm = \"bb\" with fwindow = 10."
    else if (algorithm == "bb")
      "Consider: increase maxit, increase fwindow, or try algorithm = \"cayley\"."
    else if (algorithm == "cayley")
      "Consider: increase maxit or try algorithm = \"bb\" with fwindow = 10."
    warning("Convergence not obtained in GPForth. ", maxit, " iterations used.\n",
            conv_gap, "\n", suggestion)
  }

  # --- Trim Table ---
  Table <- Table[1:(iter + 1), , drop = FALSE]

  # --- Undo normalization ---
  if (normalize)
    L <- L * W

  dimnames(L) <- dimnames(A)

  class(L) <- "loadings"
  attr(L, "covariance") <- diag(ncol(L))

  r <- list(
    loadings    = L,
    Th          = Tmat,
    Tinit       = Tinit,
    Table       = Table,
    method      = VgQ$Method,
    methodName  = method,
    methodArgs  = methodArgs,
    orthogonal  = TRUE,
    convergence = convergence,
    algorithm   = algorithm,
    fwindow     = fwindow,
    Gq          = VgQt$Gq,
    n.obs       = n.obs,
    correlation = correlation
  )
  attr(r, "A_unrotated") <- A_unrotated
  class(r) <- "GPArotation"
  r
}



GPFoblq <- function(A,
                    Tmat       = diag(ncol(A)),
                    normalize  = FALSE,
                    eps        = 1e-05,
                    maxit      = 2000,
                    method     = "quartimin",
                    methodArgs = NULL,
                    algorithm  = "bb",
                    fwindow    = 10) {
  # Gradient Projection Algorithm for oblique rotation (Bernaards & Jennrich, 2005).
  #
  # Args:
  #   A          : Unrotated factor loading matrix (p x k)
  #   Tmat       : Initial transformation matrix (default: identity)
  #   normalize  : Kaiser normalization - logical or weight matrix
  #   eps        : Convergence tolerance on gradient norm (default 1e-5)
  #   maxit      : Maximum iterations (default 2000)
  #   method     : Rotation criterion, e.g. "quartimin" (default "quartimin")
  #   methodArgs : Additional arguments passed to the vgQ criterion function
  #   algorithm  : Algorithm used for step size alpha
  #   fwindow    : width of the non-monotone line search
  #
  # Returns:
  #   A GPArotation object with:
  #     loadings    - rotated loading matrix
  #     Phi         - factor correlation matrix (t(Tmat) %*% Tmat)
  #     Th          - final transformation matrix
  #     Table       - iteration history (iter, f, log10(gradient norm), alpha)
  #     method      - method label from criterion function
  #     orthogonal  - FALSE (oblique rotation)
  #     convergence - TRUE if converged within maxit
  #     Gq          - final gradient matrix  
  #     n.obs       - n.obs from factanal;  NULL if not from factanal
  #     correlation - correlation matrix from factanal; NULL if not from factanal
  singularity_warning <- FALSE
  
  # --- Accept factanal object ---
  extracted   <- .extract_A(A)
  A           <- extracted$A
  n.obs       <- extracted$n.obs
  correlation <- extracted$correlation

  if (ncol(A) < 2)
    stop("Rotation does not make sense for single-factor models.")

  # --- Kaiser normalization ---
  if ((!is.logical(normalize)) || normalize) {
    W         <- NormalizingWeight(A, normalize = normalize)
    normalize <- TRUE
    A         <- A / W
  }

  A_unrotated <- A

  # --- Safe inverse ---
  safe_inverse <- function(X) {
    tryCatch(
      solve(X),
      error = function(e) {
        assign("singularity_warning", TRUE, envir = parent.frame(n = 2))
        s_svd <- svd(X)
        tol   <- sqrt(.Machine$double.eps)
        d_inv <- ifelse(s_svd$d > tol, 1 / s_svd$d, 0)
        s_svd$v %*% diag(d_inv) %*% t(s_svd$u)
      }
    )
  }
  
  # --- validation ---
  if (algorithm == "cayley")
    stop("algorithm = \"cayley\" is only available for orthogonal rotation. ",
         "Use GPForth or an orthogonal rotation wrapper.", call. = FALSE)
  algorithm <- tolower(algorithm)
  if (!(algorithm %in% c("legacy", "bb"))){
    # fallback to package default without crashing
    algorithm <- "bb"
    }
         
  # ensure fwindow is an integer even if the user passed NULL directly
  if (is.null(fwindow)) fwindow <- ifelse(algorithm == "legacy", 1, 10)  

  # --- Initialization ---
  Tinit    <- Tmat
  vgQfun   <- paste("vgQ", method, sep = ".")
  alpha    <- 1
  Tmat_prev <- NULL
  Gp_prev    <- NULL
  Tmat_inv <- safe_inverse(Tmat)
  L        <- A %*% t(Tmat_inv)
  VgQ      <- do.call(vgQfun, append(list(L), methodArgs))
  f        <- VgQ$f
  G        <- -t(t(L) %*% VgQ$Gq %*% Tmat_inv)
  VgQt     <- do.call(vgQfun, append(list(L), methodArgs))
  Table    <- matrix(NA_real_, nrow = maxit + 1, ncol = 4,
                     dimnames = list(NULL, c("iter", "f", "log10(s)", "alpha")))

  # --- Main loop ---
  for (iter in 0:maxit) {

    Gp <- G - Tmat %*% diag(colSums(Tmat * G))
    s  <- sqrt(sum(Gp^2))

    Table[iter + 1, ] <- c(iter, f, log10(s), alpha)

    if (s < eps) break

    # --- Barzilai-Borwein step size ---
    if (algorithm == "bb"  && !is.null(Tmat_prev)) {
      dT       <- Tmat - Tmat_prev
      dGp       <- Gp   - Gp_prev
      bb_denom <- sum(dGp^2)
      if (bb_denom > 0){
        # alpha <- abs(sum(dT * dGp) / bb_denom) # alternate BB
        alpha <- sum(dT^2) / abs(sum(dT * dGp))  # standard BB
        alpha <- max(1e-10, min(alpha, 20))}
        # BB sets alpha directly - enter line search without doubling
    } else {
      alpha <- 2 * alpha  # original behaviour
    }
   
    fwindow_range <- max(1, (iter + 1) - fwindow + 1):(iter + 1)
    target_f <- max(Table[fwindow_range ,2], na.rm = TRUE)

    # Line search
    for (i in 0:10) {
      X         <- Tmat - alpha * Gp
      v         <- 1 / sqrt(colSums(X^2))
      Tmatt     <- X %*% diag(v)
      Tmatt_inv <- safe_inverse(Tmatt)
      L         <- A %*% t(Tmatt_inv)
      VgQt      <- do.call(vgQfun, append(list(L), methodArgs))
      if ((target_f - VgQt$f) > 0.5 * s^2 * alpha) 
        break
      alpha     <- alpha / 2
    }

    Tmat_prev <- Tmat
    Gp_prev    <- Gp
    Tmat <- Tmatt
    f    <- VgQt$f
    G    <- -t(t(L) %*% VgQt$Gq %*% Tmatt_inv)
  }

  # --- Convergence check with smart diagnostics ---
  convergence <- (s < eps)
  if (iter == maxit && !convergence) {
    Phi     <- t(Tmat) %*% Tmat
    max_cor <- max(abs(Phi[lower.tri(Phi)]))
    conv_gap <- paste0("Algorithm stopped at s = ", signif(s, 4),
                       " (target eps = ", eps, ").")
    suggestion <- if (algorithm == "legacy")
      "Consider: increase maxit or try algorithm = \"bb\" with fwindow = 10."
    else if (algorithm == "bb")
      "Consider: increase maxit or increase fwindow."
    if (max_cor > 0.85) {
      warning("Convergence not obtained in GPFoblq. ", maxit, " iterations used.\n",
              conv_gap, "\n",
              "DIAGNOSTIC: Extreme factor correlations detected (Max = ",
              round(max_cor, 3), ").\n",
              "This almost always indicates factor over-extraction. ",
              "Try reducing the number of extracted factors.")
    } else {
      warning("Convergence not obtained in GPFoblq. ", maxit, " iterations used.\n",
              conv_gap, "\n", suggestion)
    }
  }

  # --- Trim Table ---
  Table <- Table[1:(iter + 1), , drop = FALSE]

  # --- Undo normalization ---
  if (normalize)
    L <- L * W

  dimnames(L) <- dimnames(A)

  # --- covariance attribute for print.loadings ---
  Phi <- t(Tmat) %*% Tmat
  class(L) <- "loadings"
  attr(L, "covariance") <- Phi

  r <- list(
    loadings    = L,
    Phi         = Phi,
    Th          = Tmat,
    Tinit       = Tinit,
    Table       = Table,
    method      = VgQ$Method,
    methodName  = method,
    methodArgs  = methodArgs,
    orthogonal  = FALSE,
    convergence = convergence,
    algorithm   = algorithm,
    fwindow     = fwindow,
    Gq          = VgQt$Gq,
    n.obs       = n.obs,
    correlation = correlation
  )
  attr(r, "A_unrotated") <- A_unrotated
  class(r) <- "GPArotation"
  r
}


###########################################################################
# Legacy functions - retained for historical reference and reproducibility.
# Original implementations from Bernaards & Jennrich (2005).
# Superseded by the updated GPForth and GPFoblq above (GPArotation 2026.4-1).
# Note: commented-out lines are preserved intentionally as part of the
# original development history.
###########################################################################


GPForth.legacy <- function(A,
                           Tmat       = diag(ncol(A)),
                           normalize  = FALSE,
                           eps        = 1e-05,
                           maxit      = 1000,
                           method     = "varimax",
                           methodArgs = NULL) {
  # Legacy implementation of GPForth (Bernaards & Jennrich, 2005).
  # Retained for historical reference and reproducibility.
  # Superseded by GPForth in GPArotation 2026.4-1.
  if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 al <- 1
 L <- A %*% Tmat
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- crossprod(A,VgQ$Gq)
 f <- VgQ$f
 Table <- NULL
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   M <- crossprod(Tmat,G)
   S <- (M + t(M))/2
   Gp <- G - Tmat %*% S
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table, c(iter, f, log10(s), al))
   if (s < eps)  break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al * Gp
     UDV <- svd(X)
     Tmatt <- UDV$u %*% t(UDV$v)
     L <- A %*% Tmatt
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     if (VgQt$f < (f - 0.5*s^2*al)) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- crossprod(A,VgQt$Gq)
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPForth. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)
 r <- list(loadings=L, Th=Tmat, Table=Table, 
        method=VgQ$Method, orthogonal=TRUE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}

GPFoblq.legacy <- function(A,
                           Tmat       = diag(ncol(A)),
                           normalize  = FALSE,
                           eps        = 1e-05,
                           maxit      = 1000,
                           method     = "quartimin",
                           methodArgs = NULL) {
  # Legacy implementation of GPFoblq (Bernaards & Jennrich, 2005).
  # Retained for historical reference and reproducibility.
  # Superseded by GPFoblq in GPArotation 2026.4-1.
  
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 al <- 1
 L <- A %*% t(solve(Tmat))
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- -t(t(L) %*% VgQ$Gq %*% solve(Tmat))
 f <- VgQ$f
 Table <- NULL
 #Table <- c(-1,f,log10(sqrt(sum(diag(crossprod(G))))),al)
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   Gp <- G - Tmat %*% diag(c(rep(1,nrow(G)) %*% (Tmat*G)))
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table,c(iter,f,log10(s),al))
   if (s < eps) break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al*Gp
     v <- 1/sqrt(c(rep(1,nrow(X)) %*% X^2))
     Tmatt <- X %*% diag(v)
     L <- A %*% t(solve(Tmatt))
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     improvement <- f - VgQt$f 
     if (improvement >  0.5*s^2*al) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- -t(t(L) %*% VgQt$Gq %*% solve(Tmatt))
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPFoblq. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)

 # N.B. renaming Lh to loadings in specificific rotations 
 #   uses fact that  Lh is first.
 r <- list(loadings=L, Phi=t(Tmat) %*% Tmat, Th=Tmat, Table=Table,
      method=VgQ$Method, orthogonal=FALSE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}