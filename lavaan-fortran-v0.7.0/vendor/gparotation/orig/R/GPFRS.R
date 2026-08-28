# Main wrapper functions for GPForth and GPFoblq with random start functionality.
# Overhauled as part of 2026.4-1

.GPA_RS_engine <- function(A, Tmat, normalize = FALSE, eps = 1e-5, maxit = 2000,
                            method = NULL, methodArgs = NULL, randomStarts = 0,
                            orthogonal = TRUE, L = NULL, algorithm	= "bb", 
                            fwindow = 10, call = NULL, ...) {
  # Internal engine implementing random-start gradient projection rotation.
  #
  # Args:
  #   A            : Unrotated factor loading matrix (p x k)
  #   Tmat         : Initial rotation/transformation matrix
  #   normalize    : Kaiser normalization - logical or weight matrix
  #   eps          : Convergence tolerance (default 1e-5)
  #   maxit        : Maximum iterations (default 2000)
  #   method       : Rotation criterion name
  #   methodArgs   : Additional arguments passed to the vgQ criterion function
  #   randomStarts : Number of random starts (0 = use Tmat as-is)
  #   orthogonal   : TRUE for orthogonal (GPForth), FALSE for oblique (GPFoblq)
  #   L            : Deprecated alias for A
  #   algorithm	   : Algorithm used for step size alpha
  #   fwindow      : fwindow length size
  #   call         : call from named wrapper function
  #
  # Returns:
  #   A GPArotation object from the best start, with optional randStartChar diagnostics

  # --- Deprecation ---
  if (!is.null(L)) {
    warning("Argument 'L' is deprecated; using 'A' instead.", call. = FALSE)
    A <- L
  }

# --- Accept factanal object ---
  extracted   <- .extract_A(A)
  A           <- extracted$A
  n.obs       <- extracted$n.obs
  correlation <- extracted$correlation

  # --- Default Tmat now that A is a plain matrix ---
  if (is.null(Tmat))
    Tmat <- diag(ncol(A))
    
  # --- Validation set fwindow length---
  if (ncol(A) <= 1)
    stop("Rotation does not make sense for single-factor models.")

  valid_algos <- c("legacy", "bb", "cayley")
  algorithm <- tolower(algorithm)
  
  if (!(algorithm %in% valid_algos)){
    stop(paste("Invalid algorithm choice. Please select from:", 
    paste(valid_algos, collapse = ",")))
  }
  
  if (algorithm == "cayley" && orthogonal == FALSE){
    stop("cayley implented only in orthogonal rotation.")
  }
  
  if (is.null(fwindow)){
    fwindow <- switch(algorithm,
  	                   "legacy" = 1,
  	                   "bb"     = 10,
  	                   "cayley" = 10)
  }
  fwindow <- max(1, round(fwindow))
  
  # --- Setup ---
  engine        <- if (orthogonal) GPForth else GPFoblq
  actual_starts <- max(1, randomStarts)
  best_res      <- NULL
  Qvalues       <- numeric(actual_starts)
  Qconverged    <- logical(actual_starts)

  # --- Execution loop ---
  for (i in 1:actual_starts) {
    current_Tmat <- if (randomStarts > 0) Random.Start(ncol(A)) else Tmat

    res <- engine(A, Tmat = current_Tmat, normalize = normalize, eps = eps,
                  maxit = maxit, method = method, methodArgs = methodArgs, 
                  algorithm = algorithm, fwindow = fwindow)

    # Use indexing instead of tail() to avoid namespace warnings
    Qvalues[i]    <- res$Table[nrow(res$Table), 2]
    Qconverged[i] <- res$convergence

    if (is.null(best_res) || Qvalues[i] < best_res$Table[nrow(best_res$Table), 2])
      best_res <- res
  }

  # --- Random start diagnostics ---
  if (randomStarts > 1) {
    Q_round    <- round(Qvalues / eps) * eps
    Qmin_round <- min(Q_round)

    best_res$randStartChar <- c(
      randomStarts = randomStarts,
      Converged    = sum(Qconverged),
      atMinimum    = sum(Q_round == Qmin_round),
      localMins    = length(unique(Q_round))
    )
  }

  # --- Finalize for S3 print/summary ---
  if (!orthogonal)
    dimnames(best_res$Phi) <- list(colnames(A), colnames(A))

  # --- Attach factanal metadata from engine-level extraction ---
  best_res$n.obs       <- n.obs
  best_res$correlation <- correlation
  best_res$call        <- call
  best_res$methodName  <- method
  best_res$methodArgs  <- methodArgs
  best_res
}

.extract_A <- function(A) {
  # Accept a factanal object or raw matrix/loadings.
  # Returns a named list: A (plain matrix), n.obs, correlation.
  if (inherits(A, "factanal")) {
    list(
      A           = unclass(A$loadings),  # strip "loadings" class -> plain matrix
      n.obs       = A$n.obs,
      correlation = A$correlation
    )
  } else {
    list(
      A           = unclass(A),           # also strips "loadings" if passed directly
      n.obs       = NULL,
      correlation = NULL
    )
  }
}

GPFRSorth <- function(A, Tmat = NULL, normalize = FALSE, eps = 1e-5,
                      maxit = 2000, method = "varimax", methodArgs = NULL,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  # Orthogonal rotation with optional random starts. Wrapper for .GPA_RS_engine.
  #
  # Args:
  #   A            : Unrotated factor loading matrix (p x k)
  #   Tmat         : Initial rotation matrix (default: identity)
  #   normalize    : Kaiser normalization - logical or weight matrix
  #   eps          : Convergence tolerance (default 1e-5)
  #   maxit        : Maximum iterations (default 2000)
  #   method       : Rotation criterion (default "varimax")
  #   methodArgs   : Additional arguments passed to the vgQ criterion function
  #   randomStarts : Number of random starts (0 = use Tmat as-is)
  #   algorithm	   : Engine used for step size alpha
  #   fwindow      : fwindow length size
  
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, normalize = normalize, eps = eps,
                 maxit = maxit, method = method, methodArgs = methodArgs,
                 randomStarts = randomStarts, orthogonal = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

GPFRSoblq <- function(A, Tmat = NULL, normalize = FALSE, eps = 1e-5,
                      maxit = 2000, method = "quartimin", methodArgs = NULL,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  # Oblique rotation with optional random starts. Wrapper for .GPA_RS_engine.
  #
  # Args:
  #   A            : Unrotated factor loading matrix (p x k)
  #   Tmat         : Initial transformation matrix (default: identity)
  #   normalize    : Kaiser normalization - logical or weight matrix
  #   eps          : Convergence tolerance (default 1e-5)
  #   maxit        : Maximum iterations (default 2000)
  #   method       : Rotation criterion (default "quartimin")
  #   methodArgs   : Additional arguments passed to the vgQ criterion function
  #   randomStarts : Number of random starts (0 = use Tmat as-is)
  #   algorithm	   : Engine used for step size alpha
  #   fwindow      : fwindow length size


  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, normalize = normalize, eps = eps,
                 maxit = maxit, method = method, methodArgs = methodArgs,
                 randomStarts = randomStarts, orthogonal = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}
