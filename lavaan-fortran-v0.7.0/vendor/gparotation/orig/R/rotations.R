###########################################################################
###########################################################################

# --- OBLIMIN ---

oblimin <- function(A, Tmat = NULL, gam = 0, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  # Oblimin oblique rotation. Special cases via gam:
  #   0   - Quartimin
  #   0.5 - Biquartimin
  #
  # Args:
  #   A            : Unrotated factor loading matrix (p x k)
  #   Tmat         : Initial transformation matrix (default: identity)
  #   gam          : Obliqueness parameter (default 0)
  #   normalize    : Kaiser normalization - logical or weight matrix
  #   randomStarts : Number of random starts (0 = use Tmat as-is)

  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "oblimin",
                 methodArgs   = list(gam = gam),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc,
                 algorithm    = algorithm, fwindow = fwindow, ...)
}

vgQ.oblimin <- function(L, gam = 0) {
  # Oblimin family of oblique rotation criteria.
  # gam controls the degree of obliqueness:
  #   0    - Quartimin
  #   0.5  - Biquartimin
  #
  # Args:
  #   L   : Factor loading matrix
  #   gam : Obliqueness parameter (default 0)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label
  X <- rowSums(L^2) - L^2
  if (gam != 0) {
        X <- sweep(X, 2, (gam / nrow(L)) * colSums(X), "-")
  }
  method <- if (gam == 0) "Oblimin Quartimin"
              else if (gam == 0.5) "Oblimin Biquartimin"
              else paste0("Oblimin gamma =", gam)
    list(Gq = L * X, f = sum(L^2 * X)/4, Method = method)
}

# --- QUARTIMIN ---

quartimin <- function(A, Tmat = NULL, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "quartimin",
                 normalize    = normalize, randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, algorithm = algorithm, 
                 fwindow = fwindow, ...)
}

vgQ.quartimin <- function(L) {
  # Quartimin oblique rotation criterion (special case of Oblimin with gam = 0).
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  k  <- ncol(L)
  L2 <- L^2
  X  <- L2 %*% (matrix(1, k, k) - diag(k))

  list(
    Gq     = L * X,
    f      = sum(L2 * X) / 4,
    Method = "Quartimin"
  )
}

# --- CRAWFORD-FERGUSON FAMILY ---

cfT <- function(A, Tmat = NULL, kappa = 0, normalize = FALSE,
                randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "cf",
                 methodArgs   = list(kappa = kappa), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

cfQ <- function(A, Tmat = NULL, kappa = 0, normalize = FALSE,
                randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "cf",
                 methodArgs   = list(kappa = kappa), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow,  ...)
}

equamax <- function(A, Tmat = NULL, kappa = NULL,
                    normalize = FALSE, randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  # Logic for default kappa if not provided by user
  if (is.null(kappa)) {
    extracted <- .extract_A(A)
    # actual calculation using extracted matrix dims
    kappa <- ncol(extracted$A) / (2 * nrow(extracted$A))
  }

  .GPA_RS_engine(A = A, Tmat = Tmat, method = "cf",
                 methodArgs   = list(kappa = kappa), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow,  ...)
}

parsimax <- function(A, Tmat = NULL, kappa = NULL,
                      normalize = FALSE, randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  if (is.null(kappa)) {
    extracted <- .extract_A(A)
    p <- nrow(extracted$A)
    m <- ncol(extracted$A)
    kappa <- (m - 1) / (m + p - 2)
  }

  .GPA_RS_engine(A = A, Tmat = Tmat, method = "cf",
                 methodArgs   = list(kappa = kappa), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow,  ...)
}

vgQ.cf <- function(L, kappa = 0) {
  # Crawford-Ferguson (CF) family of rotation criteria (1970).
  # kappa controls the complexity penalty:
  #   0              - Quartimax / Quartimin
  #   1/p            - Varimax
  #   k/(2*p)        - Equamax
  #   (k-1)/(p+k-2)  - Parsimax
  #   1              - Factor Parsimony
  #
  # Args:
  #   L     : Factor loading matrix (p x k)
  #   kappa : Complexity weight (default 0)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  p <- nrow(L)
  k <- ncol(L)
  L2 <- L^2
   
  # rowSums(L2) - L2 is the same as L2 %*% N
  # colSums(L2) - L2 is the same as M %*% L2
  rowS <- rowSums(L2)
  colS <- colSums(L2)
   
  # Efficiently broadcast the sums without building p x p matrices
  L2N <- sweep(-L2, 1, rowS, "+")
  ML2 <- sweep(-L2, 2, colS, "+")
   
  f1 <- (1 - kappa) * sum(L2 * L2N) / 4
  f2 <- kappa * sum(L2 * ML2) / 4

  method <- if (kappa == 0) 
        "Crawford-Ferguson Quartimax/Quartimin"
  else if (kappa == 1/p) 
        "Crawford-Ferguson Varimax"
  else if (kappa == k/(2 * p)) 
        "Equamax"
  else if (kappa == (k - 1)/(p + k - 2)) 
        "Parsimax"
  else if (kappa == 1) 
        "Factor Parsimony"
  else paste0("Crawford-Ferguson: kappa=", kappa)
   
  list(Gq = (1 - kappa) * L * L2N + kappa * L * ML2,
         f = f1 + f2, Method = method)
}


# --- TARGET ROTATION ---

targetT <- function(A, Tmat = NULL, Target = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  if (is.null(Target)) stop("Argument 'Target' must be specified.")

  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "target",
                 methodArgs   = list(Target = Target), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow,  ...)
}

targetQ <- function(A, Tmat = NULL, Target = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  if (is.null(Target)) stop("Argument 'Target' must be specified.")

  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "target",
                 methodArgs   = list(Target = Target), normalize    = normalize,
                 randomStarts = randomStarts, orthogonal = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.target <- function(L, Target = NULL) {
  # Target rotation criterion.
  # Minimizes the sum of squared deviations from a specified target matrix.
  # NA entries in Target are treated as unspecified (ignored in fit and gradient).
  #
  # Args:
  #   L      : Factor loading matrix
  #   Target : Target matrix (same dimensions as L; NA = unspecified element)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (zero where Target is NA)
  #     f      - objective function value
  #     Method - method label

  Gq            <- 2 * (L - Target)
  Gq[is.na(Gq)] <- 0

  list(
    Gq     = Gq,
    f      = sum((L - Target)^2, na.rm = TRUE),
    Method = "Target rotation"
  )
}

# --- PARTIALLY SPECIFIED TARGET ---

pstT <- function(A, Tmat = NULL, W = NULL, Target = NULL, normalize = FALSE,
                 randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  if (is.null(W))      stop("Argument 'W' must be specified.")
  if (is.null(Target)) stop("Argument 'Target' must be specified.")

  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "pst",
                 methodArgs   = list(W = W, Target = Target),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

pstQ <- function(A, Tmat = NULL, W = NULL, Target = NULL, normalize = FALSE,
                 randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  if (is.null(W))      stop("Argument 'W' must be specified.")
  if (is.null(Target)) stop("Argument 'Target' must be specified.")

  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "pst",
                 methodArgs   = list(W = W, Target = Target),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.pst <- function(L, W = NULL, Target = NULL) {
  # Partially specified target rotation criterion.
  # Minimizes weighted deviations from a target matrix.
  # W is a binary weight matrix specifying which elements are targeted.
  #
  # Args:
  #   L      : Factor loading matrix
  #   W      : Binary weight matrix (same dimensions as L; 0 = unspecified element)
  #   Target : Target matrix (same dimensions as L)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  Btilde <- W * Target

  list(
    Gq     = 2 * (W * L - Btilde),
    f      = sum((W * L - Btilde)^2),
    Method = "Partially specified target"
  )
}

# --- ENTROPY ---

entropy <- function(A, Tmat = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "entropy",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.entropy <- function(L) {
  # Minimum entropy oblique rotation criterion.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  L2    <- L^2
  logL2 <- log(L2 + .Machine$double.eps)

  list(
    Gq     = -(L * logL2 + L),
    f      = -sum(L2 * logL2) / 2,
    Method = "Minimum entropy"
  )
}

# --- INFOMAX ---

infomaxT <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "infomax",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

infomaxQ <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "infomax",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.infomax <- function(L) {
  # Infomax oblique rotation criterion.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  k <- ncol(L)
  p <- nrow(L)
  S <- L^2
  s <- sum(S)
  eps <- .Machine$double.eps
    
  # Proportions
  E <- S/s
  e1 <- rowSums(S)/s
  e2 <- colSums(S)/s
    
  # Entropies
  Q0 <- -sum(E * log(E + eps))
  Q1 <- -sum(e1 * log(e1 + eps))
  Q2 <- -sum(e2 * log(e2 + eps))
    
  # Gradients of Entropies
  H <- -(log(E + eps) + 1)
  h1 <- -(log(e1 + eps) + 1)
  h2 <- -(log(e2 + eps) + 1)
    
  # Centered components (E[H] logic)
  G0 <- (H - sum(S * H)/s) / s
  G1 <- (matrix(h1, p, k) - sum(rowSums(S) * h1)/s) / s
  G2 <- (matrix(h2, p, k, byrow = TRUE) - sum(colSums(S) * h2)/s) / s
    
  list(Gq = 2 * L * (G0 - G1 - G2), f = log(k) + Q0 - Q1 - Q2, Method = "Infomax")
}


# --- McCAMMON ---

mccammon <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "mccammon",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

# This vgQ updated on May 6, 2026 to align with lavaan
# Checked against numerical gradient. Only Gq needed update
# 
vgQ.mccammon <- function(L) {
  # Basic dimensions
  p <- nrow(L)
  k <- ncol(L)
 
  # Square loadings and sums
  S <- L^2
  col_sums <- colSums(S)
  total_s <- sum(S)
 
  # 1. Normalize for Q1 (Column entropy)
  # P is the proportion of a column's variance contributed by each row
  P <- sweep(S, 2, col_sums, "/")
 
  # Safe log function to handle zero loadings
  safe_log <- function(x) {
    res <- x
    res[x > 0] <- log(x[x > 0])
    res[x <= 0] <- 0
    return(res)
  }
 
  logP <- safe_log(P)
 
  # Q1: Entropy within columns (encourages simple structure per factor)
  Q1 <- -sum(P * logP)
 
  # 2. Normalize for Q2 (Across-column entropy)
  # p2 is the proportion of total variance in each column
  p2 <- col_sums / total_s
  logp2 <- safe_log(p2)
 
  # Q2: Entropy of column sums (encourages equal-sized factors)
  Q2 <- -sum(p2 * logp2)
 
  # --- Analytical Gradient Calculation ---
 
  # Component 1: Gradient of Q1 wrt S
  # H1 matches Bernaards & Jennrich logic: -(log(P) + 1)
  H1 <- -(logP + 1)
  # alpha1 is the column-wise expected value of H1
  alpha1 <- colSums(P * H1)
  # G1 is the centered gradient scaled by column sums
  G1 <- sweep(H1, 2, alpha1, "-")
  G1 <- sweep(G1, 2, col_sums, "/")
 
  # Component 2: Gradient of Q2 wrt S
  # h2 is the vector of gradients for each column sum
  h2 <- -(logp2 + 1)
  # alpha2 is the total expected value of h2
  alpha2 <- sum(p2 * h2)
  # G2 is the centered gradient scaled by the total sum of squares
  # We expand the vector h2 into a matrix that matches L
  G2_vec <- (h2 - alpha2) / total_s
  G2 <- matrix(G2_vec, p, k, byrow = TRUE)
 
  # Combine components for the total gradient Gq
  # f = log(Q1) - log(Q2)
  # d/dL f = 2 * L * ( (1/Q1)*dQ1/dS - (1/Q2)*dQ2/dS )
  Gq <- 2 * L * (G1 / Q1 - G2 / Q2)
 
  f <- log(Q1) - log(Q2)
 
  return(list(Gq = Gq, f = f, Method = "McCammon entropy"))
}



# --- GEOMIN ---

geominT <- function(A, Tmat = NULL, delta = 0.01, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "geomin",
                 methodArgs   = list(delta = delta),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

geominQ <- function(A, Tmat = NULL, delta = 0.01, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "geomin",
                 methodArgs   = list(delta = delta),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.geomin <- function(L, delta = 0.01) {
  # Geomin oblique rotation criterion.
  # delta is a small constant added for numerical stability (default 0.01).
  #
  # Args:
  #   L     : Factor loading matrix
  #   delta : Stabilizing constant (default 0.01)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  k   <- ncol(L)
  L2  <- L^2 + delta
  pro <- exp(rowMeans(log(L2)))

  list(
    Gq     = (2 / k) * (L / L2) * pro,
    f      = sum(pro),
    Method = "Geomin"
  )
}

# --- SIMPLIMAX ---

simplimax <- function(A, Tmat = NULL, k = NULL, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
    # 1. Capture the call
    mc <- match.call()
    
    # 2. Extract the matrix FIRST so we know the true nrow
    extracted <- .extract_A(A)
    actual_A  <- extracted$A

    # 3. If k wasn't provided, set it safely based on the extracted matrix
    if (is.null(k)) {
        k <- nrow(actual_A) * (ncol(actual_A) - 1)
    }

    # 4. Pass everything to the engine
    .GPA_RS_engine(A = A, Tmat = Tmat, method = "simplimax",
        methodArgs = list(k = k), normalize = normalize, randomStarts = randomStarts,
        orthogonal = FALSE, call = mc, algorithm = algorithm, fwindow = fwindow,
        ...)
}

vgQ.simplimax <- function(L, k = nrow(L) * (ncol(L)-1)) {
  # Simplimax oblique rotation criterion (Kiers, 1994).
  # Minimizes the k smallest squared loadings, targeting a simple structure
  # with k near-zero loadings per solution.
  #
  # Args:
  #   L : Factor loading matrix
  #   k : Number of near-zero loadings to target (default: nrow(L) * (ncol(L)-1))
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  L2 <- L^2
  # Create an empty logical matrix
  Imat <- matrix(FALSE, nrow(L), ncol(L))
    
  # Explicitly pick the k smallest indices to avoid threshold ties
  idx <- order(L2)[1:k]
  Imat[idx] <- TRUE
    
  list(Gq = 2 * Imat * L, f = sum(L2[Imat]), 
       Method = paste0("Simplimax (k=", k, ")"))
}

# --- BIFACTOR ---

bifactorT <- function(A, Tmat = NULL, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bifactor",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

bifactorQ <- function(A, Tmat = NULL, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bifactor",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.bifactor <- function(L) {
  # Bifactor Biquartimin oblique rotation criterion.
  # Treats the first column as a general factor (unconstrained)
  # and rotates the remaining group factor columns.
  #
  # Args:
  #   L : Factor loading matrix (first column = general factor)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (zero gradient for general factor column)
  #     f      - objective function value
  #     Method - method label

  Lt2  <- L[, -1, drop = FALSE]^2
  rowS <- rowSums(Lt2)
  Gt   <- 4 * L[, -1, drop = FALSE] * (rowS - Lt2)

  list(
    Gq     = cbind(0, Gt),
    f      = sum(rowS^2) - sum(Lt2^2),
    Method = "Bifactor Biquartimin"
  )
}

# --- BI-GEOMIN ---

bigeominT <- function(A, Tmat = NULL, delta = 0.01, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bigeomin",
                 methodArgs   = list(delta = delta),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

bigeominQ <- function(A, Tmat = NULL, delta = 0.01, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bigeomin",
                 methodArgs   = list(delta = delta),
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.bigeomin <- function(L, delta = 0.01) {
  # Bi-Geomin oblique rotation criterion.
  # Treats the first column as a general factor (unconstrained)
  # and applies Geomin rotation to the remaining group factor columns.
  #
  # Args:
  #   L     : Factor loading matrix (first column = general factor)
  #   delta : Stabilizing constant passed to vgQ.geomin (default 0.01)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (zero gradient for general factor column)
  #     f      - objective function value
  #     Method - method label

  out <- vgQ.geomin(L[, -1, drop = FALSE], delta = delta)

  list(
    Gq     = cbind(0, out$Gq),
    f      = out$f,
    Method = "Bi-Geomin"
  )
}

# --- TANDEM CRITERIA ---

tandemI <- function(A, Tmat = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "tandemI",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.tandemI <- function(L) {
  # Tandem Criterion I (Comrey, 1967)
  # Seeks factors that share high loadings on the same variables.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value (negated for minimization)
  #     Method - method label

  L2 <- L^2
  # tcrossprod(L) is L %*% t(L), but more numerically stable
  LL <- tcrossprod(L)
  LL2 <- LL^2
 
  # Pre-compute the matrix product used in both f and Gq1
  # This reduces redundant calculations and preserves precision
  M_prod <- LL2 %*% L2
 
  # f = -tr(L2' %*% LL2 %*% L2)
  # Replacing sum(diag(crossprod(...))) with sum(A * B)
  # avoids intermediate matrix noise
  f <- -sum(L2 * M_prod)
 
  # Gq1 = 4 * L * (LL2 %*% L2)
  Gq1 <- 4 * L * M_prod
 
  # Gq2 = 4 * (LL * (L2 %*% t(L2))) %*% L
  # We use tcrossprod here as well for symmetry and stability
  Gq2 <- 4 * (LL * tcrossprod(L2)) %*% L
 
  Gq <- -(Gq1 + Gq2)
 
  list(Gq = Gq, f = f, Method = "Tandem I")
}


tandemII <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "tandemII",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.tandemII <- function(L) {
  # Tandem Criterion II (Comrey, 1967)
  # Seeks factors that do NOT share high loadings on the same variables
  # (used as a second-stage refinement after Tandem I).
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  L2 <- L^2
  LL <- tcrossprod(L)
  LL2 <- LL^2
  # 1 in R will be recycled to the size of the matrix LL2
  ILL2 <- 1 - LL2
 
  M_prod <- ILL2 %*% L2
 
  # f = sum(diag(crossprod(L2, ILL2 %*% L2)))
  f <- sum(L2 * M_prod)
 
  Gq1 <- 4 * L * M_prod
  Gq2 <- 4 * (LL * tcrossprod(L2)) %*% L
 
  Gq <- Gq1 - Gq2
 
  list(Gq = Gq, f = f, Method = "Tandem II")
}

# --- OBLIMAX ---

oblimax <- function(A, Tmat = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "oblimax",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.oblimax <- function(L) {
  # Oblimax oblique rotation criterion.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value (negated for minimization)
  #     Method - method label

  s2 <- sum(L^2)
  s4 <- sum(L^4)

  list(
    Gq     = -(4 * L^3 / s4 - 4 * L / s2),
    f      = -(log(s4) - 2 * log(s2)),
    Method = "Oblimax"
  )
}

# --- BENTLER ---

bentlerT <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bentler",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

bentlerQ <- function(A, Tmat = NULL, normalize = FALSE,
                     randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "bentler",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.bentler <- function(L) {
  # Bentler's invariant pattern simplicity oblique rotation criterion.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (negated for minimization)
  #     f      - objective function value (negated for minimization)
  #     Method - method label

  L2 <- L^2
  M <- crossprod(L2)
  dM <- diag(M)
  
  # Safer log-determinant
  logDetM <- determinant(M, logarithm = TRUE)$modulus
  logDetD <- sum(log(dM))
    
  f <- -(logDetM - logDetD) / 4
    
  # Gq logic remains correct, but solve() can be replaced 
  # with solve(M, L2) for better stability in some contexts
  D_inv <- diag(1/dM)
  Gq <- -L * (L2 %*% (solve(M) - D_inv))
    
  list(Gq = Gq, f = as.numeric(f), Method = "Bentler")
}

# --- QUARTIMAX ---

quartimax <- function(A, Tmat = NULL, normalize = FALSE,
                      randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "quartimax",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.quartimax <- function(L) {
  # Quartimax orthogonal rotation criterion.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (negated for minimization)
  #     f      - objective function value (negated for minimization)
  #     Method - method label

  list(
    Gq     = -L^3,
    f      = -sum(L^4) / 4,
    Method = "Quartimax"
  )
}

# --- VARIMAX ---

Varimax <- function(A, Tmat = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "varimax",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.varimax <- function(L) {
  # Varimax orthogonal rotation criterion (Kaiser, 1958).
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix (negated for minimization)
  #     f      - objective function value (negated for minimization)
  #     Method - method label

  QL <- L^2 - rep(colMeans(L^2), each = nrow(L))

  list(
    Gq     = -L * QL,
    f      = -sum(QL^2) / 4,
    Method = "Varimax"
  )
}

# --- BINORMAMIN ---

binormamin <- function(A, Tmat = NULL, normalize = FALSE,
					   randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "binormamin",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = FALSE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)			   
					   }

vgQ.binormamin <- function(L) {
  p <- nrow(L)
  m <- ncol(L)
 
  L2    <- L^2
  u     <- colSums(L2)
  inv_u <- 1 / u
 
  # Cross-product matrix of squared loadings: N[k, q] = sum_i L_ik^2 * L_iq^2
  N <- crossprod(L2)
 
  # Outer product of column squared norms: U[k, q] = u_k * u_q
  U <- outer(u, u)
 
  # Ratio matrix C[k, q] = N[k, q] / U[k, q]
  C <- N / U
 
  # Objective value (minimization): sum of off-diagonal elements
  f <- sum(C) - sum(diag(C))
 
  # Vectorized gradient matrix calculation
  L2_over_u   <- sweep(L2, 2, inv_u, "*")
  sum_L2_u    <- rowSums(L2_over_u)
  t1_mat      <- sweep(sum_L2_u - L2_over_u, 2, inv_u, "*")
 
  N_over_u    <- sweep(N, 2, inv_u, "*")
  sum_N_k_vec <- rowSums(N_over_u) - diag(N_over_u)
  t2_vec      <- sum_N_k_vec * (inv_u^2)
 
  Gq <- 4 * L * sweep(t1_mat, 2, t2_vec, "-")
 
  list(f = f, Gq = Gq, Method = "Binormamin")
}

# --- VARIMIN ---

varimin <- function(A, Tmat = NULL, normalize = FALSE,
                    randomStarts = 0, algorithm = "bb", fwindow = 10, ...) {
  mc <- match.call()
  .GPA_RS_engine(A = A, Tmat = Tmat, method = "varimin",
                 normalize    = normalize,
                 randomStarts = randomStarts,
                 orthogonal   = TRUE, call = mc, 
                 algorithm = algorithm, fwindow = fwindow, ...)
}

vgQ.varimin <- function(L) {
  # Varimin orthogonal rotation criterion.
  # Minimizes the variance of squared loadings within factors - 
  # the complement of Varimax.
  #
  # Args:
  #   L : Factor loading matrix
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  QL <- L^2 - rep(colMeans(L^2), each = nrow(L))

  list(
    Gq     = L * QL,
    f      = sum(QL^2) / 4,
    Method = "Varimin"
  )
}

# --- LP WLS ---

vgQ.lp.wls <- function(L, W) {
  # Weighted least squares criterion for Lp rotation.
  #
  # Args:
  #   L : Factor loading matrix
  #   W : Weight matrix (same dimensions as L)
  #
  # Returns:
  #   A list with:
  #     Gq     - gradient matrix
  #     f      - objective function value
  #     Method - method label

  p <- nrow(L)

  list(
    Gq     = 2 * W * L / p,
    f      = sum(W * L^2) / p,
    Method = "Weighted least squares for Lp rotation"
  )
}