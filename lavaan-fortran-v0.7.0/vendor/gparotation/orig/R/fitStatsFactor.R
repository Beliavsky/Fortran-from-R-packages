#' Calculate Factor Simplicity Index
#'
#' Computes the Factor Simplicity Index (FSI) for a rotated factor solution,
#' following Lorenzo-Seva (2003).
#'
#' @param x A GPArotation object, a factanal object, or a factor loading matrix.
#' @param digits Integer. Number of decimal places for rounding. Default 3.
#'
#' @return A list with:
#'   \item{FSI}{numeric vector of FSI values, one per factor.}
#'   \item{FSI_mean}{mean FSI across factors.}
#'   \item{FSI_min}{minimum FSI across factors.}
#'
#' @references
#' Lorenzo-Seva, U. (2003) A factor simplicity index.
#' \emph{Psychometrika}, \bold{68}(1), 49--60.
#' \doi{10.1007/BF02296652}
#'
#' @export
calc_FSI <- function(x, digits = 3L) {

  # --- Extract loadings ---
  if (inherits(x, "GPArotation")) {
    L <- unclass(x$loadings)
  } else if (inherits(x, "factanal")) {
    L <- unclass(x$loadings)
  } else {
    L <- unclass(x)
  }

  if (!is.matrix(L))
    L <- as.matrix(L)

  p <- nrow(L)
  k <- ncol(L)

  if (p <= 1)
    stop("FSI requires more than one indicator.", call. = FALSE)

  # --- FSI per factor ---
  FSI <- numeric(k)

  for (j in 1:k) {
    l2     <- L[, j]^2
    l2_sum <- sum(l2)
    l4_sum <- sum(l2^2)

    if (l2_sum == 0) {
      FSI[j] <- 0
      next
    }

    FSI[j] <- (p * l4_sum - l2_sum^2) /
              ((p - 1) * l2_sum^2)
  }

  names(FSI) <- colnames(L)

  list(
    FSI      = round(FSI,       digits),
    FSI_mean = round(mean(FSI), digits),
    FSI_min  = round(min(FSI),  digits)
  )
}


#' Calculate AUC Simple Structure Measure
#'
#' Computes the Area Under the Curve (AUC) measure of factor loading
#' simplicity for a GPArotation object or loading matrix, following
#' Liu and Moustaki (2024).
#'
#' @param x A GPArotation object or a factor loading matrix.
#' @param digits Integer. Number of decimal places for rounding. Default 3.
#'
#' @return A list with:
#'   \item{AUC}{numeric vector of AUC values per factor}
#'   \item{AUC_mean}{mean AUC across factors}
#'   \item{AUC_min}{minimum AUC across factors}
#'
#'
#' @export
calc_AUC <- function(x, digits = 3L) {

  extracted <- .extract_A(x)
  L         <- extracted$A
  
  # --- Extract loadings ---
  if (inherits(x, "GPArotation")) {
    L <- unclass(x$loadings)
  } else {
    L <- unclass(x)
  }

  p <- nrow(L)
  k <- ncol(L)

  # --- AUC per factor ---
  AUC <- numeric(k)

  for (j in 1:k) {
    l2       <- L[, j]^2
    l2_sort  <- sort(l2, decreasing = TRUE)
    l2_sum   <- sum(l2)

    if (l2_sum == 0) {
      AUC[j] <- 0
      next
    }

    # Cumulative proportion of variance explained
    # weighted by rank position
    AUC[j] <- sum(cumsum(l2_sort / l2_sum) ) / p
  }

  # --- Baseline adjustment ---
  # Subtract diagonal baseline (0.5) so AUC = 0 means no simple structure
  # and AUC = 0.5 means perfect simple structure
  AUC_adj <- AUC - 0.5

  names(AUC)     <- colnames(L)
  names(AUC_adj) <- colnames(L)

  list(
    AUC      = round(AUC,      digits),
    AUC_adj  = round(AUC_adj,  digits),
    AUC_mean = round(mean(AUC),     digits),
    AUC_min  = round(min(AUC),      digits),
    AUC_adj_mean = round(mean(AUC_adj), digits),
    AUC_adj_min  = round(min(AUC_adj),  digits)
  )
}


#' Calculate Overall Simple Structure Indices
#'
#' Computes three overall measures of simple structure for a rotated factor
#' solution: Hoffman (1978), Gini (Kaiser, 1974), and Bentler (1977).
#' Unlike AUC and FSI which are computed per factor, these indices
#' summarize simple structure across the entire loading matrix.
#'
#' @param x A GPArotation object, a factanal object, or a factor loading matrix.
#' @param digits Integer. Number of decimal places for rounding. Default 3.
#'
#' @return A list with elements:
#'   \item{Hoffman}{Hoffman's simplicity index.}
#'   \item{Gini}{Gini simplicity index.}
#'   \item{Bentler}{Bentler's simplicity index.}
#'
#' @references
#' Bentler, P.M. (1977) Factor simplicity index and transformations.
#'   \emph{Psychometrika}, \bold{42}(2), 277--295.
#'   \doi{10.1007/BF02294054}
#'
#' Hoffman, P.J. (1978) The paramorphic representation of clinical judgment.
#'   \emph{Psychological Bulletin}, \bold{55}(2), 116--131.
#'   \doi{10.1037/h0047807}
#'
#' Kaiser, H.F. (1974) An index of factorial simplicity.
#'   \emph{Psychometrika}, \bold{39}(1), 31--36.
#'   \doi{10.1007/BF02291575}
#'
#' @export
calc_simplicity <- function(x, digits = 3L) {

  # --- Extract loadings ---
  if (inherits(x, "GPArotation")) {
    L <- unclass(x$loadings)
  } else if (inherits(x, "factanal")) {
    L <- unclass(x$loadings)
  } else {
    L <- unclass(x)
  }

  if (!is.matrix(L))
    L <- as.matrix(L)

  p <- nrow(L)
  k <- ncol(L)

  if (p <= 1)
    stop("Simplicity indices require more than one indicator.", call. = FALSE)
  if (k <= 1)
    stop("Simplicity indices require more than one factor.", call. = FALSE)

  L2 <- L^2
  L4 <- L^4

  # --- Hoffman (1978) ---
  # Measures how far loadings are from 0.5
  # H = 1 when all loadings are 0 or 1 (perfect simple structure)
  # H = 0 when all loadings equal 0.5
  H <- 1 - sum(L2 * (1 - L2)) / (p * k * 0.25)

  # --- Gini (Kaiser, 1974) ---
  # Based on Gini coefficient of absolute loadings per factor
  # then averaged across factors
  gini_j <- numeric(k)
  for (j in 1:k) {
    a       <- sort(abs(L[, j]))
    n       <- length(a)
    a_sum   <- sum(a)
    if (a_sum == 0) {
      gini_j[j] <- 0
    } else {
      ranks     <- seq_len(n)
      gini_j[j] <- 1 - (2 * sum((n - ranks + 0.5) * a)) / (n * a_sum)
    }
  }
  G <- mean(gini_j)

  # --- Bentler (1977) ---
  # Ratio of between-factor to within-factor variance of squared loadings
  col_sum_L2  <- colSums(L2)
  total_L2    <- sum(L2)
  total_L4    <- sum(L4)
  sum_col_L2_sq <- sum(col_sum_L2^2)

  numerator   <- sum_col_L2_sq / p - total_L4 / p
  denominator <- total_L2^2 / p  - total_L4 / p

  B <- if (abs(denominator) < .Machine$double.eps)
    NA
  else
    numerator / denominator

  list(
    Hoffman = round(H, digits),
    Gini    = round(G, digits),
    Bentler = round(B, digits)
  )
}



#' Calculate Hyperplane Count
#'
#' Computes the hyperplane count for a rotated factor solution.
#' The hyperplane count is the number of loadings falling within a
#' narrow band around zero, per factor and overall.
#'
#' @param x A GPArotation object, a factanal object, or a factor loading matrix.
#' @param cutoff Numeric. Half-width of the hyperplane band. Default 0.1.
#' @param digits Integer. Number of decimal places for rounding. Default 3.
#'
#' @return A list with elements:
#'   \item{HP}{integer vector of hyperplane counts per factor.}
#'   \item{HP_total}{total hyperplane count across all factors.}
#'   \item{HP_max}{theoretical maximum hyperplane count p*(k-1).}
#'   \item{HP_pct}{percentage of theoretical maximum achieved.}
#'   \item{cutoff}{the hyperplane cutoff used.}
#'
#' @references
#' Thurstone, L.L. (1947) \emph{Multiple Factor Analysis}.
#'   University of Chicago Press.
#'
#' Jennrich, R.I. (2004) Rotation to simple loadings using component loss
#'   functions: The oblique case. \emph{Psychometrika}, \bold{71}(1),
#'   173--191. 
#'
#' @export
calc_hyperplane <- function(x, cutoff = 0.1, digits = 3L) {

  # --- Extract loadings ---
  if (inherits(x, "GPArotation")) {
    L <- unclass(x$loadings)
  } else if (inherits(x, "factanal")) {
    L <- unclass(x$loadings)
  } else {
    L <- unclass(x)
  }

  if (!is.matrix(L))
    L <- as.matrix(L)

  p <- nrow(L)
  k <- ncol(L)

  if (k <= 1)
    stop("Hyperplane count requires more than one factor.", call. = FALSE)

  # --- Hyperplane count per factor ---
  HP       <- colSums(abs(L) < cutoff)
  names(HP) <- colnames(L)

  # --- Overall ---
  HP_total <- sum(HP)
  HP_max   <- p * (k - 1)  # theoretical maximum - each item in k-1 hyperplanes
  HP_pct   <- round(100 * HP_total / HP_max, digits)

  list(
    HP       = HP,
    HP_total = HP_total,
    HP_max   = HP_max,
    HP_pct   = HP_pct,
    cutoff   = cutoff
  )
}
