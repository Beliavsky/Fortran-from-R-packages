#' Estimates an unbiased VAR(1) using stochastic approximation (Bauer, Rudebusch and Wu, 2012)
#'
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param BRWinputs list. Contains the necessary inputs for the BRW model estimation:
#' \enumerate{
#'   \item \code{Cent_Measure}: "Mean" or "Median" (unbiased estimation type)
#'   \item \code{gamma}: Numeric. Adjustment parameter between 0 and 1. Default is 0.5.
#'   \item \code{N_iter}: Integer. Number of iterations for the stochastic approximation algorithm after burn-in. Default is 5000.
#'   \item \code{N_burn}: Integer. Number of burn-in iterations. Default is 15% of \code{N_iter}.
#'   \item \code{B}: Integer. Number of bootstrap samples per iteration for calculating the noisy measure of the biased estimator's mean or median. Default is 50.
#'   \item \code{check}: Logical. Indicates whether to perform a closeness check. Default is TRUE.
#'   \item \code{B_check}: Integer. Number of bootstrap samples for the closeness check. Default is 100000.
#'   \item \code{Eigen_rest}: Numeric. Restriction on the largest eigenvalue under the P-measure. Default is 1.
#' }
#' @param RiskFactors numeric matrix (\code{Td x K}). Time series of risk factors.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param FactorLabels list. Labels for all variables in the model.
#' @param GVARinputs list. Inputs for GVAR model estimation (see \code{\link{GVAR}}). Default is NULL.
#' @param JLLinputs list. Inputs for JLL model estimation (see \code{\link{JLL}}). Default is NULL.
#' @param verbose logical. Flag controlling function messaging. Default TRUE.
#'
#' @examples
#' \donttest{
#' data(RiskFacFull)
#' Factors <- t(RiskFacFull[1:7, ])
#'
#' BRWinputs <- list(
#'   Cent_Measure = "Mean", gamma = 0.4, N_iter = 1000, N_burn = 100,
#'   B = 10, check = 1, B_check = 5000
#' )
#'
#' Economies <- "China"
#' N <- 3
#' ModelType <- "JPS original"
#' FactorLabels <- NULL
#'
#' BRWpara <- Bias_Correc_VAR(ModelType, BRWinputs, Factors, Economies, FactorLabels, verbose = FALSE)
#' }
#'
#' @return Bias-corrected VAR parameters based on the framework of Bauer, Rudebusch and Wu (2012). The list contains:
#' \enumerate{
#'   \item \code{KOZ_BC}: estimated intercept (K x 1);
#'   \item \code{K1Z_BC}: estimated feedback matrix (K x K);
#'   \item \code{SSZ_BC}: estimated variance-covariance matrix (K x K);
#'   \item \code{dist}: root mean square distance (scalar);
#' }
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td} denotes the model time series dimension.
#'   \item \code{C} number of countries in the system.
#'   \item \code{K} denotes the total number of risk factors.
#' }
#'
#' @references
#' Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic Term Structure Models" \cr
#' This function offers an independent R implementation that is informed
#' by the conceptual framework outlined in Bauer, Rudebusch and Wu (2012), but adapted to the
#' present modeling context. Related Matlab routines are available on Cynthia Wu's
#' website (https://sites.google.com/view/jingcynthiawu/).
#' @export

Bias_Correc_VAR <- function(ModelType, BRWinputs, RiskFactors, Economies, FactorLabels, GVARinputs = NULL,
                            JLLinputs = NULL, verbose = TRUE) {
  if (BRWinputs$Cent_Measure == "Mean") {
    Use_Mean <- TRUE
  }

  N <- length(FactorLabels$Spanned)
  gamma <- withDefault(BRWinputs$gamma, 0.5)
  N_iter <- withDefault(BRWinputs$N_iter, 5000)
  N_burn <- withDefault(BRWinputs$N_burn, N_iter * 0.15)
  N_boot <- withDefault(BRWinputs$B, 50)
  check <- withDefault(BRWinputs$check, FALSE)
  B_check <- withDefault(BRWinputs$B_check, 100000)
  checkSigma <- withDefault(BRWinputs$checkSigma, TRUE)
  Eigen_rest <- withDefault(BRWinputs$Eigen_rest, 1)

  # Step 1: Estimation of the P-dynamics feedback matrix (Without bias correction)
  K1Z_NoBC <- Get_FeedMat_NoBC(RiskFactors, ModelType, Economies, GVARinputs, JLLinputs, FactorLabels,
    CheckInpts = FALSE
  )

  # Step 2: stochastic approximation algorithm
  K1Z_BC <- SA_algorithm(
    K1Z_NoBC, RiskFactors, BRWinputs, GVARinputs, JLLinputs, FactorLabels, Economies,
    ModelType, verbose
  )

  # Step 3: Optional closeness check
  if (check) {
    dist <- Check_comparison_NoBC(
      K1Z_BC, K1Z_NoBC, B_check, RiskFactors, GVARinputs, JLLinputs, FactorLabels,
      Economies, ModelType, Use_Mean, verbose
    )
  } else {
    dist <- NaN
  }

  # Step 4: Impose eigenvalue restriction on the feedback matrix
  K1Z_BC <- shrink_FeedMat_BC(K1Z_BC, K1Z_NoBC, Eigen_rest)

  # Step 5: Intercept and variance covariance matrix
  K0Z_BC <- (diag(ncol(RiskFactors)) - K1Z_BC) %*% matrix(colMeans(RiskFactors), ncol = 1)
  SSZ_BC <- if (checkSigma) Get_SSZ_BC(K1Z_BC, RiskFactors, GVARinputs, JLLinputs, FactorLabels, ModelType) else NaN

  output <- list(K0Z_BC = K0Z_BC, K1Z_BC = K1Z_BC, SSZ_BC = SSZ_BC, dist = dist)

  return(output)
}
############################################################################################
#' Check default value
#'
#' @param value parameter to be checked
#' @param default default value
#'
#' @keywords internal

withDefault <- function(value, default) if (is.null(value)) default else value

############################################################################################
#' Estimate feedback matrix from several models (No bias-corrected version)
#'
#' @param RiskFactors A numeric matrix (T x F) representing the time series of risk factors.
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param GVARinputs List. Inputs for GVAR model estimation.
#' @param JLLinputs List. Inputs for JLL model estimation.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param CheckInpts Check for the consistency of the provided inputs. Default is FALSE.
#'
#' @keywords internal

Get_FeedMat_NoBC <- function(RiskFactors, ModelType, Economies, GVARinputs, JLLinputs, FactorLabels,
                             CheckInpts = FALSE) {
  D_RiskFactors <- scale(RiskFactors, scale = FALSE)
  N <- length(FactorLabels$Spanned)

  if (any(ModelType == c("JPS original", "JPS global"))) {
    T_dim <- nrow(D_RiskFactors)
    Y <- t(D_RiskFactors[2:T_dim, , drop = FALSE]) # k × (T-1)
    X <- t(cbind(1, D_RiskFactors[1:(T_dim - 1), , drop = FALSE])) # (k+1) × (T-1)

    # OLS estimation: A = Y X^T (X X^T)^{-1}
    XXt <- X %*% t(X)
    Betas <- Y %*% t(X) %*% solve(XXt, tol = 1e-50)
    K1Z_NoBC <- Betas[, -1]
  } else if (ModelType == "JPS multi") {
    VARpara <- VAR(t(D_RiskFactors), VARtype = "unconstrained")
    K1Z_NoBC <- VARpara$K1Z
  } else if (any(ModelType == c("GVAR single", "GVAR multi"))) {
    GVARinputs$GVARFactors <- DataSet_BS(ModelType, t(D_RiskFactors), GVARinputs$Wgvar, Economies, FactorLabels)
    GVARpara <- GVAR(GVARinputs, N)
    K1Z_NoBC <- GVARpara$F1
  } else {
    JLLinputs$WishSigmas <- 0
    JLLPara <- JLL(t(D_RiskFactors), N, JLLinputs)
    K1Z_NoBC <- JLLPara$k1
  }

  return(K1Z_NoBC)
}
################################################################################################
#' Stochastic approximation algorithm
#'
#' @param K1Z_NoBC feedback matrix before bias-correction
#' @param RiskFactors A numeric matrix (T x F) representing the time series of risk factors.
#' @param BRWlist A list containing the necessary inputs for the BRW model estimation
#' @param GVARinputs List. Inputs for GVAR model estimation.
#' @param JLLinputs List. Inputs for JLL model estimation.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param verbose verbose Logical flag controlling function messaging.
#'
#' @keywords internal

SA_algorithm <- function(K1Z_NoBC, RiskFactors, BRWlist, GVARinputs, JLLinputs, FactorLabels, Economies,
                         ModelType, verbose) {
  if (BRWlist$Cent_Measure == "Mean") {
    Use_Mean <- TRUE
  }

  N_burn <- BRWlist$N_burn
  N_iter <- BRWlist$N_iter
  N_boot <- BRWlist$B
  gamma <- BRWlist$gamma

  N <- length(FactorLabels$Spanned)

  K1Z_All <- list()
  K1Z_All[[1]] <- K1Z_NoBC

  j <- 1 # numbers of accepted draws
  w <- 1 # index for printing on screen

  while (j <= N_burn + N_iter - 1) {
    if (j == 100 * w) {
      if (verbose) message(paste("iteration: ", j, " / ", N_burn + N_iter - 1))
      w <- w + 1
    }

    K1Z_j_M <- FeedMat_M(
      K1Z_All[[j]], N_boot, RiskFactors, GVARinputs, JLLinputs, FactorLabels, Economies,
      ModelType, Use_Mean
    )
    d <- K1Z_NoBC - K1Z_j_M
    K1Z_All[[j + 1]] <- K1Z_All[[j]] + gamma * d

    j <- j + 1
  }

  K1Z_subset <- K1Z_All[(N_burn + 1):(N_burn + N_iter)] # Exclude the first N_burn draws
  K1Z_BC <- Reduce("+", K1Z_subset) / length(K1Z_subset) # Bias-Corrected feedback matrix

  return(K1Z_BC)
}

################################################################################################
#' Computes an average or median feedback matrix across several bootstrap iterations
#'
#' @param K1Z_j Feedback matrix at the j_th iteration
#' @param N_boot Number of bootstrap samples per iteration.
#' @param RiskFactors A numeric matrix (T x F) representing the time series of risk factors.
#' @param GVARinputs List. Inputs for GVAR model estimation.
#' @param JLLinputs List. Inputs for JLL model estimation.
#' @param FactorLabels A list of character vectors with labels for all variables in the model
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param Use_Mean  Choose between the mean or the median across the bootstrap iterations.
#'
#' @keywords internal

FeedMat_M <- function(K1Z_j, N_boot, RiskFactors, GVARinputs, JLLinputs, FactorLabels, Economies,
                      ModelType, Use_Mean) {
  N <- length(FactorLabels$Spanned)

  # 1) Compute artificial series
  Art_Series <- Gen_art_series(K1Z_j, N_boot, RiskFactors)

  # 2) N_boot Bootstrap
  K1Z_AllSamples <- list()

  for (n in 1:N_boot) {
    RF_Boot <- t(Art_Series[, , n])
    K1Z_n <- Get_FeedMat_NoBC(RF_Boot, ModelType, Economies, GVARinputs, JLLinputs, FactorLabels)
    K1Z_AllSamples[[n]] <- K1Z_n
  }

  # 3) Compute the media or the mean across several bootstrap iterations
  if (Use_Mean) {
    K1Z_M <- Reduce("+", K1Z_AllSamples) / length(K1Z_AllSamples)
  } else {
    K1Z_temp <- simplify2array(K1Z_AllSamples)
    K1Z_M <- apply(K1Z_temp, c(1, 2), stats::median)
  }

  return(K1Z_M)
}
#############################################################################################################
#' Simulate N_Boot dataset from the P-dynamics
#'
#' @param K1Z_j Feedback matrix at the j_th iteration
#' @param N_Boot Number of bootstrap samples per iteration.
#' @param RFs A numeric matrix (T x F) representing the time series of risk factors.
#'
#' @keywords internal

Gen_art_series <- function(K1Z_j, N_Boot, RFs) {
  T_dim <- nrow(RFs)
  K <- ncol(RFs)
  RFs_mean <- colMeans(RFs)

  # Precompute centered RFs and residuals
  RFs_centered <- sweep(RFs, 2, RFs_mean, "-")
  Resid_org <- t(apply(RFs, 1, function(row) row - RFs_mean - K1Z_j %*% (row - RFs_mean)))

  # Sample initial values
  RFs_art <- array(NA, dim = c(K, T_dim, N_Boot))
  IDX_0 <- sample(T_dim, N_Boot, replace = TRUE)
  RFs_art[, 1, ] <- t(RFs[IDX_0, , drop = FALSE])

  # Precompute the K1Z_j %*% RFs_mean term
  K1Z_mean <- K1Z_j %*% RFs_mean

  # Simulation loop
  for (tt in 2:T_dim) {
    Idx_Random <- sample(T_dim - 1, N_Boot, replace = TRUE)

    # Get residuals for all bootstraps at once
    u_art <- t(Resid_org[Idx_Random, , drop = FALSE])

    # Extract previous period values for all bootstraps
    prev_values <- RFs_art[, tt - 1, , drop = FALSE]
    dim(prev_values) <- c(K, N_Boot)

    # Compute the VAR component for all bootstraps simultaneously
    var_component <- K1Z_j %*% prev_values

    # Update all bootstraps at once
    RFs_art[, tt, ] <- RFs_mean + K1Z_j %*% (prev_values - RFs_mean) + u_art
  }

  rownames(RFs_art) <- colnames(RFs)
  return(RFs_art)
}

########################################################################################################
#' Shrinking the largest eigenvalue
#'
#' @param K1Z_BC VAR (1) bias-corrected feedback matrix from Bauer, Rudebusch and, Wu (2012)
#' @param K1Z_NoBC VAR (1) with no bias-corrected feedback matrix from the selected ATSM
#' @param ev_restr maximum eigenvalue desired in the feedback matrix after the adjustment
#'
#' @keywords internal

shrink_FeedMat_BC <- function(K1Z_BC, K1Z_NoBC, ev_restr) {
  # Extract eigenvalues
  ev_bc <- max(Mod(eigen(K1Z_BC, only.values = TRUE)$values))
  ev_Nobc <- max(Mod(eigen(K1Z_NoBC, only.values = TRUE)$values))

  # Early return if no shrinkage needed
  if (ev_bc <= ev_restr) {
    return(K1Z_BC)
  }

  if (ev_Nobc >= ev_restr) {
    return(K1Z_NoBC)
  }

  # Shrink Feedback matrix towards no bias-correction
  Phi_diff <- K1Z_NoBC - K1Z_BC

  # Use binary search for faster convergence
  delta_min <- 0
  delta_max <- 1
  delta <- 1
  max_iter <- 200 # Prevent infinite loops
  iter <- 0

  while ((delta_max - delta_min) > 1e-6 && iter < max_iter) {
    delta <- (delta_min + delta_max) / 2
    K1Z_test <- K1Z_BC + delta * Phi_diff
    ev_test <- max(Mod(eigen(K1Z_BC, only.values = TRUE)$values))

    if (ev_test <= ev_restr) {
      delta_min <- delta # This delta works, try larger
      K1Z_BC <- K1Z_test
      ev_bc <- ev_test
    } else {
      delta_max <- delta # This delta doesn't work, try smaller
    }
    iter <- iter + 1
  }

  return(K1Z_BC)
}
##############################################################################################################
#' check how close the mean or median of the bias-corrected aproach is from the non-corrected approach
#'
#' @param K1Z_BC feedback matrix after the bias correction procedure
#' @param K1Z_NoBC feedback matrix before the bias correction procedure
#' @param B_check number of bootstrap samples used in the closeness check
#' @param RiskFactors time series of the risk factors (F x T)
#' @param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#' @param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param Use_Mean Choose between the mean or the median across the bootstrap iterations.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

Check_comparison_NoBC <- function(K1Z_BC, K1Z_NoBC, B_check, RiskFactors, GVARinputs, JLLinputs, FactorLabels,
                                  Economies, ModelType, Use_Mean, verbose) {
  K1Z_NoBC <- FeedMat_M(
    K1Z_BC, B_check, RiskFactors, GVARinputs, JLLinputs, FactorLabels, Economies,
    ModelType, Use_Mean
  )
  dist <- sqrt(mean((K1Z_BC - K1Z_NoBC)^2))
  if (verbose) message(paste("Root mean square distance: ", round(dist, 6)))
  if (verbose) message("... Done! \n")

  return(dist)
}


##############################################################################################################
#' Compute the variance-covariance matrix after the bias correction procedure
#'
#' @param K1Z_BC Feedback matrix resulting from the bias-correction procedure
#' @param RiskFactors time series of the risk factors (T x F)
#' @param GVARinputs inputs used in the estimation of the GVAR-based models (see "GVAR" function). Default is set to NULL
#' @param JLLinputs inputs used in the estimation of the JLL-based models (see "JLL" function). Default is set to NULL
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#'
#' @keywords internal

Get_SSZ_BC <- function(K1Z_BC, RiskFactors, GVARinputs, JLLinputs, FactorLabels, ModelType) {
  N <- length(FactorLabels$Spanned)

  if (any(ModelType == c("JPS original", "JPS global", "JPS multi"))) {
    T_dim <- nrow(RiskFactors)
    Xdem <- scale(RiskFactors, scale = FALSE)
    resid_tilde <- t(Xdem[2:T_dim, ]) - K1Z_BC %*% t(Xdem[1:(T_dim - 1), ])
    V_tilde <- resid_tilde %*% t(resid_tilde) / (T_dim - 1)
  } else if (any(ModelType == c("GVAR single", "GVAR multi"))) {
    Economies <- GVARinputs$Economies # necessary to ensure that the code works for the "GVAR single" model
    GVARinputs$GVARFactors <- DataSet_BS(ModelType, t(RiskFactors), GVARinputs$Wgvar, Economies, FactorLabels)
    V_tilde <- GVAR(GVARinputs, N)$Sigma_y
  } else if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    JLLinputs$WishSigmas <- 1
    JLLPara <- JLL(t(RiskFactors), N, JLLinputs)
    V_tilde <- if (any(ModelType == c("JLL original", "JLL No DomUnit"))) JLLPara$Sigmas$VarCov_NonOrtho else JLLPara$Sigmas$VarCov_Ortho
  }

  return(V_tilde)
}
