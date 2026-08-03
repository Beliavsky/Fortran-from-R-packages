#################################################################################################################
#' Compute the auxiliary parameters a.
#'
#' @param ParaValue     Parameter value
#' @param Const_Type_Full    character-based vector that describes the constraints. Constraints are:
#'                      \itemize{
#'                      \item 'Jordan' for single-country models;
#'                      \item 'Jordan; stationary' for single-country models;
#'                      \item 'Jordan MultiCountry' for multicountry models;
#'                      \item 'Jordan MultiCountry; stationary' for multicountry models;
#'                      \item 'psd' for JPS-based models;
#'                      \item 'BlockDiag' for GVAR-based models;
#'                      \item 'JLLstructure' for JLL-based models;
#'                        }
#' @param Economies       string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels   string-list based which contains the labels of all the variables present in the model
#' @param JLLinputs     list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#' @importFrom hablar s
#'
#' @keywords internal

GetAuxPara <- function(ParaValue, Const_Type_Full, Economies, FactorLabels, JLLinputs = NULL) {
  Const_Type <- Adjust_Const_Type(Const_Type_Full)

  # CASE 1 : Jordan-related constraints (K1XQ)
  if (grepl("Jordan", Const_Type)) {
    a <- Aux_Jordan(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs)
    # CASE 2: psd matrix (SSZ)
  } else if (Const_Type == "psd") {
    a <- Aux_PSD(ParaValue, Const_Type)
    # CASE 3: Block diagonal matrix (SSZ)
  } else if (Const_Type == "BlockDiag") {
    a <- Aux_BlockDiag(ParaValue, Const_Type, FactorLabels, Economies)
    # CASE 4: JLL structure of Sigma matrix (SSZ)
  } else if (Const_Type == "JLLstructure") {
    a <- Aux_JLLstruct(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs)
  }

  for (j in seq_len(length(a))) {
    a[j] <- min(max(s(c(Re(a[j]), -1e20))), 1e20)
  }
  return(a)
}

#########################################################################################################
#' Adjust the constant label
#'
#' @param Const_Type_Full Type of constraint (Full Name)
#'
#' @keywords internal

Adjust_Const_Type <- function(Const_Type_Full) {
  constraints <- c(
    "Jordan MultiCountry; stationary", "Jordan MultiCountry", "Jordan; stationary", "Jordan",
    "psd", "BlockDiag", "JLLstructure"
  )
  for (constraint in constraints) {
    if (grepl(constraint, Const_Type_Full)) {
      return(constraint)
    }
  }
}

##########################################################################################################
#' Transformation of the Jordan-related parameters (auxiliary form)
#'
#' @param ParaValue Constrained parameter
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#'
#' @keywords internal

Aux_Jordan <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs) {
  # 1) SINGLE COUNTRY SETUPS
  if (Const_Type %in% c("Jordan", "Jordan; stationary")) {
    j_aux <- Aux_jordan_OneCountry(ParaValue, Const_Type)

    # 2) MULTI COUNTRY SETUPS (Jordan MultiCountry)
  } else if (any(Const_Type == c("Jordan MultiCountry", "Jordan MultiCountry; stationary"))) {
    N <- length(FactorLabels$Spanned)
    C <- length(Economies)

    # Adjustment for JLL models
    if (!is.null(JLLinputs)) {
      ParaValue <- Jordan_JLL(ParaValue, C, N)
    }

    all_blocks <- list()
    idx0 <- 0

    for (i in seq_len(C)) {
      idx1 <- idx0 + N

      K1Q <- ParaValue[(idx0 + 1):idx1, (idx0 + 1):idx1] # Extract submatrix for economy i
      block_jordan <- Aux_jordan_OneCountry(K1Q, Const_Type) # Extract submatrix for economy i

      all_blocks[[i]] <- block_jordan
      idx0 <- idx1
    }

    # Stack all results into one matrix
    j_aux <- do.call(rbind, all_blocks)
  }

  return(j_aux)
}
###########################################################################################################
#' Auxiliary function for a single-country specification
#'
#' @param ParaValue Constrained parameter
#' @param Const_Type Type of constraint
#'
#' @references
#' Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
#' Estimation of Some Term Structure Models. EABCN Training School.\cr
#' This function offers an independent R implementation that is informed
#' by the conceptual framework outlined in Le and Singleton (2018), but adapted to the
#' present modeling context.
#'
#' @keywords internal

Aux_jordan_OneCountry <- function(ParaValue, Const_Type) {
  # a) Eigen decomposition
  eig_vals <- eigen(ParaValue, only.values = TRUE)$values

  # b) Impose stationarity restriction if requested
  if (grepl("stationary", Const_Type)) {
    eig_vals <- ImposeStat_Aux(eig_vals)
  }

  # c) Separate real and complex eigenvalues
  real_vals <- Re(eig_vals[Im(eig_vals) == 0])
  complex_vals <- eig_vals[Im(eig_vals) != 0]

  # d) Handle odd number of real eigenvalues (keep one aside)
  jordan_params <- numeric(0)

  if (length(real_vals) %% 2 == 1) {
    jordan_params <- c(jordan_params, real_vals[1])
    real_vals <- real_vals[-1] # remove the first element
  }

  # e) Process real eigenvalue pairs
  if (length(real_vals) > 0) {
    for (h in seq(1, length(real_vals), by = 2)) {
      mean_pair <- 0.5 * (real_vals[h] + real_vals[h + 1])
      diff_sq <- (0.5 * (real_vals[h] - real_vals[h + 1]))^2
      jordan_params <- c(jordan_params, mean_pair, diff_sq)
    }
  }

  # f) Process complex eigenvalue pairs
  if (length(complex_vals) > 0) {
    for (h in seq(1, length(complex_vals), by = 2)) {
      real_part <- Re(complex_vals[h])
      neg_imag_sq <- -abs(Im(complex_vals[h]))^2
      jordan_params <- c(jordan_params, real_part, neg_imag_sq)
    }
  }

  # g) Ensure result is a column matrix
  Jordan_mat <- matrix(jordan_params, ncol = 1)

  return(Jordan_mat)
}

###########################################################################################################
#' Impose stationary constraint under the risk-neutral measure
#'
#' @param yy numerical vector before imposing stationary constraint
#' @param ub upper bound
#' @param lb lower bound
#'
#' @keywords internal

ImposeStat_Aux <- function(yy, ub = 0.9999, lb = 1e-4) {
  max_val <- max(abs(yy))
  clamped_val <- min(max(max_val, lb), ub)
  bound <- (clamped_val - lb) / (ub - clamped_val)
  max_adj <- bound + log1p(-expm1(-bound) - 1)

  scale_factor <- max_adj / clamped_val
  yy_scaled <- yy * scale_factor

  if (any(is.nan(Im(yy_scaled)))) {
    yy_scaled[is.nan(Im(yy_scaled))] <- Re(yy_scaled[is.nan(Im(yy_scaled))]) + 0i # replaces the NaN in the imaginary part by 0
  }
  return(yy_scaled)
}

##########################################################################################################
#' Check for JLL models for Jordan restrictions (auxiliary form)
#'
#' @param ParaValue Constrained parameter value
#' @param C number of countries of the economic system
#' @param N number of country-specific spanned factors
#'
#' @keywords internal

Jordan_JLL <- function(ParaValue, C, N) {
  IDXMaxeigen <- seq(1, N * C, by = N) # identify the indexes of the largest eigenvalue of each country
  MaxeigenCS <- diag(ParaValue)[IDXMaxeigen]

  if (all(MaxeigenCS > 0.9999)) {
    BB <- which.min(MaxeigenCS)
    ParaValue[IDXMaxeigen[BB], IDXMaxeigen[BB]] <- 0.9998 # Replace the eigenvalue whose value is closest to 1 by 0.9998
  }

  return(ParaValue)
}
##########################################################################################################
#' Transformation of a PSD matrix (auxiliary form)
#'
#' @param ParaValue Constrained parameter value
#' @param Const_Type Type of constraint
#'
#' @keywords internal

Aux_PSD <- function(ParaValue, Const_Type) {
  # Compute auxiliary PSD matrix
  N <- dim(ParaValue)[1]

  # Compute a matrix square root such that Y=X^(1/2)*X^(1/2)
  MatSSZ <- as.matrix(ParaValue) # Useful for the case in which SSZ is a scalar

  eig <- eigen(MatSSZ)
  V <- eig$vectors
  vv <- eig$values
  if (length(vv) == 1) {
    D <- matrix(sqrt(abs(vv)), nrow = 1, ncol = 1)
  } else {
    D <- diag(sqrt(abs(vv)))
  }
  MatSSZ_Half <- V %*% D %*% solve(V) # V %*% D^(1/2) %*% V^(-1)

  # Vectorize
  MatOnes <- matrix(1, nrow = N, ncol = N)
  mat_SSZ_Half_vec <- t(t(MatSSZ_Half[which(MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))]))

  return(mat_SSZ_Half_vec)
}

######################################################################################################## "
#' Transformation of the JLL-related parameters (auxiliary form)
#'
#' @param ParaValue Constrained parameter value
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#' @keywords internal

Aux_JLLstruct <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs) {
  # Transform JLL parameters
  C <- length(Economies)
  Nspa <- length(FactorLabels$Spanned)
  Macro <- length(FactorLabels$Domestic) - Nspa
  G <- length(FactorLabels$Global)
  K <- C * (Macro + Nspa) + G

  ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$Sigma_Ye
  MatOnes <- matrix(1, ncol = K, nrow = K)
  MatOnes[ZeroIdxSigmaJLL] <- 0
  IdxNONzeroSigmaJLL <- which(MatOnes != 0 & MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))

  # Compute a matrix square root such that Y=X^(1/2)*X^(1/2).
  m <- as.matrix(ParaValue) # Useful for the case in which SSZ is a scalar
  eig <- eigen(m)
  V <- eig$vectors
  vv <- eig$values
  D <- diag(sqrt(abs(vv)))

  halfm <- V %*% D %*% solve(V) # V %*% D^(1/2) %*% V^(-1)

  a <- t(t(halfm[IdxNONzeroSigmaJLL]))

  return(a)
}

#############################################################################################################
#' Transformation of the block diagonal parameters (auxiliary form)
#'
#' @param ParaValue Constrained parameter value
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

Aux_BlockDiag <- function(ParaValue, Const_Type, FactorLabels, Economies) {
  # Extract M from Const_Type
  i <- regexpr("BlockDiag", Const_Type)
  M <- as.numeric(substr(Const_Type, start = i + 9, stop = nchar(Const_Type)))
  if (is.na(M)) M <- 1

  # Get dimensions
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)

  # Initialize storage
  a <- numeric(0)
  idx0 <- 0
  step <- c(G, rep(K, times = C))

  # Compute auxiliary block diagonal matrix
  for (i in seq_len(C + 1)) {
    idx1 <- idx0 + step[i]

    if (idx0 < idx1) {
      seq_indices <- seq(idx0 + 1, idx1)
      d <- as.matrix(ParaValue[seq_indices, seq_indices])

      if (length(d) > 0) {
        N <- nrow(d)
        Mat1s <- matrix(1, nrow = N, ncol = N)
        ZeroIdx <- which(d == 0)

        atemp <- do.call(rbind, lapply(seq_len(M), function(i) {
          MatSSZ <- as.matrix(d[, (N * (i - 1) + 1):(N * i)]) # Useful for the case in which SSZ is a scalar

          eig <- eigen(MatSSZ)
          V <- eig$vectors
          vv <- eig$values
          if (length(vv) == 1) {
            D <- matrix(sqrt(abs(vv)), nrow = 1, ncol = 1)
          } else {
            D <- diag(sqrt(abs(vv)))
          }
          MatSSZ_Half <- V %*% D %*% solve(V) # V %*% D^(1/2) %*% V^(-1)
          MatSSZ_Half[ZeroIdx] <- 0
          ab <- MatSSZ_Half[lower.tri(Mat1s, diag = TRUE)]
          matrix(ab[ab != 0], ncol = 1)
        }))

        a <- c(a, atemp)
      }
    }

    idx0 <- idx1
  }

  return(matrix(a, ncol = 1))
}
