#' Use function ML to generate the outputs from a ATSM
#'
#' @param x vector containing all the vectorized auxiliary parameters
#' @param ML_fun vector-valued objective function (function)
#' @param ListInputSet variable inputs used in the optimization (see inputs from "optimization" function)
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models (see "JLL" function)
#' @param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models (see "GVAR" function)
#' @param WithEstimation if TRUE, returns only the values of the likelihood function, else generates the entire set of outputs
#'
#'
#' @keywords internal

FunctionML_vec <- function(x, ML_fun, ListInputSet, ModelType, FactorLabels, Economies, JLLinputs,
                           GVARinputs, WithEstimation) {
  Para_Upds <- Update_ParaList(x, ModelType, FactorLabels, Economies, JLLinputs, GVARinputs, ListInputSet)

  # Perform the numerical optimization
  if (WithEstimation) {
    if (any(ModelType == c("JLL original", "JLL No DomUnit"))) {
      MLE <- ML_fun(K1XQ = Para_Upds$K1XQ$Value, ExportListOut = FALSE)
    } else {
      MLE <- ML_fun(K1XQ = Para_Upds$K1XQ$Value, SSZ = Para_Upds$SSZ$Value, ExportListOut = FALSE)
    }
  } else {
    # Gather the estimated outputs after the optimization
    if (any(ModelType == c("JLL original", "JLL No DomUnit"))) {
      ListOut <- ML_fun(K1XQ = Para_Upds$K1XQ$Value)
    } else {
      ListOut <- ML_fun(K1XQ = Para_Upds$K1XQ$Value, SSZ = Para_Upds$SSZ$Value)
    }
  }

  if (WithEstimation) {
    return(MLE)
  } else {
    return(ListOut)
  }
}

###########################################################################################################
#' Update parameters in the optimization process
#'
#' @param x  vector containing all the vectorized auxiliary parameters
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs Set of necessary inputs used in the estimation of the JLL-based models
#' @param GVARinputs Set of necessary inputs used in the estimation of the GVAR-based models
#' @param ListInputSet variable inputs used in the optimization (see "Optimization" function)
#'
#' @keywords internal

Update_ParaList <- function(x, ModelType, FactorLabels, Economies, JLLinputs = NULL, GVARinputs = NULL, ListInputSet) {
  # 1) Identify indices from K1XQ and SSZ
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Dim_K1XQ <- length(FactorLabels$Spanned)
  } else {
    C <- length(Economies)
    Dim_K1XQ <- C * length(FactorLabels$Spanned)
    if (any(ModelType == c("JLL original", "JLL No DomUnit"))) {
      Para_Idx <- Dim_K1XQ
      N_Para <- 1
    }
  }

  if (!any(ModelType == c("JLL original", "JLL No DomUnit"))) {
    Dim_SSZ <- length(x) - Dim_K1XQ
    Para_Idx <- c(Dim_K1XQ, Dim_SSZ)
    N_Para <- 2
  }

  # 2) Obtain true parameter
  j <- 0
  for (i in 1:N_Para) {
    ParaTemp <- x[(j + 1):(j + Para_Idx[i])]
    ListInputSet[[i]]$Value <- GetTruePara(
      ParaTemp, ListInputSet[[i]]$Label, FactorLabels, Economies,
      JLLinputs, GVARinputs
    )
    j <- j + Para_Idx[1]
  }

  return(ListInputSet)
}

##############################################################################################################
#' Map auxiliary (unconstrained) parameters a to constrained parameters b
#'
#' @param ParaValue unconstrained auxiliary parameter
#' @param Const_Type_Full One of the following options:
#'             \itemize{
#'             \item 'Jordan'
#'             \item 'Jordan; stationary'
#'             \item 'Jordan MultiCountry'
#'             \item 'Jordan MultiCountry; stationary'
#'             \item 'psd';
#'             \item 'BlockDiag'
#'             \item 'JLLstructure'
#'             }
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs Inputs used in the estimation of the JLL-based models
#' @param GVARinputs Inputs used in the estimation of the GVAR-based models
#'
#'
#' @keywords internal

GetTruePara <- function(ParaValue, Const_Type_Full, FactorLabels, Economies, JLLinputs = NULL, GVARinputs = NULL) {
  a <- Re(ParaValue)
  Const_Type <- Adjust_Const_Type(Const_Type_Full)

  # CASE 1 : Jordan-related constraints
  if (grepl("Jordan", Const_Type)) {
    b <- True_Jordan(a, Const_Type, FactorLabels, Economies)
    # CASE 2: psd matrix
  } else if (Const_Type == "psd") {
    b <- True_PSD(a, Const_Type)
    # CASE 3: Block diagonal matrix
  } else if (Const_Type == "BlockDiag") {
    b <- True_BlockDiag(a, Const_Type, FactorLabels, Economies, GVARinputs)
    # CASE 4: JLL structure of Sigma matrix
  } else if (Const_Type == "JLLstructure") {
    b <- True_JLLstruct(a, Const_Type, FactorLabels, Economies, JLLinputs)
  }

  return(b)
}

########################################################################################################
#' Transformation of the Jordan-related parameters (True form)
#'
#' @param ParaValue Constrained parameter
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

True_Jordan <- function(ParaValue, Const_Type, FactorLabels, Economies) {
  stationary <- grepl("stationary", Const_Type)

  # 1) SINGLE COUNTRY SETUPS
  if (Const_Type %in% c("Jordan", "Jordan; stationary")) {
    K1Q <- True_jordan_OneCountry(ParaValue, stationary)

    # 2) MULTI COUNTRY SETUPS (Jordan MultiCountry)
  } else if (Const_Type %in% c("Jordan MultiCountry", "Jordan MultiCountry; stationary")) {
    C <- length(Economies)
    N <- length(FactorLabels$Spanned)

    idx0 <- 0
    blocks <- list()

    for (j in seq_len(C)) {
      idx1 <- idx0 + N
      lQ <- ParaValue[(idx0 + 1):idx1]

      blocks[[j]] <- True_jordan_OneCountry(lQ, stationary)

      idx0 <- idx1
    }

    K1Q <- Reduce(adiag, blocks)
  }

  return(K1Q)
}
##########################################################################################################
#' True function for a single-country specification
#'
#' @param ParaValue Constrained parameter
#' @param stationary impose stationarity. Default is FALSE
#'
#' @references
#' Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
#' Estimation of Some Term Structure Models. EABCN Training School.\cr
#' This function offers an independent R implementation that is informed
#' by the conceptual framework outlined in Le and Singleton (2018), but adapted to the
#' present modeling context.
#'
#' @keywords internal

True_jordan_OneCountry <- function(ParaValue, stationary = FALSE) {
  N <- length(ParaValue)

  # Case: scalar input
  if (N == 1) {
    K1Q <- 0
  } else {
    K1Q <- rbind(rep(0, N), diag(1, N - 1, N))
  }

  offset <- 0
  if (N %% 2 != 0) {
    K1Q[1] <- ParaValue[1]
    offset <- 1
    if (!stationary) {
      ParaValue <- ParaValue[-1]
    }
  }

  if (length(ParaValue) > 0) {
    for (i in seq_len(length(ParaValue) / 2)) {
      avg_val <- ParaValue[2 * i - 1]
      diff_val <- ParaValue[2 * i]

      K1Q[offset + 2 * i - 1, offset + 2 * i - 1] <- avg_val
      K1Q[offset + 2 * i, offset + 2 * i] <- avg_val
      K1Q[offset + 2 * i - 1, offset + 2 * i] <- diff_val
    }
  }

  # Apply stationarity adjustment if requested
  if (stationary) {
    K1Q <- ImposeStat_True(ParaValue, K1Q)
  }

  return(K1Q)
}


###########################################################################################################
#' Makes sure that the stationary constraint under the risk-neutral measure is preserved
#'
#' @param params parameter of interest (scalar or matrix)
#' @param K1Q risk-neutral feedback matrix
#' @param ub upper bound. Default is 0.9999.
#' @param lb lower bound. Default is 1e-4.
#'
#' @keywords internal

ImposeStat_True <- function(params, K1Q, ub = 0.9999, lb = 1e-4) {
  # Eigenvalue scaling factor
  eigvals <- eigen(K1Q, only.values = TRUE)$values
  max_eig <- max(abs(eigvals))

  z <- pos_map(max_eig)
  scale_factor <- ((z / (1 + z)) * (ub - lb) + lb) / max_eig

  # Rescale parameters
  params <- params * scale_factor
  N <- length(params)

  # Adjust alternating entries depending on odd/even
  if (N > 1) {
    if (N %% 2 == 0) {
      params[seq(2, N, by = 2)] <- params[seq(2, N, by = 2)] * scale_factor
    } else {
      params[seq(3, N, by = 2)] <- params[seq(3, N, by = 2)] * scale_factor
    }
  }

  # Rebuild Jordan matrix
  i0 <- 0
  if (N %% 2 != 0) {
    K1Q[1] <- params[1]
    params <- params[-1]
    i0 <- 1
  }

  for (i in seq_len(length(params) / 2)) {
    mean_val <- params[2 * i - 1]
    diff_val <- params[2 * i]
    K1Q[i0 + 2 * i - 1, i0 + 2 * i - 1] <- mean_val
    K1Q[i0 + 2 * i, i0 + 2 * i] <- mean_val
    K1Q[i0 + 2 * i - 1, i0 + 2 * i] <- diff_val
  }

  return(K1Q)
}

##########################################################################################################
#' Exponential transformation
#'
#' @param x scalar (numeric)
#'
#' @keywords internal

pos_map <- function(x) {
  y <- numeric(length(x))
  neg_idx <- x < 0
  pos_idx <- !neg_idx

  if (any(neg_idx)) {
    y[neg_idx] <- log1p(expm1(x[neg_idx]) + 1)
  }
  if (any(pos_idx)) {
    y[pos_idx] <- x[pos_idx] + log1p(expm1(-x[pos_idx]) + 1)
  }
  return(y)
}

##########################################################################################################
#' Transformation of a PSD matrix (true form)
#'
#' @param ParaValue Constrained parameter value
#' @param Const_Type Type of constraint
#'
#' @keywords internal

True_PSD <- function(ParaValue, Const_Type) {
  # Dimension of covariance matrix
  N <- floor(sqrt(2 * length(ParaValue)))

  # Fill symmetric matrix with params
  Mat <- matrix(0, N, N)
  Mat[lower.tri(Mat, diag = TRUE)] <- ParaValue
  Mat <- Mat + t(Mat) - diag(diag(Mat)) # enforce symmetry

  # Construct PSD matrix
  Mat_psd <- Mat %*% t(Mat) # Note: the maximization is made with parameters of SSP^(1/2)

  return(Mat_psd)
}

#########################################################################################################
#' Transformation of the block diagonal parameters (true form)
#'
#' @param ParaValue Constrained parameter
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies  string-vector containing the names of the economies which are part of the economic system
#' @param GVARinputs Inputs used in the estimation of the GVAR-based models
#'
#' @keywords internal

True_BlockDiag <- function(ParaValue, Const_Type, FactorLabels, Economies, GVARinputs) {
  # Get the true parameter
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)

  step <- c(G * (G + 1) / 2, rep(K * (K + 1) / 2, times = C))

  # Input the zeros restrictions to the optimization vector (for the GVAR constrained models):
  if (any(GVARinputs$VARXtype == paste("constrained:", FactorLabels$Domestic))) {
    # Identify the index of the constrained variable
    zz <- nchar("constrained: ")
    VarInt <- substr(GVARinputs$VARXtype, start = zz + 1, stop = nchar(GVARinputs$VARXtype))
    IdxInt <- which(FactorLabels$Domestic == VarInt)

    # Identify the indexes of the zero restrictions within the optimization vector
    MatIntVector <- matrix(NA, nrow = K, ncol = K)
    MatIntVector[IdxInt, -IdxInt] <- 0
    MatIntVector[-IdxInt, IdxInt] <- 0

    IdxZeroRest <- MatIntVector[lower.tri(MatIntVector, diag = TRUE)] # indexes of the non-zero restrictions
    IdxNonZero <- which(is.na(IdxZeroRest)) # indexes of the zero restrictions

    # Initialize the complete vector (including the zero and the non-zero restrictions)
    stepRest <- (length(ParaValue) - step[1]) / C # number of parameters of the restricted vector
    ss <- rep(0, times = step[2]) # unrestricted vector of parameters per country
    tt <- rep(0, times = sum(step)) # complete unrestricted vector of parameters (per country + global)

    # Include the parameters of the marginal model
    idxI <- step[1]
    tt[seq_len(step[1])] <- ParaValue[seq_len(step[1])]

    # Include the parameters of the country-specific VARXs
    for (i in 1:C) {
      idxF <- idxI + stepRest
      ss[IdxNonZero] <- ParaValue[(idxI + 1):idxF]
      IdxFull <- sum(step[1:i])
      tt[(IdxFull + 1):(IdxFull + length(ss))] <- ss
      idxI <- idxF
    }

    ParaValue <- tt
  }

  idx0 <- 0
  b <- NULL

  for (j in 1:(C + 1)) {
    idx1 <- idx0 + step[j]

    if (idx0 < idx1) {
      seq_indices <- seq(idx0 + 1, idx1)
    } else {
      seq_indices <- integer(0)
    }
    d <- ParaValue[seq_indices]

    if (length(d) == 0) {
      btemp <- matrix(, nrow = 0, ncol = 0)
    } else {
      N <- floor(sqrt(2 * length(d)))
      k <- round(N * (N + 1) / 2)

      idx <- matrix(1:(N * N), c(N, N))
      MatrOne <- matrix(1, nrow = N, ncol = N)

      index1 <- t(t(idx[which(MatrOne == 1 & lower.tri(MatrOne, diag = TRUE))]))
      idx <- t(idx)
      index2 <- t(t(idx[which(MatrOne == 1 & lower.tri(MatrOne, diag = TRUE))]))
      Cvt_Mat <- matrix(0, nrow = N * N, ncol = N * (N + 1) / 2)
      Cvt_Mat[index1, ] <- diag(N * (N + 1) / 2)
      Cvt_Mat[index2, ] <- diag(N * (N + 1) / 2) # creates indexes to ensure that variance-covariance matrix is symmetric.

      btemp <- NULL
      m <- matrix(Cvt_Mat %*% d[1:k], N, N)
      btemp <- cbind(btemp, m %*% t(m)) # Recall that the maximization is made with parameters of SSP^(1/2)
    }
    if (j == 1) {
      BD_mat <- btemp
    } else {
      BD_mat <- adiag(BD_mat, btemp)
    }
    idx0 <- idx1
  }

  return(BD_mat)
}

############################################################################################################
#' Transformation of the JLL-related parameters (true form)
#'
#' @param ParaValue Constrained parameter value
#' @param Const_Type Type of constraint
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param JLLinputs Inputs used in the estimation of the JLL-based models
#'
#' @keywords internal

True_JLLstruct <- function(ParaValue, Const_Type, FactorLabels, Economies, JLLinputs) {
  # Rebuild the true parameter
  G <- length(FactorLabels$Global)
  K <- length(FactorLabels$Domestic)
  C <- length(Economies)
  Nspa <- length(FactorLabels$Spanned)
  Macro <- K - Nspa
  R <- C * K + G # Total number of factors of the model

  ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(Macro, Nspa, G, Economies, JLLinputs$DomUnit)$Sigma_Ye
  MatOnes <- matrix(1, nrow = R, ncol = R)
  MatOnes[ZeroIdxSigmaJLL] <- 0
  IdxNONzeroSigmaJLL <- which(MatOnes != 0 & MatOnes == 1 & lower.tri(MatOnes, diag = TRUE))

  # Redefine the vector of parameters to be maximized to include also the zeros restrictions
  abc <- matrix(0, R, R)
  abc[IdxNONzeroSigmaJLL] <- ParaValue # include the non-zero elements

  b <- abc %*% t(abc) # NOTE: b is not identical to SSZ, but this doesn't impact the optimization

  return(b)
}
