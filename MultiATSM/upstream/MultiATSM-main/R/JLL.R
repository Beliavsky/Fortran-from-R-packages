#' Estimates the P-dynamics from JLL-based models
#'
#' @param NonOrthoFactors numeric matrix (\code{K x Td}). Time series of risk factors before orthogonalization.
#' @param N positive integer. Number of country-specific spanned factors.
#' @param JLLinputs list. Necessary inputs to estimate JLL models:
#' \enumerate{
#'   \item \code{Economies}: character vector. Set of \code{C} economies in the system.
#'   \item \code{DomUnit}: character. Name of the dominant economy, or "None" if not assigned (for "JLL No DomUnit" model).
#'   \item \code{WishSigmas}: logical. TRUE to estimate variance-covariance matrices and Cholesky factorizations; FALSE otherwise.
#'   \item \code{SigmaNonOrtho}: NULL or F x F matrix from non-orthogonalized dynamics.
#'   \item \code{JLLModelType}: character. Permissible choices: "JLL original", "JLL joint Sigma", "JLL No DomUnit".
#' }
#' @param CheckInputs logical. Whether to perform a prior consistency check on the inputs provided in \code{JLLinputs}. Default is FALSE.
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td}: model time series dimension
#'   \item \code{C} number of countries in the system.
#'   \item \code{K}: total number of risk factors
#' }
#'
#' @examples
#' \donttest{
#' data(RiskFacFull)
#' RF_TS <- RiskFacFull
#' N <- 3
#' JLLinputs <- list(
#'   Economies = c("China", "Brazil", "Mexico", "Uruguay"), DomUnit = "China",
#'   WishSigmas = TRUE, SigmaNonOrtho = NULL, JLLModelType = "JLL original"
#' )
#' JLLPara <- JLL(RF_TS, N, JLLinputs)
#' }
#' @references
#' Jotiskhatira, P. ; Le, A. and Lundblad, C. (2015). "Why do interest rates in different currencies co-move?" (Journal of Financial Economics)
#' @return List of model parameters from both the orthogonalized and non-orthogonalized versions of the JLL-based models
#' @export

JLL <- function(NonOrthoFactors, N, JLLinputs, CheckInputs = FALSE) {
  # 0. Preliminary works/checks
  if (isTRUE(CheckInputs)) {
    CheckJLLinputs(NonOrthoFactors, JLLinputs)
  }

  # System dimension
  K <- nrow(NonOrthoFactors)
  G <- c() # Extract the number of global factors
  for (h in 1:K) {
    G[h] <- all(sapply(JLLinputs$Economies, grepl, rownames(NonOrthoFactors))[h, ] == 0)
  }
  G <- length(G[G == TRUE]) # Number of global unspnned factors
  C <- length(JLLinputs$Economies)

  # Factor labels
  Labels_All <- GetLabels_JLL(NonOrthoFactors, JLLinputs, G)

  # 1) Pre-allocation of the factors set
  Fact_NonOrtho <- Factors_NonOrtho(NonOrthoFactors, JLLinputs, Labels_All, N)

  # 2) Get coefficients from the orthogonalized regressions
  Ortho_RegSet <- OrthoReg_JLL(JLLinputs, N, Fact_NonOrtho, rownames(NonOrthoFactors), Labels_All$LabelsJLL)

  # 3) VAR(1) with orthogonalized factors
  Para_VAR_Ortho <- OrthoVAR_JLL(NonOrthoFactors, JLLinputs, Ortho_RegSet, Labels_All, N)

  # 4) Obtain the non-orthogonalized model parameters
  Para_VAR_NoOrtho <- NoOrthoVAR_JLL(Ortho_RegSet, Para_VAR_Ortho)

  # 5) Obtain sigmas/Cholesky factorizations
  Sigmas <- Get_Sigma_JLL(JLLinputs, Labels_All, Ortho_RegSet, Para_VAR_Ortho, N)

  # 6) Prepare the outputs
  outputs <- list(
    a_W = Ortho_RegSet$a_W, a_DU_CS = Ortho_RegSet$a_DU_CS, b = Ortho_RegSet$b, c = Ortho_RegSet$c,
    PIb = Ortho_RegSet$PIb, PIac = Ortho_RegSet$PIac, PI = Ortho_RegSet$PI, Ye = Para_VAR_Ortho$Ye,
    k0_e = Para_VAR_Ortho$k0_e, k1_e = Para_VAR_Ortho$k1_e, k0 = Para_VAR_NoOrtho$k0,
    k1 = Para_VAR_NoOrtho$k1, Sigmas = Sigmas
  )

  return(outputs)
}

###############################################################################################################
#' Find the indexes of zero-restrictions from the orthogonalized variance-covariance matrix from the JLL-based models
#'
#' @param M number of country-specific unspanned factors (scalar)
#' @param N number of country-specific spanned factors (scalar)
#' @param G number of global unspanned factors (scalar)
#' @param Economies Set of economies that are part of the economic system (string-vector)
#' @param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "None"
#'
#' @keywords internal
#' @return restricted version of the JLL of the Cholesky factorization (F x F)

IDXZeroRestrictionsJLLVarCovOrtho <- function(M, N, G, Economies, DomUnit) {
  C <- length(Economies)
  K <- (M + N) * C + G

  MatOnes <- matrix(1, nrow = K, ncol = K)
  MatOnes[upper.tri(MatOnes)] <- 0
  CholOrtho <- MatOnes

  if (DomUnit != "None") {
    IdxDomUnit <- which(DomUnit == Economies) # Index of the dominant country
  }

  # Zero restrictions of global variables on spanned factors
  idx0Global <- G + M
  for (h in 1:C) {
    idx1Global <- idx0Global + N
    CholOrtho[(idx0Global + 1):idx1Global, seq_len(G)] <- 0
    idx0Global <- idx1Global + M
  }

  # Zero restrictions of macro domestic variables on spanned factors
  for (i in 1:C) {
    idx0RowMacroSpanned <- G + M
    idx0ColMacroSpanned <- G + (i - 1) * (M + N)
    idx1ColMacroSpanned <- idx0ColMacroSpanned + M
    # a) For Dominant Unit
    if (DomUnit != "None") {
      for (h in 1:C) { # Fix the columns and loop through the rows
        idx1RowMacroSpanned <- idx0RowMacroSpanned + N
        CholOrtho[(idx0RowMacroSpanned + 1):idx1RowMacroSpanned, (idx0ColMacroSpanned + 1):idx1ColMacroSpanned] <- 0
        idx0RowMacroSpanned <- idx1RowMacroSpanned + M
      }
    } else {
      # b) For Non-dominant Unit
      idx1RowMacroSpanned <- idx0RowMacroSpanned + N
      CholOrtho[-((idx0ColMacroSpanned + 1):idx1ColMacroSpanned), (idx0ColMacroSpanned + 1):idx1ColMacroSpanned] <- 0
      idx0RowMacroSpanned <- idx1RowMacroSpanned + M
    }
  }

  # Zero restrictions of spanned factors on macro domestic variables
  for (i in 1:C) {
    idx0RowSpannedMacro <- G
    idx0ColSpannedMacro <- G + M + (i - 1) * (M + N)
    idx1ColSpannedMacro <- idx0ColSpannedMacro + N
    # a) For Dominant Unit
    if (DomUnit != "None") {
      for (h in 1:C) { # Fix the columns and loop through the rows
        idx1RowSpannedMacro <- idx0RowSpannedMacro + M
        CholOrtho[(idx0RowSpannedMacro + 1):idx1RowSpannedMacro, (idx0ColSpannedMacro + 1):idx1ColSpannedMacro] <- 0
        idx0RowSpannedMacro <- idx1RowSpannedMacro + N
      }
    } else {
      idx1RowSpannedMacro <- idx0RowSpannedMacro + M
      CholOrtho[-((idx0ColSpannedMacro + 1):idx1ColSpannedMacro), (idx0ColSpannedMacro + 1):idx1ColSpannedMacro] <- 0
      idx0RowSpannedMacro <- idx1RowSpannedMacro + N
    }
  }

  # Zero restrictions of Macro country i on Macro country j
  if (DomUnit != "None") {
    for (i in 1:C) {
      if (i != IdxDomUnit) {
        idx0RowMacroMacro <- G
        idx0ColMacroMacro <- G + (i - 1) * (M + N)
        idx1ColMacroMacro <- idx0ColMacroMacro + M
        for (h in 1:C) { # Fix the columns and loop through the rows
          idx1RowMacroMacro <- idx0RowMacroMacro + M
          if (i != h) {
            CholOrtho[(idx0RowMacroMacro + 1):idx1RowMacroMacro, (idx0ColMacroMacro + 1):idx1ColMacroMacro] <- 0
          }
          idx0RowMacroMacro <- idx1RowMacroMacro + N
        }
      }
    }
  }
  # Zero restrictions of Spanned factors of country i on Spanned factors country j
  if (DomUnit != "None") {
    for (i in 1:C) {
      if (i != IdxDomUnit) {
        idx0RowSpannedSpanned <- G + M
        idx0ColSpannedSpanned <- G + M + (i - 1) * (M + N)
        idx1ColSpannedSpanned <- idx0ColSpannedSpanned + N
        for (h in 1:C) { # Fix the columns and loop through the rows
          idx1RowSpannedSpanned <- idx0RowSpannedSpanned + N
          if (i != h) {
            CholOrtho[(idx0RowSpannedSpanned + 1):idx1RowSpannedSpanned, (idx0ColSpannedSpanned + 1):idx1ColSpannedSpanned] <- 0
          }
          idx0RowSpannedSpanned <- idx1RowSpannedSpanned + M
        }
      }
    }
  }

  VarCovOrtho <- CholOrtho %*% t(CholOrtho)
  IdxZerosVarCovOrtho <- which(VarCovOrtho == 0)
  IdxZerosSigma_Ye <- which(CholOrtho == 0)

  IDXzerosJLL <- list(Sigma_Ye = IdxZerosSigma_Ye, VarCovOrtho = IdxZerosVarCovOrtho)

  return(IDXzerosJLL)
}

####################################################################################################################
#' Generate the variable labels of the JLL models
#'
#' @param NonOrthoFactors Risk factors before the orthogonalization (FxT)
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#' @param G number of global unspanned factors (scalar)
#'
#' @keywords internal

GetLabels_JLL <- function(NonOrthoFactors, JLLinputs, G) {
  C <- length(JLLinputs$Economies)

  FactorsJLL <- unlist(lapply(JLLinputs$Economies, function(economy) {
    paste(rownames(NonOrthoFactors)[grepl(economy, rownames(NonOrthoFactors))], "JLL")
  }))

  FactorLabels <- lapply(JLLinputs$Economies, function(economy) {
    rownames(NonOrthoFactors)[grepl(economy, rownames(NonOrthoFactors))]
  })
  names(FactorLabels) <- JLLinputs$Economies

  FactorLabels$Global <- rownames(NonOrthoFactors)[seq_len(G)]
  LabelsJLL <- c(FactorLabels$Global, FactorsJLL)

  return(list(FactorLabels = FactorLabels, LabelsJLL = LabelsJLL))
}

####################################################################################################################
#' Makes the pre-allocation of the factors set for JLL-based models
#'
#' @param NonOrthoFactors Risk factors before the orthogonalization (FxT)
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#' @param FactorLab Variable labels from JLL-based models
#' @param N number of country-specific spanned factors (scalar)
#'
#' @keywords internal

Factors_NonOrtho <- function(NonOrthoFactors, JLLinputs, FactorLab, N) {
  M <- length(FactorLab$FactorLabels[[1]]) - N

  # Domestic factors
  FullFactorsSet <- lapply(JLLinputs$Economies, function(economy) {
    LabelofInt <- FactorLab$FactorLabels[[economy]]
    list(
      Macro = NonOrthoFactors[LabelofInt, ][1:M, ],
      Pricing = NonOrthoFactors[LabelofInt, ][(M + 1):(N + M), ]
    )
  })
  names(FullFactorsSet) <- JLLinputs$Economies

  # Global factors
  Lab_global <- FactorLab$FactorLabels$Global
  G <- length(Lab_global)
  MacroGlobal <- NonOrthoFactors[Lab_global, , drop = FALSE]

  Out <- list(MacroGlobal = MacroGlobal, FullFactorsSet = FullFactorsSet)
  return(Out)
}

####################################################################################################################
#' Compute Sigmas/Cholesky factorizations
#'
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#' @param FacSet Set of factors used in the estimation of JLL-based setups
#' @param Para_Ortho_Reg Set of parameters obtained from the JLL regressions
#' @param Para_Ortho_VAR Set of parameters obtained from the VAR(1) of JLL-based models
#' @param N number of country-specific spanned factors (scalar)
#'
#' @keywords internal

Get_Sigma_JLL <- function(JLLinputs, FacSet, Para_Ortho_Reg, Para_Ortho_VAR, N) {
  M <- length(FacSet$FactorLabels[[1]]) - N
  G <- length(FacSet$FactorLabels$Global)

  Ye <- Para_Ortho_VAR$Ye
  k0_e <- Para_Ortho_VAR$k0_e
  k1_e <- Para_Ortho_VAR$k1_e
  PI <- Para_Ortho_Reg$PI

  LabelsJLL <- FacSet$LabelsJLL

  if (JLLinputs$WishSigmas == 1) {
    # If the Variance-covariance matrix of the orthogonalized factors are NOT provided
    if (is.null(JLLinputs$SigmaNonOrtho)) {
      T_dim <- ncol(Ye)

      LHS <- Ye[, 2:T_dim]
      RHS <- Ye[, 1:(T_dim - 1)]

      et <- LHS - k0_e - k1_e %*% RHS
      SIGMA_Unres <- crossprod(t(et)) / dim(et)[2]
      # Labels
      dimnames(SIGMA_Unres) <- list(LabelsJLL, LabelsJLL)

      # If the estimation of SIGMA_Ye is necessary
      Sigma_Ye <- EstimationSigma_Ye(SIGMA_Unres, et, M, G, JLLinputs$Economies, JLLinputs$DomUnit)

      # Cholesky term (non-orthogonalized factors)
      Sigma_Y <- PI %*% Sigma_Ye

      # Variance-covariance matrices
      Sigma_Res_Ortho <- Sigma_Ye %*% t(Sigma_Ye) # Orthogonalized dynamics
      Sigma_Res_NonOrtho <- Sigma_Y %*% t(Sigma_Y) # Non-orthogonalized dynamics
    } else {
      Sigma_Res_NonOrtho <- JLLinputs$SigmaNonOrtho
      Sigma_Y <- t(chol(JLLinputs$SigmaNonOrtho))
      Sigma_Ye <- solve(PI) %*% Sigma_Y
      Sigma_Res_Ortho <- Sigma_Ye %*% t(Sigma_Ye)
    }

    ZeroIdxSigmaJLL <- IDXZeroRestrictionsJLLVarCovOrtho(M, N, G, JLLinputs$Economies, JLLinputs$DomUnit) # Identify the zero elements of the orthogonalized variance-covariance matrix
    # (useful for distinguishing real zeros from nearly zero elements later on in the code)
    Sigmas <- list(Sigma_Res_Ortho, Sigma_Res_NonOrtho, Sigma_Y, Sigma_Ye, ZeroIdxSigmaJLL)
    names(Sigmas) <- c("VarCov_Ortho", "VarCov_NonOrtho", "Sigma_Y", "Sigma_Ye", "ZeroIdxSigmaJLLOrtho")
  } else {
    Sigmas <- NULL
  }

  return(Sigmas)
}

#############################################################################################################
#' Check consistency of the inputs provided in JLL-based models
#'
#' @param RiskFactorsNonOrtho Risk factors before the orthogonalization (FxT)
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#'
#' @keywords internal

CheckJLLinputs <- function(RiskFactorsNonOrtho, JLLinputs) {
  # CHECK 1: Check whether the model type is correctly specified
  if (!JLLinputs$JLLModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    stop("JLLModelType input must be one of the following inputs: 'JLL original', 'JLL No DomUnit', 'JLL joint Sigma'.")
  }

  # CHECK 2: Check whether country names are correctly specified
  if (!all(sapply(JLLinputs$Economies, is.character))) {
    stop("All elements of the list 'Economies' must be exclusively country names")
  }

  # CHECK 3: Check for the consistency of dominant unit
  if ((grepl("JLL original", JLLinputs$JLLModelType) ||
    grepl("JLL jointSigma", JLLinputs$JLLModelType)) & JLLinputs$DomUnit == "None") {
    stop("In 'JLL original' and 'jointSigma', the  DomUnit input cannot be 'None'. One dominant country is required.")
  }

  if (grepl("JLL No DomUnit", JLLinputs$JLLModelType) & JLLinputs$DomUnit != "None") {
    stop("In 'JLL No DomUnit' DomUnit cannot cannot contain a name of a country. DomUnit must be set as 'None'.")
  }

  # CHECK 4: Check the exitence condition of global factors
  G <- c() # Extract the number of global factors
  for (h in 1:nrow(RiskFactorsNonOrtho)) {
    G[h] <- all(sapply(JLLinputs$Economies, grepl, rownames(RiskFactorsNonOrtho))[h, ] == 0)
  }
  G <- length(G[G == TRUE]) # Number of global unspnned factors

  if (G == 0) {
    stop("JLL-based models must contain at least one global factor.")
  }
}

###############################################################################################################
#' Get coefficients from the orthogonalized regressions
#'
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#' @param N number of country-specific spanned factors (scalar)
#' @param FacSet Set of factors used in the estimation of JLL-based setups
#' @param FactorLab_NonOrth Variable labels of the non-orthogonalized risk factors
#' @param FactorLab_JLL Variable labels of the orthogonalized risk factors
#'
#' @keywords internal

OrthoReg_JLL <- function(JLLinputs, N, FacSet, FactorLab_NonOrth, FactorLab_JLL) {
  # Preliminary work
  FullFactorsSet <- FacSet$FullFactorsSet
  MacroGlobal <- FacSet$MacroGlobal
  Label_DU <- JLLinputs$DomUnit # label of the dominant unit

  if (Label_DU != "None") {
    IdxDomUnit <- which(Label_DU == JLLinputs$Economies) # Index of the dominant country
  }

  Economies <- JLLinputs$Economies
  G <- nrow(FacSet$MacroGlobal)
  C <- length(Economies)
  K <- length(FactorLab_NonOrth)
  M <- (K - G) / C - N

  # 1) Orthogonalization of the pricing factors
  # Equation 6
  PricingRegressEQ6 <- lapply(JLLinputs$Economies, function(economy) {
    Pricing <- FullFactorsSet[[economy]]$Pricing
    Macro <- FullFactorsSet[[economy]]$Macro

    # Ensure Pricing is a matrix with the correct orientation
    PricingMat <- if (is.null(dim(Pricing))) {
      matrix(Pricing, ncol = length(Pricing)) # Convert to column vector if it's a numeric vector
    } else {
      Pricing
    }

    # Ensure Macro is a matrix with the correct orientation
    MacroMat <- if (is.null(dim(Macro))) {
      matrix(Macro, ncol = length(Macro)) # Convert to column vector if it's a numeric vector
    } else {
      Macro
    }

    stats::lm(t(PricingMat) ~ t(MacroMat) - 1)
  })

  b <- lapply(PricingRegressEQ6, function(model) t(model$coefficients))
  P_e <- lapply(PricingRegressEQ6, function(model) t(model$residuals))

  names(PricingRegressEQ6) <- Economies
  names(b) <- Economies
  names(P_e) <- Economies

  # Equation 10
  P_e_star <- list()
  if (Label_DU == "None") {
    c <- lapply(JLLinputs$Economies, function(economy) matrix(0, N, N))
    P_e_star <- lapply(JLLinputs$Economies, function(economy) NA)
  } else {
    PricingRegressEQ10 <- lapply(JLLinputs$Economies[-IdxDomUnit], function(economy) {
      stats::lm(t(P_e[[economy]]) ~ t(P_e[[Label_DU]]) - 1)
    })
    P_e_star <- lapply(PricingRegressEQ10, function(model) t(model$residuals))
    c <- lapply(PricingRegressEQ10, function(model) t(model$coefficients))
  }

  if (Label_DU != "None") {
    names(PricingRegressEQ10) <- Economies[-IdxDomUnit]
    names(c) <- Economies[-IdxDomUnit]
    names(P_e_star) <- Economies[-IdxDomUnit]
  }

  # 2) Orthogonalization of the macro factors
  MacroRegressEQ8 <- lapply(JLLinputs$Economies, function(economy) {
    if (is.null(dim(FullFactorsSet[[economy]]$Macro))) {
      Macro <- FullFactorsSet[[economy]]$Macro
      MacroMat <- matrix(Macro, ncol = length(Macro))
      stats::lm(t(MacroMat) ~ t(MacroGlobal) - 1)
    } else {
      stats::lm(t(FullFactorsSet[[economy]]$Macro) ~ t(MacroGlobal) - 1)
    }
  })

  a_W <- lapply(MacroRegressEQ8, function(model) t(model$coefficients))
  M_e <- lapply(MacroRegressEQ8, function(model) t(model$residuals))

  a_DU_CS <- lapply(JLLinputs$Economies, function(economy) matrix(0, M, G))
  M_e_CS <- lapply(JLLinputs$Economies, function(economy) NA)

  names(a_W) <- Economies
  names(M_e) <- Economies
  names(a_DU_CS) <- Economies
  names(M_e_CS) <- Economies

  if (Label_DU != "None") {
    MacroRegressEQ9 <- lapply(JLLinputs$Economies[-IdxDomUnit], function(economy) {
      if (is.null(dim(FullFactorsSet[[economy]]$Macro))) {
        Macro <- FullFactorsSet[[economy]]$Macro
        MacroMat <- matrix(Macro, ncol = length(Macro))
        stats::lm(t(MacroMat) ~ t(MacroGlobal) - 1)
      } else {
        stats::lm(t(FullFactorsSet[[economy]]$Macro) ~ t(MacroGlobal) + t(M_e[[Label_DU]]) - 1)
      }
    })

    a_W[JLLinputs$Economies[-IdxDomUnit]] <- lapply(MacroRegressEQ9, function(model) t(model$coefficients)[, seq_len(G)])
    if (M > 1) {
      a_DU_CS[JLLinputs$Economies[-IdxDomUnit]] <- lapply(MacroRegressEQ9, function(model) t(model$coefficients)[, (G + 1):(G + M)])
    } else {
      a_DU_CS[JLLinputs$Economies[-IdxDomUnit]] <- lapply(MacroRegressEQ9, function(model) model$coefficients)
    }
    M_e_CS[JLLinputs$Economies[-IdxDomUnit]] <- lapply(MacroRegressEQ9, function(model) t(model$residuals))
  }

  # Build the Pi matrices:
  PIb <- diag(K)
  idxRow0 <- G + M
  idxCol0 <- G
  for (i in 1:C) {
    Label_Eco <- JLLinputs$Economies[i] # label of all economies

    idxRow1 <- idxRow0 + N
    idxCol1 <- idxCol0 + M
    PIb[(idxRow0 + 1):idxRow1, (idxCol0 + 1):idxCol1] <- b[[Label_Eco]]
    idxRow0 <- idxRow1 + M
    idxCol0 <- idxCol1 + N
  }
  dimnames(PIb) <- list(FactorLab_NonOrth, FactorLab_JLL)
  # PIac
  PIac <- diag(K)
  idxRow00 <- G
  for (i in 1:C) {
    Label_Eco <- JLLinputs$Economies[i] # label of all economies

    idxRow11 <- idxRow00 + M
    PIac[(idxRow00 + 1):idxRow11, 1:G] <- a_W[[Label_Eco]]
    idxRow00 <- idxRow11 + N
  }

  if (Label_DU != "None") {
    # Place the orthogonalization of the pricing factors with respect to the dominant unit
    idxRowaDU0 <- G + M + N
    idxColaDU0 <- G
    idxColaDU1 <- idxColaDU0 + M
    for (j in 1:(C - 1)) {
      Label_No_DU <- JLLinputs$Economies[-IdxDomUnit][j] # label of the no dominant units

      idxRowaDU1 <- idxRowaDU0 + M
      PIac[(idxRowaDU0 + 1):idxRowaDU1, (idxColaDU0 + 1):idxColaDU1] <- a_DU_CS[[Label_No_DU]]
      idxRowaDU0 <- idxRowaDU1 + N
    }

    # Place the orthogonalization of the pricing factors with respect to the dominant unit
    # (c coefficients from equation 10)
    idxRowc0 <- G + M + N + M
    idxColc0 <- G + M
    idxColc1 <- idxColc0 + N
    for (j in 1:(C - 1)) {
      Label_No_DU <- JLLinputs$Economies[-IdxDomUnit][j] # label of the no dominant units

      idxRowc1 <- idxRowc0 + N
      PIac[(idxRowc0 + 1):idxRowc1, (idxColc0 + 1):idxColc1] <- c[[Label_No_DU]]
      idxRowc0 <- idxRowc1 + M
    }
  }

  dimnames(PIac) <- list(FactorLab_NonOrth, FactorLab_JLL)
  # PI
  PI <- PIb %*% PIac
  dimnames(PI) <- list(FactorLab_NonOrth, FactorLab_JLL)

  # 4) Output to export
  Times_Series <- list(P_e = P_e, P_e_star = P_e_star, M_e = M_e, M_e_CS = M_e_CS)
  Out <- list(a_W = a_W, a_DU_CS = a_DU_CS, b = b, c = c, PIb = PIb, PIac = PIac, PI = PI, TS_Factors = Times_Series)

  return(Out)
}

##############################################################################################################
#' VAR(1) with orthogonalized factors (JLL models)
#'
#' @param NonOrthoFactors Risk factors before the orthogonalization (FxT)
#' @param JLLinputs List of necessary inputs to estimate JLL-based setups
#' @param Ortho_Set Set of orthogonalized risk factors
#' @param FactLabels Variable labels of the orthogonalized risk factors
#' @param N number of country-specific spanned factors (scalar)
#'
#' @keywords internal

OrthoVAR_JLL <- function(NonOrthoFactors, JLLinputs, Ortho_Set, FactLabels, N) {
  # 1) Preliminary work
  LabelsJLL <- FactLabels$LabelsJLL
  Label_DU <- JLLinputs$DomUnit
  if (Label_DU != "None") {
    IdxDomUnit <- which(Label_DU == JLLinputs$Economies) # Index of the dominant country
  }

  K <- nrow(NonOrthoFactors)
  T_dim <- ncol(NonOrthoFactors)

  C <- length(JLLinputs$Economies)
  M <- length(FactLabels$FactorLabels[[1]]) - N
  G <- length(FactLabels$FactorLabels$Global)

  M_e <- Ortho_Set$TS_Factors$M_e
  P_e <- Ortho_Set$TS_Factors$P_e
  M_e_CS <- Ortho_Set$TS_Factors$M_e_CS
  P_e_star <- Ortho_Set$TS_Factors$P_e_star

  # 2) Build orthogonalized domestic risk factors
  AllDomFactorsOrtho <- do.call(rbind, lapply(JLLinputs$Economies, function(economy) {
    rbind(M_e[[economy]], P_e[[economy]])
  }))

  if (Label_DU != "None") {
    DomUnitOrtho <- rbind(M_e[[Label_DU]], P_e[[Label_DU]])
    NoDomUnitOrtho <- do.call(rbind, lapply(JLLinputs$Economies[-IdxDomUnit], function(economy) {
      rbind(M_e_CS[[economy]], P_e_star[[economy]])
    }))
    AllDomFactorsOrtho <- rbind(DomUnitOrtho, NoDomUnitOrtho)
  }

  # Add global factors
  MacroGlobal <- NonOrthoFactors[seq_len(G), ]
  Ye <- rbind(MacroGlobal, AllDomFactorsOrtho)
  rownames(Ye) <- LabelsJLL

  # Build the vector of non-orthogonalized factors
  if (Label_DU == "None") {
    Y <- NonOrthoFactors
  } else {
    Y <- matrix(NA, nrow = K, ncol = T_dim)
    rownames(Y) <- rownames(NonOrthoFactors)
    Y[seq_len(G), ] <- MacroGlobal # Global factors
    Y[(G + 1):(G + M + N), ] <- NonOrthoFactors[FactLabels$FactorLabels[[IdxDomUnit]], ] # Dominant country

    COUNTER0 <- G + M + N
    for (i in 1:(C - 1)) { # Non-dominant countries
      Eco_NoDU <- JLLinputs$Economies[-IdxDomUnit][i]

      COUNTER1 <- N + M + COUNTER0
      Y[(COUNTER0 + 1):COUNTER1, ] <- NonOrthoFactors[FactLabels$FactorLabels[[Eco_NoDU]], ]
      COUNTER0 <- COUNTER1
    }
  }

  # 3.a) Set the constraints on the feedback matrix
  Bcon_Mat <- FeedbackMatrixRestrictionsJLL(Label_DU, K, G, M, N)
  # 3.b) Estimate the VAR(1) with the orthogonalized variables
  intercept <- rep(1, times = T_dim - 1)
  RHS <- rbind(intercept, Ye[, 1:(T_dim - 1)])
  LHS <- Ye[, 2:T_dim]

  Coeff <- Est_RestOLS(LHS, RHS, Bcon_Mat)
  k0_e <- Coeff[, 1]
  k1_e <- Coeff[, 2:(K + 1)]

  # Add Labels to k1_e
  dimnames(k1_e) <- list(LabelsJLL, LabelsJLL)

  # Output to export
  Out <- list(k0_e = k0_e, k1_e = k1_e, Ye = Ye, Bcon = Bcon_Mat)
  return(Out)
}

#########################################################################################################
#' Obtain the non-orthogonalized model parameters
#'
#' @param Para_Ortho_Reg Set of parameters obtained from the JLL regressions
#' @param Para_Ortho_VAR Set of parameters obtained from the VAR(1) of JLL-based models
#'
#' @keywords internal

NoOrthoVAR_JLL <- function(Para_Ortho_Reg, Para_Ortho_VAR) {
  PI <- Para_Ortho_Reg$PI
  k0_e <- Para_Ortho_VAR$k0_e
  k1_e <- Para_Ortho_VAR$k1_e
  Bcon <- Para_Ortho_VAR$Bcon

  K <- nrow(PI)

  # Obtain the non-orthogonalized factors
  k0 <- PI %*% k0_e
  k1 <- PI %*% k1_e %*% solve(PI)
  # Ensures that the almost zero elements of k1 are actually zeros (this procedure avoids weird IRFs)
  idxZEROS <- which(Bcon[, 2:(K + 1)] == 0)
  k1[idxZEROS] <- 0

  Out <- list(k0 = k0, k1 = k1)
  return(Out)
}

#########################################################################################################
#' Build the log-likelihood function of the P-dynamics from the JLL-based models
#'
#' @param VecPara vector that contains all the non-zero entries of the lower-triangular part of the Cholesky factorization
#' @param res residuals from the VAR of the JLL model (K x T)
#' @param IdxNONzero vector that contains indexes of the matrix of the non-zero entries of the Cholesky factorization
#' @param K dimensions of the variance-covariance matrix (scalar)
#'
#' @keywords internal
#' @return value of the log-likelihood function (scalar)

llk_JLL_Sigma <- function(VecPara, res, IdxNONzero, K) {
  Se <- matrix(0, K, K)
  Se[IdxNONzero] <- VecPara # restricted Se matrix
  Sigma_Res <- Se %*% t(Se)
  y <- GaussianDensity(res, Sigma_Res)

  llk <- -mean(y)

  return(llk)
}

#################################################################################################################
#' Impose the zero-restrictions on the Cholesky-factorization from JLL-based models.
#'
#' @param SigmaUnres unrestricted variance-covariance matrix (K X K)
#' @param M number of domestic unspanned factors per country (scalar)
#' @param G number of global unspanned factors (scalar)
#' @param N number of country-specific spanned factors (scalar)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#'
#' @keywords internal
#' @return restricted version the Cholesky factorization matrix from JLL-based models (K x K)

CholRestrictionsJLL <- function(SigmaUnres, M, G, N, Economies, DomUnit) {
  C <- length(Economies)
  if (DomUnit != "None") {
    IdxDomUnit <- which(DomUnit == Economies) # Index of the dominant country
  }

  # Transform the matrix to be the Cholesky form
  SigmaUnres <- t(chol(SigmaUnres))

  # i) Zero restrictions of global variables on spanned factors
  idx0Global <- G + M
  for (h in 1:C) {
    idx1Global <- idx0Global + N
    SigmaUnres[(idx0Global + 1):idx1Global, 1:G] <- 0
    idx0Global <- idx1Global + M
  }

  # ii) Zero restrictions of macro domestic variables on spanned factors
  for (i in 1:C) {
    idx0RowMacroSpanned <- G + M
    idx0ColMacroSpanned <- G + (i - 1) * (M + N)
    idx1ColMacroSpanned <- idx0ColMacroSpanned + M
    # a) With a dominant unit
    if (DomUnit != "None") {
      for (h in 1:C) { # Fix the columns and loop through the rows
        idx1RowMacroSpanned <- idx0RowMacroSpanned + N
        SigmaUnres[(idx0RowMacroSpanned + 1):idx1RowMacroSpanned, (idx0ColMacroSpanned + 1):idx1ColMacroSpanned] <- 0
        idx0RowMacroSpanned <- idx1RowMacroSpanned + M
      }
    } else {
      # b) No dominant unit
      idx1RowMacroSpanned <- idx0RowMacroSpanned + N
      SigmaUnres[-((idx0ColMacroSpanned + 1):idx1ColMacroSpanned), (idx0ColMacroSpanned + 1):idx1ColMacroSpanned] <- 0
      idx0RowMacroSpanned <- idx1RowMacroSpanned + M
    }
  }

  # iii) Zero restrictions of spanned factors on macro domestic variables
  for (i in 1:C) {
    idx0RowSpannedMacro <- G
    idx0ColSpannedMacro <- G + M + (i - 1) * (M + N)
    idx1ColSpannedMacro <- idx0ColSpannedMacro + N
    # a) With a dominant unit
    if (DomUnit != "None") {
      for (h in 1:C) { # Fix the columns and loop through the rows
        idx1RowSpannedMacro <- idx0RowSpannedMacro + M
        SigmaUnres[(idx0RowSpannedMacro + 1):idx1RowSpannedMacro, (idx0ColSpannedMacro + 1):idx1ColSpannedMacro] <- 0
        idx0RowSpannedMacro <- idx1RowSpannedMacro + N
      }
    } else {
      # b) No dominant unit
      idx1RowSpannedMacro <- idx0RowSpannedMacro + M
      SigmaUnres[-((idx0ColSpannedMacro + 1):idx1ColSpannedMacro), (idx0ColSpannedMacro + 1):idx1ColSpannedMacro] <- 0
      idx0RowSpannedMacro <- idx1RowSpannedMacro + N
    }
  }
  # iv) Zero restrictions of Macro country i on Macro country j
  if (DomUnit != "None") {
    for (i in 1:C) {
      if (i != IdxDomUnit) {
        idx0RowMacroMacro <- G
        idx0ColMacroMacro <- G + (i - 1) * (M + N)
        idx1ColMacroMacro <- idx0ColMacroMacro + M
        for (h in 1:C) { # Fix the columns and loop through the rows
          idx1RowMacroMacro <- idx0RowMacroMacro + M
          if (i != h) {
            SigmaUnres[(idx0RowMacroMacro + 1):idx1RowMacroMacro, (idx0ColMacroMacro + 1):idx1ColMacroMacro] <- 0
          }
          idx0RowMacroMacro <- idx1RowMacroMacro + N
        }
      }
    }
  }

  # v) Zero restrictions of Spanned factors of country i on Spanned factors country j
  if (DomUnit != "None") {
    for (i in 1:C) {
      if (i != IdxDomUnit) {
        idx0RowSpannedSpanned <- G + M
        idx0ColSpannedSpanned <- G + M + (i - 1) * (M + N)
        idx1ColSpannedSpanned <- idx0ColSpannedSpanned + N
        for (h in 1:C) { # Fix the columns and loop through the rows
          idx1RowSpannedSpanned <- idx0RowSpannedSpanned + N
          if (i != h) {
            SigmaUnres[(idx0RowSpannedSpanned + 1):idx1RowSpannedSpanned, (idx0ColSpannedSpanned + 1):idx1ColSpannedSpanned] <- 0
          }
          idx0RowSpannedSpanned <- idx1RowSpannedSpanned + M
        }
      }
    }
  }

  CholJLL <- SigmaUnres

  return(CholJLL)
}

###########################################################################################################
#' Estimate numerically the Cholesky-factorization from the JLL-based models
#'
#' @param SigmaUnres unrestricted variance-covariance matrix (K x K)
#' @param res residuals from the VAR of the JLL model (K x T)
#' @param M number of domestic unspanned factors per country (scalar)
#' @param G number of global unspanned factors (scalar)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#'
#' @keywords internal
#' @return Cholesky-factorization after the maximization (K x K)

EstimationSigma_Ye <- function(SigmaUnres, res, M, G, Economies, DomUnit) {
  # SIGMA_Ye
  K <- nrow(SigmaUnres)
  C <- length(Economies)
  N <- (K - G - M * C) / C

  # Set the constraints in the Sigma matrix
  Se <- CholRestrictionsJLL(SigmaUnres, M, G, N, Economies, DomUnit)
  IdxNONzeroSigmaJLL <- which(Se != 0)
  x <- Se[IdxNONzeroSigmaJLL] # vector containing the initial guesses

  MLfunction <- function(...) llk_JLL_Sigma(..., res = res, IdxNONzero = IdxNONzeroSigmaJLL, K = K)

  res <- stats::optim(
    par = x,
    fn = function(par) ML_stable(x, MLfunction),
    method = "Nelder-Mead",
    control = list(
      maxit = 200000 * length(x),
      reltol = 1e-2,
      trace = 0
    )
  )

  Xmax <- res$par

  SIGMA_Ye <- matrix(0, K, K)
  SIGMA_Ye[IdxNONzeroSigmaJLL] <- Xmax # Cholesky term (orthogonalized factors)

  # Labels
  rownames(SIGMA_Ye) <- rownames(SigmaUnres)
  colnames(SIGMA_Ye) <- rownames(SigmaUnres)

  return(SIGMA_Ye)
}

#######################################################################################################
#' Set the zero-restrictions on the feedback matrix of JLL's P-dynamics
#'
#' @param DomUnit Name of the economy which is assigned as the dominant unit. \cr
#'               If no dominant unit is assigned, then this variable is defined as "none"
#' @param K Total number of risk factors of the economic system (scalar)
#' @param G Number of global unspanned factors (scalar)
#' @param M Number of country-specific unspanned factors (scalar)
#' @param N Number of country-specific spanned factors (scalar)
#'
#' @keywords internal
#' @return matrix containing the restrictions of the feedback matrix (K x K)

FeedbackMatrixRestrictionsJLL <- function(DomUnit, K, G, M, N) {
  C <- (K - G) / (M + N) # number of countries of the system

  Bcon <- matrix(0, nrow = K, ncol = K + 1) # Includes all variables and the intercept
  Bcon[, 1] <- NaN # Intercept
  Bcon[, 2:(G + 1)] <- NaN # Global

  if (DomUnit == "None") {
    IDXrow000 <- G
    IDXcol000 <- G + 1

    for (i in 1:C) {
      IDXrow111 <- IDXrow000 + M + N
      IDXcol111 <- IDXcol000 + M + N

      Bcon[(IDXrow000 + 1):IDXrow111, (IDXcol000 + 1):IDXcol111] <- NaN

      IDXrow000 <- IDXrow111
      IDXcol000 <- IDXcol111
    }
  } else { # With a dominant unit
    Bcon[, (G + 1 + 1):(G + M + N + 1)] <- NaN # Dominant unit

    IDXrow000 <- G + M + N
    IDXcol000 <- G + M + N + 1

    for (i in 1:(C - 1)) {
      IDXrow111 <- IDXrow000 + M + N
      IDXcol111 <- IDXcol000 + M + N

      Bcon[(IDXrow000 + 1):IDXrow111, (IDXcol000 + 1):IDXcol111] <- NaN

      IDXrow000 <- IDXrow111
      IDXcol000 <- IDXcol111
    }
  }

  return(Bcon)
}
