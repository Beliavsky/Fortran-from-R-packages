#' Estimates a GVAR(1) and VARX(1,1,1) models
#'
#' @param GVARinputs list. Inputs for GVAR model estimation:
#' \enumerate{
#'        \item  \code{Economies}:  character vector. Contains the \code{C} names of the economies included in the system.
#'        \item  \code{GVARFactors}: list. All variables used in the estimation of the VARX model \cr
#'                (see e.g. \code{GVARFactors} file for details);
#'        \item \code{VARXtype}: Permissible:
#'  \itemize{
#'        \item \code{'unconstrained'}: model is estimated without constraints (each equation is estimated individually by ordinary least square);
#'        \item \code{'constrained: Spanned Factors'}: The model is estimated with the restriction that foreign pricing factors do NOT
#'                                affect (i) domestic economic variables and (ii) domestic pricing factors (estimation via restricted least squares).
#'        \item \code{'constrained : [factor_name]'}: The model is estimated with the restriction that the specified risk factor
#'              is influenced only by its own lagged values and the lagged values of its corresponding star variables.
#'        (estimation via restricted least squares.)
#'          }
#'          \item \code{Wgvar}: The GVAR transition matrix (\code{C x C}) used in the model solution. \cr
#'                                (See the output from the \code{\link{Transition_Matrix}} function.).
#' }
#' @param N positive integer. Number of country-specific spanned factors.
#' @param CheckInputs logical. Whether to perform a prior consistency check on the inputs provided in \code{GVARinputs}. Default is FALSE.
#'
#' @return list. Contains:
#' \enumerate{
#'   \item parameters of the country-specific VARX(1,1,1):
#'     \itemize{
#'       \item intercept (M + N x 1)
#'       \item phi_1   (M + N x M + N)
#'       \item phi_1* (M + N x M + N)
#'       \item phi_g (M + N x M + N)
#'       \item Sigma (M + N x G)
#'     }
#'   \item parameters of the GVAR:
#'     \itemize{
#'       \item F0 (K x K)
#'       \item F1 (K x K)
#'       \item Sigma_y (K x K)
#'     }
#' }
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{C}: number of countries in the system
#'   \item \code{G}: number of global unspanned factors
#'   \item \code{M}: number of country-specific unspanned factors
#'   \item \code{N}: number of country-specific spanned factors
#'   \item \code{K}: total number of risk factors (K = C x (N + M) + G)
#' }
#'
#'
#' @examples
#' data(GVARFactors)
#'
#' GVARinputs <- list(
#'   Economies = c("China", "Brazil", "Mexico", "Uruguay"),
#'   GVARFactors = GVARFactors, VARXtype = "unconstrained"
#' )
#'
#' GVARinputs$Wgvar <- matrix(c(
#'   0, 0.83, 0.86, 0.38,
#'   0.65, 0, 0.13, 0.55,
#'   0.32, 0.12, 0, 0.07,
#'   0.03, 0.05, 0.01, 0
#' ), nrow = 4, ncol = 4)
#' N <- 3
#'
#' GVARPara <- GVAR(GVARinputs, N)
#'
#' @references
#' Chudik, A. and Pesaran, M. H. (2016). "Theory and Practice of GVAR modelling" (Journal of Economic Surveys)
#' @export

GVAR <- function(GVARinputs, N, CheckInputs = FALSE) {
  # 0) Check whether there are inconsistency in the specification of GVARinputs
  if (isTRUE(CheckInputs)) {
    CheckInputsGVAR(GVARinputs, N)
  }

  # 1) Preliminary work
  # Labels of group of variables
  DomAndStarLabels <- names(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors)
  L <- length(DomAndStarLabels)
  DomLabels <- DomAndStarLabels[1:(L / 2)]
  StarLabels <- DomAndStarLabels[(L / 2 + 1):L]
  GlobalLabels <- names(GVARinputs$GVARFactors$Global)

  # 2) Prepare variables to be used in the estimation
  Factors_Clean <- GVAR_PrepFactors(GVARinputs, DomLabels, StarLabels, GlobalLabels, N)

  # 3) Estimate the country-specific VARX(1,1,1) models
  ParaVARX <- VARX(GVARinputs, Factors_Clean, DomLabels, StarLabels, GlobalLabels, N)

  # 4) Estimate the dynamics of the global variables: VAR(1)
  ParaGlobal <- MarginalModelPara(GVARinputs)

  # 5) Build a GVAR(1)
  GVARoutputs <- BuildGVAR(ParaVARX, ParaGlobal, GVARinputs, DomLabels, GlobalLabels, N)

  return(GVARoutputs)
}

#####################################################################################################################
#####################################################################################################################
#' Estimate numerically the variance-covariance matrix from the GVAR-based models
#'
#' @param SigmaUnres Unrestricted variance-covariance matrix (K x K)
#' @param res residuals from the VAR of a GVAR model (K x T)
#' @param IdxVarRest index of the variable that is selected as strictly exogenous
#'
#' @keywords internal
#' @return  restricted version of the variance-covariance matrix a GVAR model (K x K)


EstimationSigma_GVARrest <- function(SigmaUnres, res, IdxVarRest) {
  # Choleski Factor
  Se <- t(chol(SigmaUnres))
  Se[IdxVarRest, -IdxVarRest] <- 0
  Se[-IdxVarRest, IdxVarRest] <- 0

  K <- nrow(SigmaUnres)

  # Set the constraints in the Sigma matrix
  IdxNONzeroGVAR <- which(Se != 0)
  x <- Se[IdxNONzeroGVAR] # vector containing the initial guesses


  MLfunction <- function(...) llk_JLL_Sigma(..., res = res, IdxNONzero = IdxNONzeroGVAR, K = K)

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
  SIGMA_Ye[IdxNONzeroGVAR] <- Xmax # Cholesky term (orthogonalized factors)
  SIGMA_Res <- SIGMA_Ye %*% t(SIGMA_Ye)

  # Labels
  rownames(SIGMA_Res) <- rownames(SigmaUnres)
  colnames(SIGMA_Res) <- rownames(SigmaUnres)

  return(SIGMA_Res)
}


################################################################################################################
#' Prepare risk factors for the estimation of the GVAR model
#'
#' @param GVARinputs List of inputs for GVAR-based models
#' @param DomLabels string-based vector containing label of the domestic risk factors
#' @param StarLabels string-based vector containing label of the star domestic risk factors
#' @param GlobalLabels string-based vector containing label of the global risk factors
#' @param N number of country-specific spanned factors (scalar)
#'
#'
#' @keywords internal


GVAR_PrepFactors <- function(GVARinputs, DomLabels, StarLabels, GlobalLabels, N) {
  T_dim <- length(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors[[1]])
  C <- length(GVARinputs$Economies)
  G <- length(GlobalLabels)
  M <- length(DomLabels) - N

  # 1) Prepare variables to be used in the estimation
  # a) Z.t:
  Z.t <- lapply(GVARinputs$Economies, function(economy) {
    X <- sapply(1:(M + N), function(j) GVARinputs$GVARFactors[[economy]]$Factors[[j]][2:T_dim])
    t(X)
  })
  names(Z.t) <- GVARinputs$Economies

  # b) Z lagged (Z.Lt)
  Z.Lt <- lapply(GVARinputs$Economies, function(economy) {
    X <- sapply(1:(M + N), function(j) GVARinputs$GVARFactors[[economy]]$Factors[[j]][1:(T_dim - 1)])
    t(X)
  })
  names(Z.Lt) <- GVARinputs$Economies

  # c) Z star lagged (Zstar.Lt)
  idx1 <- M + N
  Zstar.Lt <- lapply(GVARinputs$Economies, function(economy) {
    X <- sapply(1:(M + N), function(j) GVARinputs$GVARFactors[[economy]]$Factors[[idx1 + j]][1:(T_dim - 1)])
    t(X)
  })
  names(Zstar.Lt) <- GVARinputs$Economies

  # d) Global lagged (Global.Lt)
  Global.Lt <- do.call(rbind, lapply(seq_len(G), function(j) GVARinputs$GVARFactors$Global[[j]][1:(T_dim - 1)]))
  rownames(Global.Lt) <- GlobalLabels

  # Export outputs
  Outputs <- list(Z.t = Z.t, Z.Lt = Z.Lt, Zstar.Lt = Zstar.Lt, Global.Lt = Global.Lt)
  return(Outputs)
}

##############################################################################################################
#' Estimate a VARX(1,1,1)
#'
#' @param GVARinputs List of inputs for GVAR-based models
#' @param Factors_GVAR list containing the set of eisk factors used in the estimation of the VARX models
#' @param DomLabels string-based vector containing label of the domestic risk factors
#' @param StarLabels string-based vector containing label of the star domestic risk factors
#' @param GlobalLabels string-based vector containing label of the global risk factors
#' @param N number of country-specific spanned factors (scalar)
#'
#' @keywords internal


VARX <- function(GVARinputs, Factors_GVAR, DomLabels, StarLabels, GlobalLabels, N) {
  # Preliminary work
  Z.t <- Factors_GVAR$Z.t
  Z.Lt <- Factors_GVAR$Z.Lt
  Zstar.Lt <- Factors_GVAR$Zstar.Lt
  Global.Lt <- Factors_GVAR$Global.Lt

  C <- length(GVARinputs$Economies)
  T_dim <- length(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors[[1]])
  G <- length(GlobalLabels)
  M <- length(DomLabels) - N

  # Prepare regressor set in LHS and RHS
  LHS <- lapply(GVARinputs$Economies, function(economy) Z.t[[economy]])
  RHS <- lapply(GVARinputs$Economies, function(economy) {
    rbind(rep(1, times = T_dim - 1), Z.Lt[[economy]], Zstar.Lt[[economy]], Global.Lt)
  })
  names(LHS) <- GVARinputs$Economies
  names(RHS) <- GVARinputs$Economies

  # Set the foreign contemporaneous matrices to be equal to zero
  phi0_star <- lapply(GVARinputs$Economies, function(economy) matrix(0, nrow = M + N, ncol = M + N))
  names(phi0_star) <- GVARinputs$Economies

  # 1) Estimate the VARX(1,1,1)
  # a) unconstrained system:
  if (GVARinputs$VARXtype == "unconstrained") {
    ModelEstimates <- lapply(GVARinputs$Economies, function(economy) {
      stats::lm(t(LHS[[economy]]) ~ t(RHS[[economy]]) - 1)
    })
    Coeff <- lapply(ModelEstimates, function(model) t(model$coefficients))
    et <- lapply(ModelEstimates, function(model) t(model$residuals))
    Sigma <- lapply(et, function(residuals) crossprod(t(residuals)) / T_dim)
    names(Coeff) <- GVARinputs$Economies
    names(Sigma) <- GVARinputs$Economies

    #  b) constrained system: zero restrictions for the spanned factors in the feedback matrix
  } else if (GVARinputs$VARXtype == "constrained: Spanned Factors") {
    idxIntercept <- 1
    idxM <- idxIntercept + M
    idxP <- idxM + N
    idxM_star <- idxP + M
    idxP_star <- idxM_star + N
    idx_global <- idxP_star + G

    Bcon <- lapply(GVARinputs$Economies, function(economy) {
      Bcon <- matrix(NaN, nrow = nrow(LHS[[economy]]), ncol = nrow(RHS[[economy]]))
      Bcon[, (idxM_star + 1):idxP_star] <- 0
      Bcon
    })
    names(Bcon) <- GVARinputs$Economies

    Coeff <- lapply(GVARinputs$Economies, function(economy) {
      Est_RestOLS(LHS[[economy]], RHS[[economy]], Bcon[[economy]])
    })
    names(Coeff) <- GVARinputs$Economies
    et <- lapply(GVARinputs$Economies, function(economy) LHS[[economy]] - Coeff[[economy]] %*% RHS[[economy]])
    Sigma <- lapply(et, function(residuals) crossprod(t(residuals)) / T_dim)
    names(Sigma) <- GVARinputs$Economies

    # c) constrained system: one variable of the system is only affected by its own lags and the star counterparts
  } else if (any(GVARinputs$VARXtype == paste("constrained:", DomLabels))) {
    VARXLabs <- c("Intercept", DomLabels, StarLabels, GlobalLabels)
    zz <- nchar("constrained: ")
    VarInt <- substr(GVARinputs$VARXtype, start = zz + 1, stop = nchar(GVARinputs$VARXtype))

    idxIntercept <- 1
    idxCol <- which(grepl(VarInt, VARXLabs))
    idxRow <- which(grepl(VarInt, DomLabels))

    # c.1) Identify the zero-restrictions:
    Bcon <- lapply(GVARinputs$Economies, function(economy) {
      Bcon <- matrix(NaN, nrow = nrow(LHS[[economy]]), ncol = nrow(RHS[[economy]]))
      rownames(Bcon) <- DomLabels
      colnames(Bcon) <- VARXLabs
      Bcon[idxRow, -c(idxIntercept, idxCol)] <- 0
      Bcon
    })
    names(Bcon) <- GVARinputs$Economies

    Coeff <- lapply(GVARinputs$Economies, function(economy) {
      Est_RestOLS(LHS[[economy]], RHS[[economy]], Bcon[[economy]])
    })
    names(Coeff) <- GVARinputs$Economies
    et <- lapply(GVARinputs$Economies, function(economy) LHS[[economy]] - Coeff[[economy]] %*% RHS[[economy]])

    names(et) <- GVARinputs$Economies
    # Initial guess
    SigmaUnrest <- lapply(et, function(residuals) crossprod(t(residuals)) / T_dim)

    # Estimate sigma with restrictions
    Sigma <- lapply(GVARinputs$Economies, function(economy) {
      EstimationSigma_GVARrest(SigmaUnrest[[economy]], et[[economy]], idxRow)
    })
    names(Sigma) <- GVARinputs$Economies
  }

  # 2) Prepare outputs:
  idxPhi0 <- 1
  idxPhi1 <- idxPhi0 + (N + M)
  idxPhi1_star <- idxPhi1 + (N + M)
  idxPhi_global <- idxPhi1_star + G

  if (G == 0) {
    Idx_G <- c()
  } else {
    Idx_G <- (idxPhi1_star + 1):idxPhi_global
  }

  ParaVARX <- lapply(GVARinputs$Economies, function(economy) {
    list(
      Phi0 = Coeff[[economy]][, idxPhi0],
      Phi1 = Coeff[[economy]][, (idxPhi0 + 1):idxPhi1],
      Phi1_star = Coeff[[economy]][, (idxPhi1 + 1):idxPhi1_star],
      Phi_global = Coeff[[economy]][, Idx_G],
      Sigma = Sigma[[economy]],
      Phi0_star = phi0_star[[economy]]
    )
  })
  names(ParaVARX) <- GVARinputs$Economies

  return(ParaVARX)
}

##################################################################################################################
#' Build the GVAR(1) from the country-specific VARX(1,1,1)
#'
#' @param ParaVARX Set of VARX model parameters
#' @param GlobalPara Set of marginal model parameters
#' @param GVARinputs List of inputs for GVAR-based models
#' @param DomLabels string-based vector containing label of the domestic risk factors
#' @param GlobalLabels string-based vector containing label of the global risk factors
#' @param N number of country-specific spanned factors (scalar)
#'
#' @importFrom magic adiag
#'
#' @keywords internal


BuildGVAR <- function(ParaVARX, GlobalPara, GVARinputs, DomLabels, GlobalLabels, N) {
  C <- length(GVARinputs$Economies)
  T_dim <- nrow(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors[[1]])
  G <- length(GlobalLabels)
  M <- length(DomLabels) - N

  #  1) Build a0, Ai0, Ai1 and Wi
  # a) a0
  a0 <- Get_a0(GVARinputs, ParaVARX)
  # b) Ai0 and Ai1:
  Ai0 <- lapply(GVARinputs$Economies, function(economy) {
    cbind(diag(N + M), ParaVARX[[economy]]$Phi0_star)
  })
  Ai1 <- lapply(GVARinputs$Economies, function(economy) {
    cbind(ParaVARX[[economy]]$Phi1, ParaVARX[[economy]]$Phi1_star)
  })
  names(Ai0) <- GVARinputs$Economies
  names(Ai1) <- GVARinputs$Economies
  # c) Wi (link matrices)
  Wi <- BuildLinkMat(GVARinputs, N, M)

  # 2) Compute G0, G1 and Sigma
  Gs_Sigma <- Get_G0G1Sigma(ParaVARX, GVARinputs, Ai0, Ai1, Wi)

  # 3) Compute Gy.0,  Gy.1
  # a) Gy.0:
  Gy.0 <- adiag(diag(G), Gs_Sigma$G0)
  # b) Gy.1:
  Gy.1 <- Get_Gy1(ParaVARX, GVARinputs, Gs_Sigma$G1, GlobalPara$Phi_w1)

  # 4) Build the GVAR(1): y_t = F0 + F1* y_{t-1} + (Gy.0)^(-1)ey_t ( equation 19)
  # a) F0 and F1:
  Phi0VARX <- do.call(rbind, lapply(GVARinputs$Economies, function(econ) {
    matrix(ParaVARX[[econ]]$Phi0, ncol = 1)
  }))

  F0 <- solve(Gy.0) %*% rbind(GlobalPara$Phi_w0, Phi0VARX)
  F1 <- solve(Gy.0) %*% Gy.1

  # b) Sigma_y:
  Sigma_y <- adiag(GlobalPara$Sigma_w, Gs_Sigma$Sigma)

  # 5) Prepare labels of the tables
  labelsDomVar <- unlist(lapply(GVARinputs$Economies, function(econ) {
    paste(DomLabels, econ)
  }))

  labelsTables <- c(GlobalLabels, labelsDomVar)

  dimnames(F1) <- list(labelsTables, labelsTables)
  dimnames(Sigma_y) <- list(labelsTables, labelsTables)
  rownames(F0) <- labelsTables

  GVARoutputs <- list(VARX = ParaVARX, Gy.0 = Gy.0, F0 = F0, F1 = F1, Sigma_y = Sigma_y)
  return(GVARoutputs)
}
##########################################################################################################
#' Build country-specific link matrices
#'
#' @param GVARinputs  List of inputs for GVAR-based models
#' @param N number of country-specific spanned factors (scalar)
#' @param M number of country-specific unspanned factors (scalar)
#'
#' @keywords internal

BuildLinkMat <- function(GVARinputs, N, M) {
  C <- length(GVARinputs$Economies)

  bottomWi <- lapply(1:C, function(i) matrix(NA, nrow = N + M, ncol = C * (N + M)))

  # Bottom part of Wi
  for (j in 1:C) {
    count0 <- 0
    for (i in 1:C) {
      count1 <- count0 + (N + M)
      bottomWi[[j]][, (count0 + 1):count1] <- GVARinputs$Wgvar[j, i] * diag(N + M)
      count0 <- count1
    }
  }

  names(bottomWi) <- GVARinputs$Economies

  # Top part of Wi
  # Initialize topWi as a list of NA-filled matrices
  topWi <- lapply(1:C, function(i) matrix(NA, nrow = N + M, ncol = C * (N + M)))

  IndexPos <- diag(C)
  for (j in 1:C) {
    count0 <- 0
    for (i in 1:C) {
      count1 <- count0 + (N + M)
      topWi[[j]][, (count0 + 1):count1] <- IndexPos[j, i] * diag(N + M)
      count0 <- count1
    }
  }
  names(topWi) <- GVARinputs$Economies

  # Concatenate topWi and bottomWi into Wi
  Wi <- lapply(1:C, function(i) rbind(topWi[[i]], bottomWi[[i]]))
  names(Wi) <- GVARinputs$Economies

  return(Wi)
}

##############################################################################################################
#' Estimate the marginal model for the global factors
#'
#' @param GVARinputs List of inputs for GVAR-based models
#'
#' @keywords internal

MarginalModelPara <- function(GVARinputs) {
  T_dim <- nrow(GVARinputs$GVARFactors[[GVARinputs$Economies[1]]]$Factors[[1]])
  G <- length(GVARinputs$GVARFactors$Global)

  X <- do.call(rbind, lapply(GVARinputs$GVARFactors$Global, matrix, ncol = T_dim))
  if (length(X) != 0) {
    X <- t(X) # useful command for the case in which G <- 1
    RHS <- as.matrix(X[2:T_dim, ])
    LHS <- as.matrix(X[1:(T_dim - 1), ])

    Phi_w1 <- matrix(NA, nrow = G, ncol = G)
    Phi_w0 <- matrix(NA, nrow = G, ncol = 1)
    for (i in seq_len(G)) {
      Phi_w0[i, ] <- stats::lm(LHS[, i] ~ RHS)$coefficients[1]
      Phi_w1[i, ] <- stats::lm(LHS[, i] ~ RHS)$coefficients[seq_len(max(G, 0)) + 1]
    }

    eta_t <- t(LHS) - matrix(Phi_w0, ncol = T_dim - 1, nrow = G) - Phi_w1 %*% t(RHS)

    Sigma_w <- (eta_t %*% t(eta_t)) / T_dim
  } else {
    Phi_w0 <- c()
    Phi_w1 <- matrix(, ncol = 0, nrow = 0)
    Sigma_w <- matrix(, ncol = 0, nrow = 0)
  }

  ParaExport <- list(Phi_w0 = Phi_w0, Phi_w1 = Phi_w1, Sigma_w = Sigma_w)

  return(ParaExport)
}

##########################################################################################################
#' Get the intercept, feedback matrix and the variance-covariance matrix from GVAR without global factors
#'
#' @param ParaVARX Set of VARX model parameters
#' @param GVARinputs List of inputs for GVAR-based models
#' @param Ai0 list containing the country-specific intercepts
#' @param Ai1 list containing the country-specific feedback matrices
#' @param Wi list containing the country-specific link matrices
#'
#' @keywords internal


Get_G0G1Sigma <- function(ParaVARX, GVARinputs, Ai0, Ai1, Wi) {
  C <- length(Ai1)
  MN <- nrow(Ai1[[1]])

  # a) G0 and G1
  G0 <- do.call(rbind, lapply(1:C, function(i) as.matrix(Ai0[[i]] %*% Wi[[i]], ncol = C * MN)))
  G1 <- do.call(rbind, lapply(1:C, function(i) as.matrix(Ai1[[i]] %*% Wi[[i]], ncol = C * MN)))

  # b) Sigma
  Sigma <- matrix(0, ncol = C * (MN), nrow = C * (MN))
  count0 <- 0
  for (i in 1:C) {
    count1 <- count0 + (MN)
    Sigma[(count0 + 1):count1, (count0 + 1):count1] <- ParaVARX[[GVARinputs$Economies[i]]]$Sigma
    count0 <- count1
  }

  out <- list(G0 = G0, G1 = G1, Sigma = Sigma)

  return(out)
}
############################################################################################################
#' Obtain the country-specific a0
#'
#' @param GVARinputs List of inputs for GVAR-based models
#' @param ParaVARX List containing the set of VARX model parameters
#'
#' @keywords internal


Get_a0 <- function(GVARinputs, ParaVARX) {
  C <- length(ParaVARX)
  MN <- length(ParaVARX[[1]]$Phi0)

  a0 <- matrix(NA, nrow = C * (MN), ncol = 1)

  count0 <- 0
  for (j in 1:C) {
    count1 <- count0 + MN
    a0[(count0 + 1):count1] <- do.call(rbind, lapply(ParaVARX[[GVARinputs$Economies[j]]]$Phi0, matrix, ncol = 1))
    count0 <- count1
  }

  return(a0)
}


########################################################################################################
#' Compute the feedback matrix from a GVAR model with global factors
#'
#' @param ParaVARX List containing the set of VARX model parameters
#' @param GVARinputs  List of inputs for GVAR-based models
#' @param G1 feedback matrix from a GVAR without global variables
#' @param Phi_w1 feedback matrix from a marginal model
#'
#' @keywords internal

Get_Gy1 <- function(ParaVARX, GVARinputs, G1, Phi_w1) {
  G <- length(GVARinputs$GVARFactors$Global)
  C <- length(GVARinputs$Economies)
  MN <- length(ParaVARX[[1]]$Phi0)

  if (G != 0) {
    D1 <- list()
    for (i in 1:C) {
      D1[[i]] <- ParaVARX[[GVARinputs$Economies[[i]]]]$Phi_global
    }
    D1 <- do.call(rbind, lapply(D1, matrix, ncol = G))
  } else {
    D1 <- c()
  }

  topGy.1 <- matrix(0, nrow = G, ncol = C * (MN) + G)
  topGy.1[seq_len(G), seq_len(G)] <- Phi_w1

  bottomGy.1 <- cbind(D1, G1)
  Gy.1 <- rbind(topGy.1, bottomGy.1)

  return(Gy.1)
}
#############################################################################################################
#' Check consistency of the inputs provided in GVARinputs
#'
#' @param GVARinputs List of inputs for GVAR-based models
#' @param N number of country-specific spanned factors (scalar)
#'
#' @keywords internal


CheckInputsGVAR <- function(GVARinputs, N) {
  # CHECK 1: correctly specified number of spanned factors
  if (!(N %in% 1:8)) {
    stop("N, the number of country-specific spanned factors, must be an integer between 1 and 8.")
  }

  # CHECK 2: Check the consistency of the names of the lists in GVARinputs
  if (!all(c("Economies", "GVARFactors", "VARXtype", "Wgvar") %in% names(GVARinputs))) {
    stop("The list elements of GVARinputs must be named 'Economies', 'GVARFactors', 'VARXtype', 'Wgvar'")
  }

  # CHECK 3: Check whether country names are correctly specified
  if (!(all(sapply(GVARinputs$Economies, is.character)))) {
    stop("All elements of the list 'Economies' must be exclusively country names")
  }

  # CHECK 4: Check for the consistency of VARXtype
  if (!(grepl("^constrained", GVARinputs$VARXtype) || GVARinputs$VARXtype == "unconstrained")) {
    stop("GVARinputs$VARXtype must be 'unconstrained'or'constrained'.")
  }
}
