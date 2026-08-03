#' Compute the maximum likelihood function of all models
#'
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param SSZ variance-covariance matrix (F x F)
#' @param K0Z intercept from the P-dynamics (F x 1)
#' @param K1Z feedback matrix from the P-dynamics (F x F)
#' @param se Variance of the portfolio of yields observed with error (scalar). Default is set to NULL.
#' @param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param P complete set of spanned factors (N x T or CN x T)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param WpcaFull composite matrix of weights the portfolios observed with and without errors
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param GVARinputs if the model chosen is the "GVAR single" or "GVAR multi", the "GVARinputs" should be specified (see "GVAR" function)
#' @param JLLinputs if the model chosen is JLL-based. "JLLinputs" should contain (i) DomUnit, (ii) WishSigmas, (iii) SigmaNonOrtho, (iv) JLLModelType (See JLL function)
#' @param BS_outputs Generates simplified output list in the bootstrap setting. Default is set to FALSE.
#' @param ExportListOut export the complete ATSM outputs. Default is TRUE.
#'
#' @references
#' \itemize{
#'  \item Candelon, C. and Moura, R. (2024). “A Multicountry Model of the Term Structures of Interest Rates with a GVAR.”
#'  Journal of Financial Econometrics 22 (5): 1558–87.
#'  \item Jotikasthira, C; Le, A. and Lundblad, C (2015). “Why Do Term Structures in Different Currencies Co-Move?”
#'  Journal of Financial Economics 115: 58–83.
#'  \item Joslin, S,; Priebsch, M. and Singleton, K. (2014). “Risk Premiums in Dynamic Term Structure Models with Unspanned Macro Risks.”
#'  Journal of Finance 69 (3): 1197–1233.
#'  \item Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on Gaussian dynamic term structure models".
#'  The Review of Financial Studies.
#'  \item Le, A. and Singleton, K. (2018). "A Small Package of Matlab Routines for the Estimation of Some Term Structure Models."
#'  Euro Area Business Cycle Network Training School - Term Structure Modelling.
#' }
#'
#' @keywords internal

MLEdensity <- function(K1XQ, r0, SSZ, K0Z, K1Z, se, Gy.0, mat, Y, Z, P, Wpca, We, WpcaFull, dt, Economies,
                       FactorLabels, ModelType, GVARinputs = NULL, JLLinputs = NULL, BS_outputs = FALSE,
                       ExportListOut = TRUE) {
  # 0) Initialize some variables
  N <- length(FactorLabels$Spanned)
  if (ModelType == "JLL joint Sigma") SSZ <- Update_SSZ_JLL(SSZ, Z, N, JLLinputs)

  # 1) Compute loadings:
  # Slope coefficients (Bn's)
  LoadBs <- Get_Bs(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economies, ModelType)
  # Intercepts (An's):
  if (is.null(r0) || length(r0) == 0) r0 <- Get_r0(Y, P, N, mat, dt, LoadBs, Wpca, We, Economies, ModelType)
  LoadAs <- Get_As(LoadBs, Wpca, r0, dt, Economies, ModelType)

  # 2) Build Log-likelihood density function
  llkOut <- Get_llk(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType)
  y <- llkOut$y
  se <- llkOut$se

  # 3) Variance-yields
  VarYields <- Get_SigmaYields(Y, N, mat, WpcaFull, se, ModelType)

  # 4) Output to export
  if (ExportListOut) {
    out <- OptOutputs(
      Y, Z, mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
      y, GVARinputs, JLLinputs, Economies, ModelType, BS_outputs
    )
  }

  if (!ExportListOut) {
    return(y)
  } else {
    return(out)
  }
}
##############################################################################################################
#' Extract the indexes related to the spanned factors in the variance-covariance matrix
#'
#' @param G number of global unspanned factors (scalar)
#' @param M number of domestic unspanned factors per country (scalar)
#' @param N number of domestic spanned factors per country (scalar)
#' @param C number of countries of the economic system (scalar)
#'
#' @keywords internal

IdxSpanned <- function(G, M, N, C) {
  K <- C * (N + M) + G
  vector <- 1:K

  idxA <- seq(from = G + M, by = N + M, length.out = C)
  IDX_F <- unlist(lapply(1:C, function(i) vector[(idxA[i] + 1):(idxA[i] + N)]))

  return(IDX_F)
}
##############################################################################################################
#' Compute long-run risk neutral mean (r0) for the various models
#'
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param P complete set of spanned factors (N x T or CN x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#' @param B_list list containing the B loadings
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_r0 <- function(Y, P, N, mat, dt, B_list, Wpca, We, Economies, ModelType) {
  # General procedure:
  # (i) A0 = (I - Bx(W*Bx)^(-1)*W)*Ax0;
  # (ii) A1 = (I - Bx(W*Bx)^(-1)*W)*Axr;
  # (iii) APer= We * A1
  # (iv) vt = We*(Yt - Bx(W*Bx)^(-1)*Pt - A0)
  # Recall that BnP= Bx(W*Bx)^(-1)

  J <- length(mat) # number of country-specific yields used in estimation;
  T_dim <- dim(Y)[2]
  tt <- 2:T_dim

  BnP <- B_list$BnP
  B_adj <- B_list$B_adj

  # 1) Get r0 for the models estimated on a country-by-country basis
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    A0 <- (diag(J) - BnP %*% Wpca) %*% B_adj
    A1 <- (diag(J) - BnP %*% Wpca) %*% matrix(1, J, 1) / dt

    #  r0 = APer'*vt/(APer'*APer)
    A_per <- t(A1) %*% (t(We) %*% We)
    num <- A_per %*% A1 # Numerator from equation r0 equation
    den <- A_per %*% (rowMeans(Y[, tt] - BnP %*% P[, tt]) - A0) # Denominator from r0 equation

    check_numeric(num, "num")
    check_numeric(den, "den")

    r0 <- solve(num, den)
  } else {
    # 2) Get r0 for the models estimated on a jointly basis
    C <- length(Economies)
    r0 <- numeric(C)

    idxA <- 0
    idxB <- 0
    idxC <- 0

    for (i in 1:C) {
      idxAA <- idxA + N
      idxBB <- idxB + J
      idxCC <- idxC + J - N

      BnPCS <- matrix(BnP[(idxB + 1):idxBB, (idxA + 1):idxAA], nrow = J)
      WpcaCS <- Wpca[(idxA + 1):idxAA, (idxB + 1):idxBB]
      B_adjCS <- B_adj[(idxB + 1):idxBB]
      WeCS <- We[(idxC + 1):idxCC, (idxB + 1):idxBB]
      YCS <- Y[(idxB + 1):idxBB, ]
      PCS <- P[(idxA + 1):idxAA, ]

      A0 <- (diag(J) - BnPCS %*% WpcaCS) %*% B_adjCS
      A1 <- (diag(J) - BnPCS %*% WpcaCS) %*% matrix(1, J, 1) / dt

      # r0 = APer'*vt/(APer'*APer)
      A_per <- t(A1) %*% (t(WeCS) %*% WeCS)
      PCS_adj <- if (N == 1) PCS[tt] else PCS[, tt]

      num <- A_per %*% A1 # Numerator from r0 equation
      den <- A_per %*% (rowMeans(YCS[, tt] - BnPCS %*% PCS_adj) - A0) # Denominator from r0 equation

      check_numeric(num, "num")
      check_numeric(den, "den")

      r0[i] <- solve(num, den)

      check_numeric(r0, "r0")

      idxA <- idxAA
      idxB <- idxBB
      idxC <- idxCC
    }
  }

  return(r0)
}
##############################################################################################################
#' Build the B loadings
#'
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param SSZ variance-covariance matrix (F x F)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economy string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_Bs <- function(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economy, ModelType) {
  Lab_SingleQ <- c("JPS original", "JPS global", "GVAR single")

  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N

  # 1) Get BnX and BnP
  # Yields are affine function of the vector Pt,i.e., Y(t) = AnP +BnP*P(t).
  # Further, we define Z(t) as an affine function of P(t) such that: Z(t) = phi0+ phi1*P(t).
  # As such, we can write P(t) = phi1^(-1)*(Z(t) - phi0).
  BnX_Unadj <- Get__BnXAnX(round(mat / dt), K1XQ, ModelType, Economies = Economy)$BnX
  BnX <- BnX_Unadj / dt
  check_numeric(BnX, "BnX")

  # We multiply by "1/dt" to produce annualized results.

  if (any(ModelType == Lab_SingleQ & N == 1)) {
    BnP <- matrix(BnX) %*% InvMat_Robust(Wpca %*% BnX)
  } else {
    BnP <- BnX %*% InvMat_Robust(Wpca %*% BnX)
  }
  check_numeric(BnP, "BnP")

  # 2) Covariance matrix of latent factors X:
  if (any(ModelType == Lab_SingleQ)) {
    AllLabels <- GetLabels_sepQ(Economy, ModelType, FactorLabels)
    dimnames(SSZ) <- list(AllLabels, AllLabels)
    LabelSpannedCS <- c(FactorLabels$Tables[[Economy]][-(1:M)])
    idxSpa <- match(LabelSpannedCS, AllLabels)
  } else {
    G <- length(FactorLabels$Global)
    C <- length(Economy)
    idxSpa <- IdxSpanned(G, M, N, C)
  }

  WBX <- Wpca %*% BnX # WBX = W * BnX
  SSP <- SSZ[idxSpa, idxSpa] # SSP -> Covariance matrix from observed (spanned) factors

  if (is.null(dim(SSP))) {
    SSX <- matrix(Inf, nrow = 1, ncol = 1) # Initialize with Inf
  } else {
    SSX <- matrix(Inf, nrow = nrow(SSP), ncol = ncol(SSP)) # Try to assign the correct value
  }

  SSX <- safe_solve(WBX, SSP) %*% safe_solve(t(WBX))
  check_numeric(SSX, "SSX")

  # 3) Optimal estimate of r0. NOTE: We set r0=0, because when r0=0, An = B_adj.
  if (any(ModelType == Lab_SingleQ)) {
    r0 <- 0
  } else {
    r0 <- rep(0, times = C)
  }

  B_adj <- Get__BnXAnX(round(mat / dt), K1XQ, ModelType, r0, SSX, Economies = Economy)$B_adj
  B_adj <- t(t(B_adj / dt))
  check_numeric(B_adj, "B_adj")

  out <- list(BnP = BnP, BnX = BnX, B_adj = B_adj, SSP = SSP, SSX = SSX, idxSpa = idxSpa)

  return(out)
}
##############################################################################################################
#' Compute the A loadings
#'
#' @param LoadBs list containing the B loadings
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param dt time interval unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_As <- function(LoadBs, Wpca, r0, dt, Economies, ModelType) {
  B_adj <- LoadBs$B_adj
  BnP <- LoadBs$BnP

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    AnX <- B_adj + as.numeric(r0 / dt)
  } else {
    C <- length(Economies)
    J <- nrow(B_adj) / C

    Q <- do.call(rbind, lapply(1:C, function(i) matrix(rep(r0[i], times = J) / dt, nrow = J, ncol = 1)))
    AnX <- B_adj + Q
  }

  AnP <- AnX - BnP %*% (Wpca %*% AnX) # AnP = AnX - BnX*(W*BnX)^(-1)*W*AnX

  return(list(AnX = AnX, AnP = AnP))
}
##############################################################################################################
#' Generate the factor labels for models estimated on a country-by-country basis
#'
#' @param Economy string containing the names of the economy to be estimated
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param FactorLabels list containing the factor labels
#'
#' @keywords internal

GetLabels_sepQ <- function(Economy, ModelType, FactorLabels) {
  if (ModelType == "JPS original") {
    AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", "GVAR single"))) {
    AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  return(AllLabels)
}
##############################################################################################################
#' Compute the log-likelihood function
#'
#' @param P time-series of spanned factors (N x T or CN x T)
#' @param Y time-series of yields (J x T or CJ x T)
#' @param Z time-series of risk factors (F x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param We matrix of weights of the portfolios observed with errors ((J-N) x J or C(J-N) x CJ)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x CJ)
#' @param K0Z matrix of intercepts (P-dynamics)
#' @param K1Z feedback matrix (P-dynamics)
#' @param SSZ variance-covariance matrix (P-dynamics)
#' @param LoadBs list containing the B loadings
#' @param LoadAs list containing the A loadings
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_llk <- function(P, Y, Z, N, mat, We, Wpca, K0Z, K1Z, SSZ, LoadBs, LoadAs, ModelType) {
  T_dim <- dim(P)[2]
  t <- 2:T_dim
  J <- length(mat)

  Peo <- We %*% Y # portfolio observed WITH errors
  BnP <- LoadBs$BnP
  AnP <- LoadAs$AnP

  # 1) density of yields pricing errors:
  MatOnes <- matrix(1, ncol = ncol(P), nrow = 1)
  Pe <- We %*% (BnP %*% P + AnP %*% MatOnes)

  eQ <- Peo[, t] - Pe[, t]

  # 2) Standard deviation of the measurement error.
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    se <- sqrt(mean(as.vector(eQ^2)))
  } else {
    C <- nrow(P) / N
    se <- numeric(C)
    idx0 <- 0
    for (i in 1:C) {
      idx1 <- idx0 + J - N
      se[i] <- sqrt(mean(as.vector(eQ[(idx0 + 1):idx1, ]^2))) # Standard deviation of the measurement error (scalar).
      idx0 <- idx1
    }
  }

  # 3) The log-likelihood function:
  # Cross-sectional density (i.e. density for the portfolios observed with measurement error)
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    y <- GaussianDensity(eQ, se^2 * diag(J - N))
  } else {
    aa <- se
    idx0 <- 0
    for (h in 1:C) {
      idx1 <- idx0 + J - N
      se[(idx0 + 1):idx1] <- rep(aa[h], times = J - N)
      idx0 <- idx1
    }

    se <- se[seq(1, C * J - C * N, by = J - N)] # Recast se
    y <- GaussianDensity(eQ, se^2 * diag(C * J - C * N))
  }



  # Time series density:
  MatOne_Z <- matrix(1, nrow = 1, ncol = ncol(Z) - 1)
  eP <- Z[, -1] - K1Z %*% Z[, 1:(T_dim - 1)] - K0Z %*% MatOne_Z

  y <- y + GaussianDensity(eP, SSZ) # Cross-sectional density + time-series density (final likelihood function, except for the Jacobian term)

  # Jacobian: (so that the density is for the observed yields (and the non-yields variables Z)
  #           and not the yields portfolios which depend on W).
  y <- y + 0.5 * log(abs(det(Wpca %*% t(Wpca))))
  y <- -t(as.vector(y)) # necessary for minimization

  return(list(y = y, se = se))
}
##############################################################################################################
#' computes the density function of a gaussian process
#'
#' @param res matrix of residuals (N x T)
#' @param SSZ covariance matrices or array of covariance matrices (N x N)
#'
#' @keywords internal
#'
#' @return y vector of density (1 x T)

GaussianDensity <- function(res, SSZ) {
  N <- nrow(res)

  check_numeric(SSZ, "SSZ")

  SSres <- chol2inv(chol(SSZ)) %*% res
  logabsdetSSZ <- 0.5 * (2 * sum(log(abs(diag(chol(SSZ %*% t(SSZ)))))))

  y <- -0.5 * N * log(2 * pi) - 0.5 * logabsdetSSZ - 0.5 * abs(colSums(res * SSres))

  return(y)
}
##############################################################################################################
#' Compute the variance-covariance matrix of the bond yields
#'
#' @param YieldsTS matrix of yields used in estimation (J x T or CJ x T)
#' @param N number of country-specific spanned factors
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param WpcaFull composite matrix of weights the portfolios observed with and without errors
#' @param se Variance of the portfolio of yields observed with error (scalar). Default is set to NULL
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

Get_SigmaYields <- function(YieldsTS, N, mat, WpcaFull, se, ModelType) {
  J <- length(mat)

  # Single-country models
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Sigma_e_Row <- c(rep(0, times = N), rep(se, times = J - N)^2) # Variances of portfolio WITHOUT and WITH errors
    Sigma_e <- Sigma_e_Row * diag(J)
    SIGMA_yields <- solve(WpcaFull) %*% Sigma_e %*% t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)
  } else {
    # Multicountry models
    C <- length(se)

    Sigma_e_AllRows <- rep(0, times = N * C) # Pre-allocate the vector with N*C zeros
    for (i in 1:C) {
      Sigma_e_CS <- c(rep(0, N), rep(se[i]^2, J - N)) # Variances of portfolio WITHOUT and WITH errors
      Sigma_e_AllRows[((i - 1) * J + 1):(i * J)] <- Sigma_e_CS # Update the correct segment
    }

    Sigma_e_AllCountries <- diag(C * J) * Sigma_e_AllRows
    SIGMA_yields <- solve(WpcaFull) %*% Sigma_e_AllCountries %*% t(solve(WpcaFull))
    VarYields <- diag(SIGMA_yields)

    names(VarYields) <- rownames(YieldsTS)
    VarYields <- t(t(VarYields))
  }

  return(VarYields)
}
##############################################################################################################
#' Update the variance-covariance matrix from the "JLL joint Sigma" model. Necessary for optimization
#'
#' @param SSZ Variance-covariance matrix from JLL model
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param N number of country-specific spanned factors
#' @param JLLinputs List of inputs from JLL models
#'
#' @keywords internal

Update_SSZ_JLL <- function(SSZ, Z, N, JLLinputs) {
  Sigma_Ye <- t(chol(SSZ)) # in this case, SSZ is the Variance-covariance matrix of the ORTHOGONALIZED dynamics
  JLLinputs$WishSigmas <- 0 # Avoid recomputing the variance-covariance matrices
  JLLPara <- JLL(Z, N, JLLinputs)
  Sigma_Y <- JLLPara$PI %*% Sigma_Ye
  SSZ <- Sigma_Y %*% t(Sigma_Y) # Variance-covariance matrix non-orthogonalized

  return(SSZ)
}
##############################################################################################################
#' Prepare outputs to export after the model optimization
#'
#' @param Y matrix of yields used in estimation (J x T or CJ x T)
#' @param Z complete set of spanned and unspanned factors (F x T)
#' @param mat vector of maturities (in years) of yields used in estimation (J x 1)
#' @param N number of country-specific spanned factors
#' @param dt time interval unit of the model (scalar)
#' @param Wpca matrix of weights of the portfolios observed without errors (N x J or CN x J)
#' @param K1XQ risk-neutral feedback matrix (N x N or CN x CN)
#' @param SSZ variance-covariance matrix (F x F)
#' @param LoadAs list containing the A loadings
#' @param LoadBs list containing the B loadings
#' @param r0 long-run interest rate (scalar or vector with length C)
#' @param se Variance of the portfolio of yields observed with error (scalar).
#' @param K0Z intercept from the P-dynamics (F x 1)
#' @param K1Z feedback matrix from the P-dynamics (F x F)
#' @param Gy.0 matrix of contemporaneous terms from the P-dynamics (F x F)
#' @param VarYields variance-covariance matrix of the bond yields
#' @param y likelihood of each time series (Tx1)
#' @param GVARinputs List of inputs from GVAR models
#' @param JLLinputs List of inputs from JLL models
#' @param Economies string containing the names of the economy to be estimated
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BS_out Bootstrap output. Default is FALSE.
#'
#' @keywords internal

OptOutputs <- function(Y, Z, mat, N, dt, Wpca, K1XQ, SSZ, LoadAs, LoadBs, r0, se, K0Z, K1Z, Gy.0, VarYields,
                       y, GVARinputs, JLLinputs, Economies, ModelType, BS_out = FALSE) {
  # a) Build LIST 1
  # a.1) List "input"
  Inputs <- list(Y = Y, AllFactors = Z, mat = mat, N = N, dt = dt, Wpca = Wpca)
  if (ModelType %in% c("GVAR single", "GVAR multi")) {
    Inputs$Wgvar <- GVARinputs$Wgvar
  } else if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    Inputs$DomUnit <- JLLinputs$DomUnit
  }

  # a.2) List for the full set of parameter estimates
  Max_llk <- mean(-y)
  ModEst <- list(Max_llk = Max_llk)

  ModEst$Q <- list(K1XQ = K1XQ, r0 = r0, se = se, VarYields = VarYields)
  ModEst$P <- list(SSZ = SSZ, K0Z = K0Z, K1Z = K1Z, Gy.0 = Gy.0)

  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    JLLinputs$WishSigmas <- 1
    JLLinputs$SigmaNonOrtho <- SSZ
    JLLoutcomesOrtho <- JLL(Z, N, JLLinputs)
    # Remove the non-orthogonalized outputs
    JLLoutcomesOrtho$k0 <- NULL
    JLLoutcomesOrtho$k1 <- NULL
    if (BS_out) {
      JLLoutcomesOrtho$a_W <- NULL
      JLLoutcomesOrtho$b <- NULL
      JLLoutcomesOrtho$c <- NULL
      JLLoutcomesOrtho$a_DU_CS <- NULL
      JLLoutcomesOrtho$PIb <- NULL
      JLLoutcomesOrtho$PIac <- NULL
      JLLoutcomesOrtho$Ye <- NULL
    }

    # Add the zeros to the variance-covariance matrix (to avoid rounding problems)
    ZeroIdxVarCov <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$VarCovOrtho
    ZeroIdxSigma_Ye <- JLLoutcomesOrtho$Sigmas$ZeroIdxSigmaJLLOrtho$Sigma_Ye
    JLLoutcomesOrtho$Sigmas$VarCov_Ortho[ZeroIdxVarCov] <- 0
    JLLoutcomesOrtho$Sigmas$Sigma_Ye[ZeroIdxSigma_Ye] <- 0

    PI <- JLLoutcomesOrtho$PI
    JLLoutcomesOrtho$Sigmas$Sigma_Y <- PI %*% JLLoutcomesOrtho$Sigmas$Sigma_Ye
    ModEst$P$JLLoutcomes <- JLLoutcomesOrtho
  }


  # a.3) List of pricing loadings
  ModEst$Q$Load$X <- list(B = LoadBs$BnX, A = LoadAs$AnX, SS = LoadBs$SSX)

  # Q dynamics: Pt as risk factors: % PCN = U0 + U1*X, where U0 = W*AnX and U1 = W*BnX
  U1 <- Wpca %*% LoadBs$BnX
  U0 <- Wpca %*% LoadAs$AnX

  ModEst$Q$Load$P <- Rotate_Lat_Obs(ModEst$Q$Load$X, U1, U0)

  row.names(ModEst$Q$Load$P$A) <- row.names(Y)
  row.names(ModEst$Q$Load$P$B) <- row.names(Y)
  idxSpanned <- LoadBs$idxSpa
  colnames(ModEst$Q$Load$P$B) <- row.names(Z)[idxSpanned]

  if (BS_out) {
    ModEst$Q$Load$X <- NULL # Delete the parameters related to the latent factors
    # (Not needed for the Bootstrap, save some time and pc memory)
  }

  # Summary list to export
  Out <- list(Inputs = Inputs, ModEst = ModEst)
  return(Out)
}
###################################################################################################
#' Rotate latent states to observed ones
#'
#' @param X_List    list of risk-neutral parameters: (i) intercept (A_X); slope (B_X) and volatility matrix (SS_X)
#' @param U1  rotation matrix (N x N)
#' @param U0  rotation vector (N x 1)
#'
#' @keywords internal

Rotate_Lat_Obs <- function(X_List, U1, U0) {
  N <- dim(U1)[1]
  U1_inv <- solve(U1)
  Obs_List <- list()

  # 1) Get volatility matrix
  M <- dim(X_List$SS)[[2]] / N - 1
  SSi <- array(X_List$SS, c(N, N, M + 1))
  SSi <- U1 %*% drop(SSi) %*% t(U1) # Sigma(Z) = U1*Sigma(X)*(U1)'
  Obs_List$SS <- array(SSi, c(N, N * (M + 1)))

  # 2) Get loadings A and B
  Obs_List$B <- if (N == 1) matrix(X_List$B) %*% solve(U1) else X_List$B %*% U1_inv # BX* BZ = U1 --> BZ = BX*U1^(-1)
  Obs_List$A <- X_List$A - Obs_List$B %*% U0 # AX = AZ + BZ*U0 --> AZ = AX - BZ*U0

  return(Obs_List)
}

##########################################################################################################
#' Safe matrix inversion with conditioning check
#'
#' @param A matrix
#' @param B matrix
#' @param reg  Numeric scalar (default: \code{1e-12}). Small ridge term added to the diagonal of \code{A} when the system is ill-conditioned.
#' @param rcond_tol Numeric, default = 1e-12. Tolerance threshold for the reciprocal condition number.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal

safe_solve <- function(A, B = NULL, reg = 1e-12, rcond_tol = 1e-12, verbose = FALSE) {
  A <- as.matrix(A)
  n <- ncol(A)

  # 0) Quick check: size/NA
  if (any(is.na(A))) {
    if (verbose) message("safe_solve: matrix contains NA")
    return(NULL)
  }

  # 1) Try direct solve first
  direct <- tryCatch(
    {
      if (is.null(B)) solve(A) else solve(A, B)
    },
    error = function(e) NULL,
    warning = function(w) NULL
  )

  if (!is.null(direct) && all(is.finite(direct))) {
    return(direct)
  }

  # 2) Compute SVD to estimate conditioning
  sv <- svd(A)
  if (length(sv$d) == 0) {
    return(NULL)
  }
  rcond_est <- sv$d[length(sv$d)] / sv$d[1]

  if (verbose) message(sprintf("safe_solve: rcond_est = %.3e", rcond_est))

  # If well conditioned enough, attempt direct solve again
  if (!is.na(rcond_est) && rcond_est > rcond_tol) {
    out <- tryCatch(
      {
        if (is.null(B)) solve(A) else solve(A, B)
      },
      error = function(e) NULL
    )
    if (!is.null(out)) {
      return(out)
    }
  }

  # 3) Try tiny regularization: (A + reg*I)^{-1}
  regA <- A + reg * diag(n)
  reg_try <- tryCatch(
    {
      if (is.null(B)) solve(regA) else solve(regA, B)
    },
    error = function(e) NULL
  )

  if (!is.null(reg_try) && all(is.finite(reg_try))) {
    if (verbose) message("safe_solve: used regularized solve")
    return(reg_try)
  }

  # 4) Fallback: SVD pseudo-inverse
  tol_svd <- max(dim(A)) * max(sv$d) * .Machine$double.eps
  d_inv <- ifelse(sv$d > tol_svd, 1 / sv$d, 0)
  A_pinv <- sv$v %*% diag(d_inv, length(d_inv)) %*% t(sv$u)

  if (is.null(B)) {
    return(A_pinv)
  } else {
    return(A_pinv %*% B)
  }
}
###########################################################################################
#' Robust method for matrix inversion
#'
#' @param M squared matrix
#'
#' @keywords internal

InvMat_Robust <- function(M) {
  rc <- tryCatch(rcond(M), error = function(e) 0)
  if (any(!is.finite(M)) || rc < 1e-12) {
    # Fallback: add tiny ridge before inversion
    eps <- 1e-10 * mean(diag(M))
    invM <- solve(M + eps * diag(ncol(M)))
  } else {
    invM <- solve(M)
  }

  return(invM)
}
############################################################################################
#' Check for presence of NAs and infinite in numeric variables
#'
#' @param x variable to be tested (numeric)
#' @param name variable labels
#'
#' @keywords internal

check_numeric <- function(x, name) {
  if (any(!is.finite(x))) {
    stop(sprintf("Non-finite values detected in '%s'.", name))
  }
}
