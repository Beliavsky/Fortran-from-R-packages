#' Constructs the model numerical outputs (model fit, IRFs, GIRFs, FEVDs, GFEVDs, and term premia)
#'
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param ModelPara list. Point estimates of the model parameters. See outputs from \code{\link{Optimization}}
#' @param InputsForOutputs list. Inputs for generating IRFs, GIRFs, FEVDs, GFEVDs, and Term Premia.
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param Folder2save Folder path where the outputs will be stored. Default option saves the outputs in a temporary directory.
#' @param verbose Logical flag controlling function messaging. Default is TRUE.
#'
#' @examples
#' data("ParaSetEx")
#' data("InpForOutEx")
#' # Adjust inputs according to the loaded features
#' ModelType <- "JPS original"
#' Economy <- "Brazil"
#' FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)
#'
#' NumOut <- NumOutputs(ModelType, ParaSetEx, InpForOutEx, FacLab, Economy,
#'   Folder2save = NULL, verbose = FALSE
#' )
#'
#' @section Available methods:
#' - `autoplot(object, type)`
#'
#' @returns
#' An object of class 'ATSMNumOutputs' containing the following keys elements:
#' \enumerate{
#' \item Model parameter estimates
#' \item Model fit of bond yields
#' \item IRFs
#' \item FEVDs
#' \item GIRFs
#' \item GFEVDs
#' \item Bond yield decomposition
#' }
#'
#' @details
#' Both IRFs and FEVDs are computed using the Cholesky decomposition method. The risk factors are ordered as follows: (i) global unspanned factors, and (ii) domestic unspanned and spanned factors for each country. The order of countries follows the sequence defined in the \code{Economies} vector.
#'
#' @references
#' Pesaran, H. Hashem, and Shin, Yongcheol. "Generalized impulse response analysis in linear multivariate models." Economics letters 58.1 (1998): 17-29.
#' @export

NumOutputs <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, Folder2save = NULL,
                       verbose = TRUE) {
  if (verbose) message("2.2) Computing numerical outputs")

  AllNumOutputs <- list()
  AllNumOutputs <- OutputConstruction(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, verbose)

  # Save relevant numerical outputs
  PEoutputs <- list(ModelPara, AllNumOutputs)
  names(PEoutputs) <- c("Model Parameters", "Numerical Outputs")

  FolderPath <- if (is.null(Folder2save)) tempdir() else Folder2save
  saveRDS(PEoutputs, paste(FolderPath, "/PEoutputs_", InputsForOutputs$"Label Outputs", ".rds", sep = ""))

  # Generate graphs, if previously selected
  GraphicalOutputs(
    ModelType, ModelPara, AllNumOutputs, InputsForOutputs, Economies, FactorLabels,
    FolderPath, verbose
  )

  # Create model output object
  attr(AllNumOutputs, "NumOuts") <- list(
    ModelType = ModelType,
    ModelPara = ModelPara,
    NumOut = AllNumOutputs,
    Inputs = InputsForOutputs,
    Economies = Economies,
    FactorLabels = FactorLabels
  )

  # Return the structured Outputs object
  return(structure(AllNumOutputs, class = "ATSMNumOutputs"))
}

######################################################################################################
######################################################################################################
#' Numerical outputs (variance explained, model fit, IRFs, GIRFs, FEVDs, GFEVDs, and risk premia decomposition)
#' for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelPara list of model parameter estimates (See the "Optimization" function)
#' @param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param verbose Logical flag controlling function messaging. Default is TRUE.
#'
#' @keywords internal

OutputConstruction <- function(ModelType, ModelPara, InputsForOutputs, FactorLabels, Economies, verbose = TRUE) {
  # Output summary
  # Total Variance Explained and Model fit
  if (verbose) message(" ** Total Variance explained")
  Total_Var_exp <- VarianceExplained(ModelType, ModelPara, FactorLabels, Economies)

  if (verbose) message(" ** Model Fit")
  ModFit <- YieldsFit(ModelType, ModelPara, FactorLabels, Economies)

  # IRF and GIRF
  if (verbose) message(" ** IRFs and GIRFs")
  IRFout <- IRFandGIRF(ModelType, ModelPara, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)

  # FEVD and GFEVD
  if (verbose) message(" ** FEVDs and GFEVDs")
  FEVDout <- FEVDandGFEVD(ModelType, ModelPara, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # Risk Premia Decomposition
  if (verbose) message(" ** Term Premia")
  TermPremia <- TermPremiaDecomp(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies)


  NumericalOutputs <- list(Total_Var_exp, ModFit, IRFout$IRFs, FEVDout$FEVDs, IRFout$GIRFs, FEVDout$GFEVDs, TermPremia)
  names(NumericalOutputs) <- c("PC var explained", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremiaDecomp")

  return(NumericalOutputs)
}

######################################################################################################
####################### 1) Total Variance Explained #########################################
######################################################################################################
#' Percentage explained by the spanned factors of the variations in the set of observed yields for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelPara List of model parameter estimates (see the "Optimization" function)
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

VarianceExplained <- function(ModelType, ModelPara, FactorLabels, Economies) {
  N <- length(FactorLabels$Spanned)
  Total_Var_exp <- list()

  # Models estimated individually
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    for (i in 1:length(Economies)) {
      H <- eigen(stats::cov(t(ModelPara[[ModelType]][[Economies[i]]]$Inputs$Y)))$values
      percentages_explained <- cumsum(H) / sum(H)
      Total_Var_exp[[i]] <- percentages_explained[1:N]
    }
  } else {
    # Models estimated jointly
    J <- length(ModelPara[[ModelType]]$Inputs$mat)
    idx0 <- 0
    for (i in 1:length(Economies)) {
      idx1 <- idx0 + J
      H <- eigen(stats::cov(t(ModelPara[[ModelType]]$Inputs$Y[(idx0 + 1):idx1, ])))$values
      percentages_explained <- cumsum(H) / sum(H)
      Total_Var_exp[[i]] <- percentages_explained[1:N]
      idx0 <- idx1
    }
  }
  names(Total_Var_exp) <- Economies

  return(Total_Var_exp)
}

######################################################################################################
########################################### 2) Model Fit #############################################
######################################################################################################
#' Computes two measures of model fit for bond yields (all models)
#'
#' @param ModelType a string-vector containing the label of the model to be estimated
#' @param ModelPara List of model parameter estimates (See the "Optimization" function)
#' @param FactorLabels a string-list based which contains the labels of all the variables present in the model
#' @param Economies a string-vector containing the names of the economies which are part of the economic system
#'
#' @details
#' "Model-implied yields" is the measure of fit based exclusively on the risk-neutral parameters, whereas the
#' "Model-Fit" takes into account both the risk-neutral and the physical parameters.
#'
#' @references
#' See, for instance, Jotiskhatira, Le and Lundblad (2015). "Why do interest rates in different currencies co-move?" (Journal of Financial Economics)
#'
#' @keywords internal

YieldsFit <- function(ModelType, ModelPara, FactorLabels, Economies) {
  # Extract dimensions
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)
  Output <- list()

  # Helper function to compute spanned factors
  get_spanned_factors <- function(Z, LabelSpannedCS, AllLabels) {
    IdxSpanned <- match(LabelSpannedCS, AllLabels)
    return(Z[IdxSpanned, ])
  }

  # Helper function to compute model fit
  compute_model_fit <- function(Afull, Bspanned, P, J, T, YieldData) {
    return(Y_Fit(Afull, Bspanned, P, J, T, dimnames(YieldData)))
  }

  # Helper function to compute model-implied yields
  compute_model_implied_yields <- function(Afull, Bfull, K0Z, K1Z, Z, J, T, YieldData) {
    return(Y_ModImp(Afull, Bfull, K0Z, K1Z, Z, J, T, dimnames(YieldData)))
  }

  # I) Models estimated individually
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    mat <- ModelPara[[ModelType]][[Economies[1]]]$Inputs$mat
    J <- length(mat)

    for (i in 1:C) {
      # Extract necessary data
      InfoSet <- ModelPara[[ModelType]][[Economies[i]]]
      YieldData <- InfoSet$Inputs$Y
      Z <- InfoSet$Inputs$AllFactors
      T_dim <- ncol(Z)
      Afull <- InfoSet$ModEst$Q$Load$P$A
      Bspanned <- InfoSet$ModEst$Q$Load$P$B
      K0Z <- InfoSet$ModEst$P$K0Z
      K1Z <- InfoSet$ModEst$P$K1Z

      # Determine AllLabels and LabelSpannedCS
      if (ModelType == "JPS original") {
        AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])
      } else {
        AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
      }
      LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])

      # Get spanned factors
      P <- get_spanned_factors(Z, LabelSpannedCS, AllLabels)

      # Compute model fit and implied yields
      Yieldfit <- compute_model_fit(Afull, Bspanned, P, J, T_dim, YieldData)
      Bfull <- BUnspannedAdapSep(G, M, ModelPara[[ModelType]], Economies, Economies[i], ModelType)
      dimnames(Bfull) <- list(rownames(YieldData), rownames(Z))
      YieldModelImplied <- compute_model_implied_yields(Afull, Bfull, K0Z, K1Z, Z, J, T_dim, YieldData)

      # Store results
      Output[[ModelType]][[Economies[i]]] <- list(
        "Yield Fit" = Yieldfit,
        "Yield Model Implied" = YieldModelImplied
      )
    }
  } else {
    # II) Models estimated jointly
    InfoSet <- ModelPara[[ModelType]]
    Z <- InfoSet$Inputs$AllFactors
    YieldData <- InfoSet$Inputs$Y
    mat <- InfoSet$Inputs$mat
    Afull <- InfoSet$ModEst$Q$Load$P$A
    Bspanned <- InfoSet$ModEst$Q$Load$P$B
    K0Z <- InfoSet$ModEst$P$K0Z
    K1Z <- InfoSet$ModEst$P$K1Z
    J <- length(mat)
    T_dim <- ncol(Z)

    # Compute IdxSpanned
    IdxSpanned <- c()
    idxSpa0 <- G + M
    for (j in 1:C) {
      idxSpa1 <- idxSpa0 + N
      IdxSpanned <- c(IdxSpanned, (idxSpa0 + 1):idxSpa1)
      idxSpa0 <- idxSpa1 + M
    }

    # Get spanned factors
    P <- Z[IdxSpanned, ]

    # Compute model fit and implied yields
    Yieldfit <- compute_model_fit(Afull, Bspanned, P, C * J, T_dim, YieldData)
    Bfull <- BUnspannedAdapJoint(G, M, N, C, J, Bspanned)
    dimnames(Bfull) <- list(rownames(YieldData), rownames(Z))
    YieldModelImplied <- compute_model_implied_yields(Afull, Bfull, K0Z, K1Z, Z, C * J, T_dim, YieldData)

    Output <- list(
      "Yield Fit" = Yieldfit,
      "Yield Model Implied" = YieldModelImplied
    )
  }

  return(Output)
}

######################################################################################################
########################################### 3) IRFs and GIRFs ########################################
######################################################################################################
#' IRFs and GIRFs for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelPara list of model parameter estimates (See the "Optimization" function)
#' @param IRFhoriz single numerical vector containing the desired horizon of analysis for the IRFs
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @details
#' The Structural shocks from the IRFs are identified via Cholesky decomposition
#'
#' @keywords internal

IRFandGIRF <- function(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies) {
  # Pre-allocation
  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  IRFoutputs <- list()
  GIRFoutputs <- list()

  # Helper function to compute IRFs and GIRFs for a single country
  compute_single_country <- function(ModelType, Para_Set, Econ, Economies) {
    J <- length(Para_Set[[Econ]]$Inputs$mat)
    K <- nrow(Para_Set[[Econ]]$Inputs$AllFactors)
    YieldsLabel <- rownames(Para_Set[[Econ]]$Inputs$Y) # Yield labels

    # Summarize inputs for the IRFs
    SIGMA <- Para_Set[[Econ]]$ModEst$P$SSZ # KxK (variance-covariance matrix)
    K1Z <- Para_Set[[Econ]]$ModEst$P$K1Z # KxK (feedback matrix)
    B <- BUnspannedAdapSep(G, M, Para_Set, Economies, Econ, ModelType)

    # Compute IRFs
    IRFs <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Econ)
    IRFoutputs[[ModelType]][[Econ]] <- IRFs # Store Country specific IRFs

    # Compute GIRFs
    G0.y <- Para_Set[[Econ]]$ModEst$P$Gy.0
    GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, Econ)
    GIRFoutputs[[ModelType]][[Econ]] <- GIRFs # Store Country specific GIRFs

    return(list(IRFs = IRFs, GIRFs = GIRFs))
  }

  # Helper function to compute IRFs and GIRFs for joint country models
  compute_joint_country <- function(ModelType, Para_Set, FactorLabels, Economies) {
    JLL_Label <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")
    J <- length(Para_Set$Inputs$mat)
    K <- nrow(Para_Set$Inputs$AllFactors)
    YieldsLabel <- rownames(Para_Set$Inputs$Y) # Yield labels

    # Summarize inputs for the IRFs
    if (ModelType %in% JLL_Label) {
      SIGMA <- Para_Set$ModEst$P$JLLoutcomes$Sigmas$Sigma_Y # For JLL models, we selected the cholesky factor, which won't be computed inside "ComputeIRFs"
    } else {
      SIGMA <- Para_Set$ModEst$P$SSZ # KxK (variance-covariance matrix)
    }
    K1Z <- Para_Set$ModEst$P$K1Z # KxK (feedback matrix)
    BSpanned <- Para_Set$ModEst$Q$Load$P$B
    B <- BUnspannedAdapJoint(G, M, N, C, J, BSpanned)

    # Compute IRFs
    IRFoutputs[[ModelType]] <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType)

    # Compute GIRFs
    G0.y <- Para_Set$ModEst$P$Gy.0
    if (ModelType %in% JLL_Label) {
      SIGMA <- ModelPara[[ModelType]]$ModEst$P$JLLoutcomes$Sigmas$VarCov_NonOrtho
    }
    GIRFs <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType)
    GIRFoutputs[[ModelType]] <- GIRFs # Store Country specific GIRFs

    # Compute orthogonalized IRFs and GIRFs for JLL-based models
    if (ModelType %in% JLL_Label) {
      K1Ze <- Para_Set$ModEst$P$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- Para_Set$ModEst$P$JLLoutcomes$PI
      Se <- Para_Set$ModEst$P$JLLoutcomes$Sigmas$Sigma_Ye
      SIGMA_e <- Para_Set$ModEst$P$JLLoutcomes$Sigmas$VarCov_Ortho

      # Compute IRFs orthogonalized
      IRFOrtho <- list()
      IRFOrtho[[ModelType]] <- ComputeIRFs(Se, K1Ze, B, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType,
        PI = PI, Mode = "Ortho"
      )

      # Gather Outputs
      IRFoutputs[[ModelType]]$Factors <- list(NonOrtho = IRFoutputs[[ModelType]]$Factors, Ortho = IRFOrtho[[ModelType]]$Factors)
      IRFoutputs[[ModelType]]$Yields <- list(NonOrtho = IRFoutputs[[ModelType]]$Yields, Ortho = IRFOrtho[[ModelType]]$Yields)

      # Compute GIRFs orthogonalized
      GIRFsOrtho <- list()
      GIRFsOrtho[[ModelType]] <- ComputeGIRFs(SIGMA_e, K1Ze, B, G0.y, FactorLabels, K, C * J, IRFhoriz, YieldsLabel,
        ModelType,
        PI = PI, Mode = "Ortho"
      )

      # Gather Outputs
      GIRFoutputs[[ModelType]]$Factors <- list(NonOrtho = GIRFoutputs[[ModelType]]$Factors, Ortho = GIRFsOrtho[[ModelType]]$Factors)
      GIRFoutputs[[ModelType]]$Yields <- list(NonOrtho = GIRFoutputs[[ModelType]]$Yields, Ortho = GIRFsOrtho[[ModelType]]$Yields)
    }

    return(list(IRFoutputs = IRFoutputs, GIRFoutputs = GIRFoutputs))
  }


  # 1) SINGLE COUNTRY MODELS
  Para_Set <- ModelPara[[ModelType]]
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (i in 1:C) {
      IRFset_CS <- compute_single_country(ModelType, Para_Set, Economies[i], Economies)
      IRFoutputs[[ModelType]][[Economies[i]]] <- IRFset_CS$IRFs
      GIRFoutputs[[ModelType]][[Economies[i]]] <- IRFset_CS$GIRFs
    }
  } else {
    # 2) JOINT COUNTRY MODELS
    IRFset <- compute_joint_country(ModelType, Para_Set, FactorLabels, Economies)
    IRFoutputs <- IRFset$IRFoutputs
    GIRFoutputs <- IRFset$GIRFoutputs
  }

  Out <- list(IRFs = IRFoutputs, GIRFs = GIRFoutputs)

  return(Out)
}

#########################################################################################################
################################### 4) FEVD and GFEVD ###################################################
#########################################################################################################
#' FEVDs and GFEVDs for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelPara list of model parameter estimates (see the "Optimization" function)
#' @param FEVDhoriz single numerical vector containing the desired horizon of analysis for the FEVDs and GFEVDs
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @details
#' Structural shocks are identified via Cholesky decomposition
#'
#' @keywords internal

FEVDandGFEVD <- function(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies) {
  # Pre-allocation
  C <- length(Economies)
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  FEVDoutputs <- list()
  GFEVDoutputs <- list()

  # Helper function to compute FEVDs and GFEVDs for a single country
  compute_single_country <- function(ModelType, Para_Set, Econ, Economies) {
    K <- nrow(Para_Set[[Econ]]$Inputs$AllFactors)
    J <- length(Para_Set[[Econ]]$Inputs$mat)
    YieldsLabel <- rownames(Para_Set[[Econ]]$Inputs$Y) # Yield labels

    # Summarize inputs for the FEVDs
    SIGMA <- Para_Set[[Econ]]$ModEst$P$SSZ
    K1Z <- Para_Set[[Econ]]$ModEst$P$K1Z
    G0 <- Para_Set[[Econ]]$ModEst$P$Gy.0
    B <- BUnspannedAdapSep(G, M, Para_Set, Economies, Econ, ModelType)

    # Compute FEVDs
    FEVDs <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Econ)

    # Compute GFEVDs
    GFEVDs <- ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, Econ)

    # Return results
    return(list(FEVDs = FEVDs, GFEVDs = GFEVDs))
  }

  # Helper function to compute FEVDs and GFEVDs for joint country models
  compute_joint_country <- function(ModelType, ParaSet, FactorLabels, Economies) {
    J <- length(ParaSet$Inputs$mat)
    K <- nrow(ParaSet$Inputs$AllFactors)
    YieldsLabel <- rownames(ParaSet$Inputs$Y) # Yield labels

    # Summarize inputs for the FEVD
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
      Chol_Fac_JLL <- ParaSet$ModEst$P$JLLoutcomes$Sigmas$Sigma_Y
    }

    K1Z <- ParaSet$ModEst$P$K1Z # KxK (feedback matrix)
    SIGMA <- ParaSet$ModEst$P$SSZ # KxK (variance-covariance matrix)
    BSpanned <- ParaSet$ModEst$Q$Load$P$B
    B <- BUnspannedAdapJoint(G, M, N, C, J, BSpanned)
    G0 <- ParaSet$ModEst$P$Gy.0

    # Compute FEVD
    FEVDoutputs[[ModelType]] <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C * J, FEVDhoriz, YieldsLabel,
      ModelType,
      CholFac_JLL = Chol_Fac_JLL
    )

    # Compute GFEVD
    GFEVDoutputs[[ModelType]] <- ComputeGFEVDs(
      SIGMA, K1Z, G0, B, FactorLabels, K, C * J, FEVDhoriz, YieldsLabel,
      ModelType
    )

    # Compute orthogonalized FEVDs and GFEVDs for JLL-based models
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
      K1Ze <- ParaSet$ModEst$P$JLLoutcomes$k1_e # KxK (feedback matrix)
      PI <- ParaSet$ModEst$P$JLLoutcomes$PI
      SIGMA_Ortho <- ParaSet$ModEst$P$JLLoutcomes$Sigmas$VarCov_Ortho
      Chol_Fac_JLL_Ortho <- ParaSet$ModEst$P$JLLoutcomes$Sigmas$Sigma_Ye

      # Compute FEVDs orthogonalized
      FEVDOrtho <- list()
      FEVDOrtho[[ModelType]] <- ComputeFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C * J, FEVDhoriz, YieldsLabel,
        ModelType,
        CholFac_JLL = Chol_Fac_JLL_Ortho, PI = PI, Mode = "Ortho"
      )

      # Gather Outputs
      FEVDoutputs[[ModelType]]$Factors <- list(NonOrtho = FEVDoutputs[[ModelType]]$Factors, Ortho = FEVDOrtho[[ModelType]]$Factors)
      FEVDoutputs[[ModelType]]$Yields <- list(NonOrtho = FEVDoutputs[[ModelType]]$Yields, Ortho = FEVDOrtho[[ModelType]]$Yields)

      # Compute GFEVDs orthogonalized
      GFEVDsOrtho <- list()
      GFEVDsOrtho[[ModelType]] <- ComputeGFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C * J, FEVDhoriz, YieldsLabel,
        ModelType,
        PI = PI, Mode = "Ortho"
      )

      # Gather Outputs
      GFEVDoutputs[[ModelType]]$Factors <- list(NonOrtho = GFEVDoutputs[[ModelType]]$Factors, Ortho = GFEVDsOrtho[[ModelType]]$Factors)
      GFEVDoutputs[[ModelType]]$Yields <- list(NonOrtho = GFEVDoutputs[[ModelType]]$Yields, Ortho = GFEVDsOrtho[[ModelType]]$Yields)
    }

    # Return results
    return(list(FEVDoutputs = FEVDoutputs, GFEVDoutputs = GFEVDoutputs))
  }

  # 1) SINGLE COUNTRY MODELS
  Para_Set <- ModelPara[[ModelType]]
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    for (i in 1:C) {
      FEVDset_CS <- compute_single_country(ModelType, Para_Set, Economies[i], Economies)
      FEVDoutputs[[ModelType]][[Economies[i]]] <- FEVDset_CS$FEVDs
      GFEVDoutputs[[ModelType]][[Economies[i]]] <- FEVDset_CS$GFEVDs
    }
  } else {
    # 2) JOINT COUNTRY MODELS
    FEVDset <- compute_joint_country(ModelType, Para_Set, FactorLabels, Economies)
    FEVDoutputs <- FEVDset$FEVDoutputs
    GFEVDoutputs <- FEVDset$GFEVDoutputs
  }

  Out <- list(FEVDs = FEVDoutputs, GFEVDs = GFEVDoutputs)
  return(Out)
}

######################################################################################################
####################### 5) Risk Premia Decomposition #################################################
######################################################################################################
#' Decomposition of yields into the average of expected future short-term interest rate and risk premia for all models
#'
#' @param ModelPara list of model parameter estimates (see the "Optimization" function)
#' @param FactorLabels  string-list based which contains all the labels of all the variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#' @param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#'
#'
#' @keywords internal

TermPremiaDecomp <- function(ModelPara, FactorLabels, ModelType, InputsForOutputs, Economies) {
  WishFP <- InputsForOutputs$ForwardPremia

  # 1) Compute the country-specific expected components
  EP_list <- ExpectedComponent(ModelPara, InputsForOutputs, ModelType, Economies, FactorLabels, WishFP)

  # 2) Compute Term Premium
  OutputTP <- TermPremia(ModelPara, EP_list$avexp, ModelType, Economies)

  # 3) Forward Premia (if desired)
  if (WishFP) {
    OutputFP <- ForwardPremia(ModelPara, EP_list$avexpFP, ModelType, FactorLabels, InputsForOutputs, Economies)
  } else {
    OutputFP <- NA
  }

  Output <- list(RiskPremia = OutputTP, ForwardPremia = OutputFP)

  return(Output)
}

######################################################################################################
######################################## AUXILIARY FUNCTIONS #########################################
######################################################################################################
#' Model-implied yields (cross-section)
#'
#' @param ALoad  A loadings
#' @param BLoad B loadings
#' @param Spa_TS time series of spanned factors
#' @param MatLength length of the vector of maturities
#' @param TDim Time-series dimension
#' @param YieldLab Label of yields
#'
#' @keywords internal

Y_Fit <- function(ALoad, BLoad, Spa_TS, MatLength, TDim, YieldLab) {
  Yieldfit <- matrix(NA, nrow = MatLength, ncol = TDim)
  if (is.null(dim(Spa_TS))) {
    for (h in 1:TDim) {
      Yieldfit[, h] <- ALoad + BLoad %*% Spa_TS[h]
    }
  } else {
    for (h in 1:TDim) {
      Yieldfit[, h] <- ALoad + BLoad %*% Spa_TS[, h]
    }
  }
  dimnames(Yieldfit) <- YieldLab

  return(Yieldfit)
}

##############################################################
#' Model-implied yields (P-dynamics)
#'
#' @param ALoad  A loadings
#' @param BLoad B loadings
#' @param K0Z intercept from the P-dynamics
#' @param K1Z feedback matrix from the P-dynamics
#' @param PdynFact time series of the risk-factors spanned factors
#' @param MatLength length of the vector of maturities
#' @param TDim Time-series dimension
#' @param YieldLab Label of yields
#'
#'
#' @keywords internal

Y_ModImp <- function(ALoad, BLoad, K0Z, K1Z, PdynFact, MatLength, TDim, YieldLab) {
  YieldModelImplied <- matrix(NA, nrow = MatLength, ncol = TDim)
  for (h in 2:TDim) { #  first observation is discarded
    YieldModelImplied[, h] <- ALoad + BLoad %*% (K0Z + K1Z %*% PdynFact[, h - 1])
  }
  dimnames(YieldModelImplied) <- YieldLab

  return(YieldModelImplied)
}

#####################################################################################################
#' Transform B_spanned into B_unspanned for jointQ models
#'
#' @param G number of global unspanned factors
#' @param M number of domestic unspanned factors
#' @param N number of domestic spanned factors
#' @param C number of economies of the economic system
#' @param J number of country-specific observed bond yields
#' @param BSpanned B that accomodates only the map to the spanned factors only
#'
#' @keywords internal

BUnspannedAdapJoint <- function(G, M, N, C, J, BSpanned) {
  K <- C * (N + M) + G
  CJ <- C * J

  BUnspanned <- matrix(0, nrow = CJ, ncol = K)

  idxA <- 0
  idxB <- G + M
  idxC <- 0

  for (i in 1:C) {
    idxAA <- idxA + J
    idxBB <- idxB + N
    idxCC <- idxC + N
    BUnspanned[(idxA + 1):idxAA, (idxB + 1):idxBB] <- BSpanned[(idxA + 1):idxAA, (idxC + 1):idxCC]
    idxA <- idxAA
    idxB <- idxBB + M
    idxC <- idxCC
  }

  return(BUnspanned)
}

#################################################################################################
#' Transform B_spanned into B_unspanned for sepQ models
#'
#' @param G number of global unspanned factors
#' @param M number of domestic unspanned factors per country
#' @param ModelPara_Short list of model parameter estimates (See the "Optimization" function)
#' @param Economies complet set of economies of the economic system
#' @param Economy  specific economy under study
#' @param ModelType a string-vector containing the label of the model to be estimated
#'
#' @keywords internal

BUnspannedAdapSep <- function(G, M, ModelPara_Short, Economies, Economy, ModelType) {
  i <- match(Economy, Economies)
  C <- length(Economies)
  N <- ModelPara_Short[[Economy]]$Inputs$N
  J <- length(ModelPara_Short[[Economy]]$Inputs$mat)


  if (ModelType == "JPS original") {
    K <- N + M + G
    BUnspanned <- matrix(0, nrow = J, ncol = K)
    BSpanned <- ModelPara_Short[[Economies[i]]]$ModEst$Q$Load$P$B
    BUnspanned[, (K - N + 1):K] <- BSpanned
  } else if (ModelType %in% c("JPS global", "GVAR single")) {
    K <- C * (N + M) + G
    BUnspanned <- matrix(0, nrow = J, ncol = K)
    BSpanned <- ModelPara_Short[[Economies[i]]]$ModEst$Q$Load$P$B

    IDX <- list()
    idx0 <- G + M
    for (h in 1:C) {
      idx1 <- idx0 + N
      IDX[[h]] <- (idx0 + 1):idx1
      idx0 <- idx1 + M
    }

    BUnspanned[, IDX[[i]]] <- BSpanned
  }

  return(BUnspanned)
}
#####################################################################################################
#' Compute IRFs of all models
#'
#' @param SIGMA Variance-covariance matrix
#' @param K1Z Loading As
#' @param BLoad Loading Bs
#' @param FactorLabels List containing the label of factors
#' @param FacDim Dimension of the P-dynamics
#' @param MatLength Length of the maturity vector
#' @param IRFhoriz Horizon of the analysis
#' @param YieldsLabel Label of bond yields
#' @param ModelType Desired model type
#' @param Economy specific economy under study
#' @param PI matrix PI for JLL-based models
#' @param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#' @keywords internal

ComputeIRFs <- function(SIGMA, K1Z, BLoad, FactorLabels, FacDim, MatLength, IRFhoriz, YieldsLabel, ModelType,
                        Economy = NULL, PI = NULL, Mode = FALSE) {
  # 1) Initialization of IRFs of interest
  tempFactors <- array(0, c(FacDim, FacDim, IRFhoriz))
  tempYields <- array(0, c(MatLength, FacDim, IRFhoriz))


  # 2) Compute the IRFs
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    AdjTerm <- PI
  } else {
    AdjTerm <- diag(FacDim)
  }

  # Choleski term
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    S <- SIGMA
  } else {
    S <- t(chol(SIGMA))
  }

  # Shock at t=0:
  tempFactors[, , 1] <- S
  tempYields[, , 1] <- BLoad %*% AdjTerm %*% S
  # Shock at t=1:
  for (r in 2:IRFhoriz) {
    if (r == 2) {
      A1h <- K1Z
    } else {
      A1h <- A1h %*% K1Z
    }
    tempFactors[, , r] <- A1h %*% S # IRF (t+h) = A1^h*S
    tempYields[, , r] <- BLoad %*% AdjTerm %*% A1h %*% S
  }
  IRFRiskFactors <- aperm(tempFactors, c(3, 1, 2))
  IRFYields <- aperm(tempYields, c(3, 1, 2))

  Horiz <- t(t(0:(IRFhoriz - 1))) # Add a column for horizon of interest

  # 3) Adjust the variable labels
  # Factor Labels
  if (ModelType == "JPS original") {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  } else {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  dimnames(IRFRiskFactors) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(IRFYields) <- list(Horiz, YieldsLabel, AllFactorsLabels)

  Out <- list(Factors = IRFRiskFactors, Yields = IRFYields)

  return(Out)
}

##########################################################################################################
#' Compute GIRFs for all models
#'
#' @param Sigma.y Variance-covariance matrix
#' @param F1 Feedback matrix
#' @param BLoad Loading Bs
#' @param G0.y matrix of contemporaneous terms
#' @param FactorLabels List containing the labels of the factors
#' @param FacDim Dimension of the P-dynamics
#' @param MatLength Length of the maturity vector
#' @param GIRFhoriz Horizon of the analysis
#' @param YieldsLabel Label o yields
#' @param ModelType desired Model type
#' @param Economy Economy under study
#' @param PI matrix PI for JLL-based models
#' @param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#' #' @references
#' \itemize{
#' \item This function is partially based on the version of the "irf" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#' @keywords internal

ComputeGIRFs <- function(Sigma.y, F1, BLoad, G0.y, FactorLabels, FacDim, MatLength, GIRFhoriz, YieldsLabel,
                         ModelType, Economy = NULL, PI = NULL, Mode = FALSE) {
  # 1) Dynamic multiplier:
  Ry.h <- array(NA, c(FacDim, FacDim, GIRFhoriz))
  Ry.h[, , 1] <- diag(FacDim) # dynamic multiplier at t=0

  for (w in 2:GIRFhoriz) {
    Ry.h[, , w] <- F1 %*% Ry.h[, , w - 1]
  }

  # 2) Build the vector containing the one unit-shock for each variable of the system
  ey.j <- diag(FacDim)

  # 3) GIRFs:
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    AdjTerm <- PI
  } else {
    AdjTerm <- diag(FacDim)
  }

  # 3.1) Factors
  AllResponsesToAllShocksFactors <- array(NA, c(FacDim, GIRFhoriz, FacDim))
  AllResponsesToShockOfOneVariableFactors <- matrix(NA, ncol = GIRFhoriz, nrow = FacDim)
  for (g in 1:FacDim) {
    for (w in 1:GIRFhoriz) {
      numFactors <- AdjTerm %*% (Ry.h[, , w] %*% solve(G0.y) %*% Sigma.y %*% ey.j[, g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
      demFactors <- 1 / sqrt((t(ey.j[, g]) %*% Sigma.y %*% ey.j[, g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
      AllResponsesToShockOfOneVariableFactors[, w] <- numFactors * drop(demFactors)
    }
    AllResponsesToAllShocksFactors[, , g] <- AllResponsesToShockOfOneVariableFactors
  }

  GIRFFactors <- aperm(AllResponsesToAllShocksFactors, c(2, 1, 3))

  # 3.2) Yields
  AllResponsesToAllShocksYields <- array(NA, c(MatLength, GIRFhoriz, FacDim))
  AllResponsesToShockOfOneVariableYields <- matrix(NA, ncol = GIRFhoriz, nrow = MatLength)
  for (g in 1:FacDim) {
    for (w in 1:GIRFhoriz) {
      numYields <- BLoad %*% AdjTerm %*% (Ry.h[, , w] %*% solve(G0.y) %*% Sigma.y %*% ey.j[, g]) # numerator from equation at the bottom of the page 22 (PS, 1998)
      demYields <- 1 / sqrt((t(ey.j[, g]) %*% Sigma.y %*% ey.j[, g])) # denominator from equation at the bottom of the page 22 (PS, 1998)
      AllResponsesToShockOfOneVariableYields[, w] <- numYields * drop(demYields)
    }
    AllResponsesToAllShocksYields[, , g] <- AllResponsesToShockOfOneVariableYields
  }

  GIRFYields <- aperm(AllResponsesToAllShocksYields, c(2, 1, 3))

  # 4) Prepare labels for the output
  if (ModelType == "JPS original") {
    labelsGIRF <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", "GVAR single"))) {
    labelsGIRF <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  } else {
    labelsGIRF <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  # 4.1) Add columns containing the horizons
  Horiz <- t(t(0:(GIRFhoriz - 1))) # Add a column for horizon of interest

  # 4.2) Labels
  dimnames(GIRFFactors) <- list(Horiz, labelsGIRF, labelsGIRF)
  dimnames(GIRFYields) <- list(Horiz, YieldsLabel, labelsGIRF)

  GIRFoutputs <- list(Factors = GIRFFactors, Yields = GIRFYields)

  return(GIRFoutputs)
}

#######################################################################################################
#' Compute FEVDs for all models
#'
#' @param SIGMA Variance-covariance matrix
#' @param K1Z Loading As
#' @param G0 contemporaneous terms
#' @param BLoad Loading Bs
#' @param FactorLabels List containing the label of factors
#' @param FacDim Dimension of the P-dynamics
#' @param MatLength Length of the maturity vector
#' @param FEVDhoriz Horizon of the analysis
#' @param YieldsLabel Label of bond yields
#' @param ModelType Desired model type
#' @param Economy specific economy under study
#' @param CholFac_JLL Cholesky factorization term from JLL models
#' @param PI matrix PI for JLL-based models
#' @param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#' @references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#' @keywords internal

ComputeFEVDs <- function(SIGMA, K1Z, G0, BLoad, FactorLabels, FacDim, MatLength, FEVDhoriz, YieldsLabel,
                         ModelType, Economy = NULL, CholFac_JLL = NULL, PI = NULL, Mode = FALSE) {
  # 1) Dynamic multipliers
  Ry.h <- array(NA, c(FacDim, FacDim, FEVDhoriz))
  Ry.h[, , 1] <- diag(FacDim) # dynamic multiplier at t=0

  for (l in 2:FEVDhoriz) {
    Ry.h[, , l] <- K1Z %*% Ry.h[, , l - 1]
  }

  # 2) Initialization
  vslct <- diag(FacDim)
  eslct <- diag(FacDim)

  # 2.1) Minor preliminary work
  invG <- diag(nrow(G0)) / G0
  invG[!is.finite(invG)] <- 0
  invGSigmau <- solve(G0) %*% SIGMA

  # Choleski term
  if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    P <- CholFac_JLL
  } else {
    P <- t(chol(invGSigmau))
  }

  # Adjustment term
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    AdjTerm <- PI
  } else {
    AdjTerm <- diag(FacDim)
  }

  scale <- 1

  # 3) FEVD
  # 3.1) Factors
  FEVDresFactors <- array(NA, c(nrow = FacDim, ncol = FacDim, FEVDhoriz))
  num <- matrix(0, nrow = FacDim, ncol = FacDim)
  den <- rep(0, times = FacDim)

  for (l in 1:FEVDhoriz) {
    acc1 <- (eslct %*% Ry.h[, , l] %*% P %*% vslct)^2
    num <- num + acc1
    acc2 <- diag(eslct %*% Ry.h[, , l] %*% invGSigmau %*% t(invG) %*% t(Ry.h[, , l]) %*% eslct)
    den <- den + acc2
    FEVDresFactors[, , l] <- scale * num / den
  }

  FEVDFactors <- aperm(FEVDresFactors, c(3, 2, 1))

  # 3.2) Yields
  eslctCJ <- diag(MatLength)
  vslctCJ <- diag(FacDim)

  FEVDresYields <- array(NA, c(nrow = MatLength, ncol = FacDim, FEVDhoriz))
  num <- matrix(0, nrow = MatLength, ncol = FacDim)
  den <- matrix(0, nrow = MatLength, ncol = FacDim)

  for (l in 1:FEVDhoriz) {
    acc1 <- (eslctCJ %*% BLoad %*% AdjTerm %*% Ry.h[, , l] %*% P %*% vslctCJ)^2
    num <- num + acc1
    acc2 <- diag(eslctCJ %*% BLoad %*% AdjTerm %*% Ry.h[, , l] %*% invGSigmau %*% t(invG) %*% t(Ry.h[, , l]) %*% t(AdjTerm) %*% t(BLoad) %*% eslctCJ)
    den <- den + acc2
    FEVDresYields[, , l] <- scale * num / den
  }

  FEVDYields <- aperm(FEVDresYields, c(3, 2, 1))

  # 4) Prepare labels
  if (ModelType == "JPS original") {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  } else {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  # 5) Export outputs
  Horiz <- 1:(FEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error

  dimnames(FEVDFactors) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(FEVDYields) <- list(Horiz, AllFactorsLabels, YieldsLabel)

  Out <- list(Factors = FEVDFactors, Yields = FEVDYields)

  return(Out)
}

########################################################################################################
#' Compute GFEVDs for all models
#'
#' @param SIGMA Variance-covariance matrix
#' @param K1Z Loading As
#' @param G0 contemporaneous terms
#' @param BLoad Loading Bs
#' @param FactorLabels List containing the label of factors
#' @param FacDim Dimension of the P-dynamics
#' @param MatLength Length of the maturity vector
#' @param GFEVDhoriz Horizon of the analysis
#' @param YieldsLabel Label of bond yields
#' @param ModelType Desired model type
#' @param Economy specific economy under study
#' @param PI matrix PI for JLL-based models
#' @param Mode allows for the orthogonalized version in the case of JLL-based models
#'
#' @references
#' \itemize{
#' \item This function is a modified and extended version of the "fevd" function from
#' Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0, available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.
#'
#' \item Pesaran and Shin, 1998. "Generalized impulse response analysis in linear multivariate models" (Economics Letters)
#' }
#'
#' @keywords internal

ComputeGFEVDs <- function(SIGMA, K1Z, G0, BLoad, FactorLabels, FacDim, MatLength, GFEVDhoriz, YieldsLabel,
                          ModelType, Economy, PI = NULL, Mode = FALSE) {
  # 1) Dynamic multipliers
  Ry.h <- array(NA, c(FacDim, FacDim, GFEVDhoriz))
  Ry.h[, , 1] <- diag(FacDim) # dynamic multiplier at t=0

  for (l in 2:GFEVDhoriz) {
    Ry.h[, , l] <- K1Z %*% Ry.h[, , l - 1]
  }

  # 2) Initialization/ Minor preliminary work
  GFEVDresFac <- array(NA, c(nrow = FacDim, ncol = GFEVDhoriz, FacDim))
  vslct <- diag(FacDim)
  eslct <- diag(FacDim)

  invG <- diag(nrow(G0)) / G0
  invG[!is.finite(invG)] <- 0
  invGSigmau <- solve(G0) %*% SIGMA

  scale <- 1 / diag(SIGMA)

  # Adjustment term
  if (Mode == "Ortho" & any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    AdjTerm <- PI
  } else {
    AdjTerm <- diag(FacDim)
  }


  # 3) GFEVD
  # 3.1) Factors
  for (h in 1:FacDim) {
    n <- 1
    num <- matrix(0, nrow = FacDim, ncol = GFEVDhoriz)
    den <- matrix(0, nrow = FacDim, ncol = GFEVDhoriz)
    while (n <= GFEVDhoriz) {
      for (l in 1:n) {
        acc1 <- t((eslct[, h] %*% AdjTerm %*% Ry.h[, , l] %*% invGSigmau %*% vslct)^2) # Contribution of all j variables to explain i
        num[, n] <- num[, n] + acc1
        acc2 <- eslct[, h] %*% AdjTerm %*% Ry.h[, , l] %*% invGSigmau %*% t(invG) %*% t(Ry.h[, , l]) %*% eslct[, h]
        den[, n] <- den[, n] + matrix(1, nrow = FacDim) %*% acc2
      }
      GFEVDresFac[, n, h] <- t(t(scale * num[, n])) / den[, n]
      n <- n + 1
    }
  }

  GFEVDFactors <- aperm(GFEVDresFac, c(2, 1, 3)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

  #  Normalization of the GFEVD for the factors
  # (Make sure that the sum of the errors equal to one in each period)
  DEM <- array(NA, c(nrow = GFEVDhoriz, ncol = 1, FacDim))
  GFEVDFactorsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol = FacDim, FacDim))

  for (h in 1:FacDim) {
    for (n in 1:GFEVDhoriz) {
      DEM[n, 1, h] <- sum(GFEVDFactors[n, , h])
      GFEVDFactorsNormalized[n, , h] <- GFEVDFactors[n, , h] / DEM[n, , h]
    }
  }

  # 3.2) Yields
  # Initialization
  GFEVDresYie <- array(NA, c(nrow = MatLength, ncol = FacDim, GFEVDhoriz))
  vslctYie <- diag(FacDim)
  eslctYie <- diag(MatLength)

  num <- matrix(0, nrow = MatLength, ncol = FacDim)
  den <- matrix(0, nrow = MatLength, ncol = FacDim)

  for (l in 1:GFEVDhoriz) {
    acc1 <- (eslctYie %*% BLoad %*% AdjTerm %*% Ry.h[, , l] %*% invGSigmau %*% vslctYie)^2
    num <- num + acc1
    acc2 <- diag(eslctYie %*% BLoad %*% AdjTerm %*% Ry.h[, , l] %*% invGSigmau %*% t(invG) %*% t(Ry.h[, , l]) %*% t(AdjTerm) %*% t(BLoad) %*% eslctYie)
    den <- den + acc2
    for (q in 1:FacDim) {
      GFEVDresYie[, q, l] <- scale[q] * (num / den)[, q] # note: unlike the GFEVD of the factors, note that the "scale" variable is now at the acc1
    }
  }

  GFEVDYields <- aperm(GFEVDresYie, c(3, 2, 1)) # Non-normalized GFEVD (i.e. rows need not sum up to 1)

  #  Normalization of the GFEVD for the factors
  # (Make sure that the sum of the errors equal to one in each period)
  DEM <- array(NA, c(nrow = GFEVDhoriz, ncol = 1, MatLength))
  GFEVDYieldsNormalized <- array(NA, c(nrow = GFEVDhoriz, ncol = FacDim, MatLength))

  for (h in 1:MatLength) {
    for (n in 1:GFEVDhoriz) {
      DEM[n, 1, h] <- sum(GFEVDYields[n, , h])
      GFEVDYieldsNormalized[n, , h] <- GFEVDYields[n, , h] / DEM[n, , h]
    }
  }

  # 4) Prepare labels
  if (ModelType == "JPS original") {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economy]])
  } else if (any(ModelType == c("JPS global", "GVAR single"))) {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  } else {
    AllFactorsLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
  }

  # 5) Outputs to export
  Horiz <- 1:(GFEVDhoriz) # # We don't subtract 1, because there is no contemporaneous effect for the FORECAST error

  dimnames(GFEVDFactorsNormalized) <- list(Horiz, AllFactorsLabels, AllFactorsLabels)
  dimnames(GFEVDYieldsNormalized) <- list(Horiz, AllFactorsLabels, YieldsLabel)

  Out <- list(Factors = GFEVDFactorsNormalized, Yields = GFEVDYieldsNormalized)

  return(Out)
}

########################################################################################################
#' Fit yields for all maturities of interest
#'
#' @param MatInt numerical vector containing the fit maturities of interest
#' @param ModelPara List of model parameter estimates (See the "Optimization" function)
#' @param FactorLabels a string-list based which contains all the labels of all the variables present in the model
#' @param ModelType a string-vector containing the label of the model to be estimated
#' @param Economies a string-vector containing the names of the economies which are part of the economic system
#' @param YLab Label of yields ("Months" or "Yields")
#'
#' @keywords internal

YieldsFitAll <- function(MatInt, ModelPara, FactorLabels, ModelType, Economies, YLab) {
  # 0) Auxliary functions
  # a) Function to compute fitted yields
  compute_fitted_yields <- function(MatInt, MatAll, K1XQ, ModelType, r0, SSX, X, T_dim, dt, Economies,
                                    Time_Labels, Yield_Labels) {
    LoadingsLat <- Get__BnXAnX(MatAll, K1XQ, ModelType, r0, SSX, Economies)
    AnXAll <- LoadingsLat$AnX / dt
    BnXAll <- LoadingsLat$BnX / dt

    if (ModelType %in% Joint_Lab) {
      C <- length(Economies)
      FitLat <- matrix(NA, nrow = C * length(MatInt), ncol = T_dim)
    } else {
      FitLat <- matrix(NA, nrow = length(MatInt), ncol = T_dim)
    }

    IdxMatInt <- sort(unlist(lapply(MatInt, function(h) seq(h, length(AnXAll), by = max(MatAll)))))
    AnXInt <- AnXAll[IdxMatInt]
    BnXInt <- BnXAll[IdxMatInt, ]
    FitLat <- matrix(AnXInt, nrow = nrow(FitLat), ncol = T_dim) + BnXInt %*% X

    if (ModelType %in% Joint_Lab) {
      colnames(FitLat) <- Time_Labels
    } else {
      dimnames(FitLat) <- list(Yield_Labels, Time_Labels)
    }

    return(FitLat)
  }

  # b) Function to compute time series of the latent factors
  compute_X <- function(ZZ, Wpca, BnX, AnX, T_dim) {
    PP <- ZZ
    X <- matrix(NA, nrow = nrow(PP), ncol = T_dim)
    for (t in 1:T_dim) {
      X[, t] <- solve(Wpca %*% BnX, tol = 1e-50) %*% (PP[, t] - Wpca %*% AnX)
    }
    return(X)
  }

  # 1) Preliminar work
  Joint_Lab <- c("JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma")

  C <- length(Economies)
  Mod_ParaSet <- ModelPara[[ModelType]]
  N <- if (ModelType %in% Joint_Lab) Mod_ParaSet$Inputs$N else Mod_ParaSet[[Economies[1]]]$Inputs$N
  dt <- if (ModelType %in% Joint_Lab) Mod_ParaSet$Inputs$dt else Mod_ParaSet[[Economies[1]]]$Inputs$dt
  mat <- if (ModelType %in% Joint_Lab) Mod_ParaSet$Inputs$mat else Mod_ParaSet[[Economies[1]]]$Inputs$mat
  J <- length(mat)
  M <- length(FactorLabels$Domestic) - N
  G <- length(FactorLabels$Global)

  T_dim <- if (ModelType %in% Joint_Lab) {
    ncol(Mod_ParaSet$Inputs$AllFactors)
  } else {
    ncol(Mod_ParaSet[[Economies[1]]]$Inputs$AllFactors)
  }

  # 2) Compute results for joint or separate models
  # a) JointQ models
  if (ModelType %in% Joint_Lab) {
    BnX <- Mod_ParaSet$ModEst$Q$Load$X$B
    AnX <- Mod_ParaSet$ModEst$Q$Load$X$A
    K1XQ <- Mod_ParaSet$ModEst$Q$K1XQ
    SSX <- Mod_ParaSet$ModEst$Q$Load$X$SS
    r0 <- Mod_ParaSet$ModEst$Q$r0
    Wpca <- Mod_ParaSet$Inputs$Wpca
    ZZ <- Mod_ParaSet$Inputs$AllFactors

    b <- IdxSpanned(G, M, N, C)
    PP <- ZZ[b, ]
    X <- compute_X(PP, Wpca, BnX, AnX, T_dim)

    MatAll <- 1:max(mat / dt)
    Time_Labels <- colnames(Mod_ParaSet$Inputs$AllFactors)
    Yield_Labels <- paste(MatInt, YLab, sep = "")
    FittedYieldsPerMat <- compute_fitted_yields(
      MatInt, MatAll, K1XQ, ModelType, r0, SSX, X, T_dim, dt, Economies,
      Time_Labels, Yield_Labels
    )

    h <- length(MatInt)

    FittedYieldsCS <- lapply(1:C, function(i) {
      start_row <- (i - 1) * h + 1
      end_row <- i * h
      FittedYieldsPerMat[start_row:end_row, , drop = FALSE]
    })

    names(FittedYieldsCS) <- Economies

    return(FittedYieldsCS)
  } else {
    # b) SepQ models
    FittedYieldsPerMat <- list()

    for (i in 1:C) {
      BnX <- Mod_ParaSet[[Economies[i]]]$ModEst$Q$Load$X$B
      AnX <- Mod_ParaSet[[Economies[i]]]$ModEst$Q$Load$X$A
      K1XQ <- Mod_ParaSet[[Economies[i]]]$ModEst$Q$K1XQ
      SSX <- Mod_ParaSet[[Economies[i]]]$ModEst$Q$Load$X$SS
      r0 <- Mod_ParaSet[[Economies[i]]]$ModEst$Q$r0
      Wpca <- Mod_ParaSet[[Economies[i]]]$Inputs$Wpca
      ZZ <- Mod_ParaSet[[Economies[i]]]$Inputs$AllFactors

      if (ModelType == "JPS original") {
        AllLabels <- c(FactorLabels$Global, FactorLabels$Tables[[Economies[i]]])
      } else {
        AllLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
      }

      LabelSpannedCS <- c(FactorLabels$Tables[[Economies[i]]][-(1:M)])
      b <- match(LabelSpannedCS, AllLabels)
      PP <- ZZ[b, ]

      X <- compute_X(PP, Wpca, BnX, AnX, T_dim)
      MatAll <- 1:max(mat / dt)
      col_names <- colnames(Mod_ParaSet[[Economies[1]]]$Inputs$AllFactors)
      row_names <- paste(MatInt, YLab, sep = "")
      FittedYieldsPerMat[[i]] <- compute_fitted_yields(MatInt, MatAll, K1XQ, ModelType, r0, SSX, X, T_dim, dt, Economies, col_names, row_names)
    }

    names(FittedYieldsPerMat) <- Economies
    return(FittedYieldsPerMat)
  }
}

#################################################################################################
#' Adjust vector of maturities
#'
#' @param mat vector of maturities (J x 1)
#' @param UnitYields Available options: "Month" and "Year"
#'
#' @keywords internal

MatAdjusted <- function(mat, UnitYields) {
  if (UnitYields == "Month") {
    k <- 12
  } else if (UnitYields == "Year") {
    k <- 1
  }
  matAdjUnit <- mat * k

  return(matAdjUnit)
}
########################################################################################################
#' Get the expected component of all models
#'
#' @param ModelPara list of model parameter estimates
#' @param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#' @param ModelType desired model type
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param WishFP If users wants to compute the forward premia. Default is FALSE.
#'
#' @keywords internal

ExpectedComponent <- function(ModelPara, InputsForOutputs, ModelType, Economies, FactorLabels, WishFP = FALSE) {
  # 0) Preliminary work
  N <- length(FactorLabels$Spanned)
  UnitYields <- InputsForOutputs$UnitMatYields

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    mat <- ModelPara[[ModelType]][[1]]$Inputs$mat
  } else {
    mat <- ModelPara[[ModelType]]$Inputs$mat
  }

  matAdjUnit <- MatAdjusted(mat, UnitYields)

  if (WishFP) {
    matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
    matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]
  }

  # 1) Compute risk-neutral parameters
  rhos <- rhoParas(ModelPara, N, ModelType, Economies)

  # 2) Compute the expectation component
  # Pure expected component
  avexp <- Compute_EP(ModelPara, ModelType, UnitYields, matAdjUnit, N, rhos, Economies, FactorLabels)

  # expected component related to the forward premia (if desired)
  if (WishFP) {
    avexpFP <- Compute_EP(
      ModelPara, ModelType, UnitYields, matAdjUnit, N, rhos, Economies, FactorLabels,
      WishFP, matMIN, matMAX
    )
  } else {
    avexpFP <- list()
  }

  return(list(avexp = avexp, avexpFP = avexpFP))
}

#################################################################################################
#' Compute  risk-neutral intercept and slope
#'
#' @param ModelPara list of model parameter estimates
#' @param N number of country-specific spanned factors
#' @param ModelType desired model type
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

rhoParas <- function(ModelPara, N, ModelType, Economies) {
  # Compute the intercept and slope coefficients of the short rate expressed as a function of the spanned factors
  # By definition: r_t = r0 + rho1_X* X_t
  # But X_t = (W*BX)^(-1)(P- WAx)
  # so r_t = rho0_PP + rho1_PP*P_t
  # where (i) rho0_PP = r0 - rho1_X*(W*BX)^(-1)W*AX and (ii) rho1_PP = rho1_X (W*BX)^(-1)

  # 1) Models estimated for countries SEPARETLY
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    rho0_PP <- vector(mode = "list", length = length(Economies))
    rho1_PP <- vector(mode = "list", length = length(Economies))
    names(rho0_PP) <- Economies
    names(rho1_PP) <- Economies

    dt <- ModelPara[[ModelType]][[Economies[1]]]$Inputs$dt

    for (i in 1:length(Economies)) {
      Para_Set <- ModelPara[[ModelType]][[Economies[i]]]
      BnX <- Para_Set$ModEst$Q$Load$X$B
      AnX <- Para_Set$ModEst$Q$Load$X$A
      Wpca <- Para_Set$Inputs$Wpca
      r0 <- Para_Set$ModEst$Q$r0

      # Compute rhos
      rho1_X <- rep(1, N)
      rho0_PP[[i]] <- as.numeric((r0 - rho1_X %*% solve(Wpca %*% BnX, tol = 1e-50) %*% Wpca %*% AnX) / dt)
      rho1_PP[[i]] <- (rho1_X %*% solve(Wpca %*% BnX, tol = 1e-50)) / dt
    }
  } else {
    # 2) Models estimated for countries JOINTLY
    Para_Set <- ModelPara[[ModelType]]
    dt <- Para_Set$Inputs$dt
    BnX <- Para_Set$ModEst$Q$Load$X$B
    AnX <- Para_Set$ModEst$Q$Load$X$A
    Wpca <- Para_Set$Inputs$Wpca
    r0 <- Para_Set$ModEst$Q$r0

    # Compute rhos
    rho1_X_CS <- rep(1, N)
    rho1_X <- matrix(0, nrow = length(Economies), ncol = N * length(Economies))

    idx0 <- 0
    for (j in 1:length(Economies)) {
      idx1 <- idx0 + N
      rho1_X[j, (idx0 + 1):idx1] <- rho1_X_CS
      idx0 <- idx1
    }

    rho0_PP <- (r0 - rho1_X %*% solve(Wpca %*% BnX, tol = 1e-50) %*% Wpca %*% AnX) / dt
    rho1_PP <- (rho1_X %*% solve(Wpca %*% BnX, tol = 1e-50)) / dt
  }

  return(list(rho0_PP = rho0_PP, rho1_PP = rho1_PP))
}

###################################################################################################
#' Compute the expected component for all models
#'
#' @param ModelPara list of model parameter estimates
#' @param ModelType Desired model type
#' @param UnitYields Available options: "Month" and "Year"
#' @param matAdjUnit Adjusted vector of matutities
#' @param N number of country-specific spanned factors
#' @param rhoList List of risk-neutral parameters
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels List of factor labels
#' @param WishFP If users wants to compute the forward premia. Default is FALSE.
#' @param matMIN For the forward premia, the shortest maturity of the remium of interest
#' @param matMAX For the forward premia, the longest maturity of the remium of interest
#'
#' @keywords internal

Compute_EP <- function(ModelPara, ModelType, UnitYields, matAdjUnit, N, rhoList, Economies, FactorLabels,
                       WishFP = FALSE, matMIN = FALSE, matMAX = FALSE) {
  # 0) Preliminary work
  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")
  M <- length(FactorLabels$Domestic) - N

  dt <- if (ModelType %in% SepQ_Lab) {
    ModelPara[[ModelType]][[Economies[1]]]$Inputs$dt
  } else {
    ModelPara[[ModelType]]$Inputs$dt
  }

  k <- if (UnitYields == "Month") 12 else 1
  YLab <- if (UnitYields == "Month") "M" else "Y"
  J <- length(matAdjUnit)

  if (WishFP) { # Forward premia features
    EP_Lab <- "FP_"
    ExpecCompLength <- 1
    matMINAdj <- round((matMIN / k) / dt)
    matMAXAdj <- round((matMAX / k) / dt)
  } else {
    EP_Lab <- "RP_"
    ExpecCompLength <- round((matAdjUnit / k) / dt)
  }

  # 1) Compute the expected component

  ParaSet <- ModelPara[[ModelType]]
  # jointQ models
  if (!(ModelType %in% SepQ_Lab)) {
    K0Z <- ParaSet$ModEst$P$K0Z
    K1Z <- ParaSet$ModEst$P$K1Z
    ZZ <- ParaSet$Inputs$AllFactors
  }

  avexp <- list()
  for (i in 1:length(Economies)) {
    econ <- Economies[i]
    # a) Pre-allocation
    if (ModelType %in% SepQ_Lab) { # SepQ models
      K0Z <- ParaSet[[econ]]$ModEst$P$K0Z
      K1Z <- ParaSet[[econ]]$ModEst$P$K1Z
      ZZ <- ParaSet[[econ]]$Inputs$AllFactors
    }

    # b) Extract spanned factors from the list of unspanned factors
    IdxSpanned <- if (ModelType %in% SepQ_Lab) {
      all_labels <- if (ModelType == "JPS original") {
        c(FactorLabels$Global, FactorLabels$Tables[[econ]])
      } else {
        c(FactorLabels$Global, FactorLabels$Tables$AllCountries)
      }

      match(FactorLabels$Tables[[econ]][-(1:M)], all_labels)
    } else {
      IdxAllSpanned(ModelType, FactorLabels, Economies)
    }

    # c) Get the expected component
    avexpCS <- matrix(0, ncol(ZZ), length(ExpecCompLength))
    dimnames(avexpCS) <- list(colnames(ZZ), paste(EP_Lab, ExpecCompLength, YLab, sep = ""))

    for (h in 1:length(ExpecCompLength)) { # per bond maturity
      # Drop the 1-period maturity (If available), since by definition the expected component explains 100%
      # of the 1-period maturity yield
      if (h == 1 && matAdjUnit[1] == 1) {

        if (ModelType %in% SepQ_Lab) {
        avexpCS[ , 1] <- ModelPara[[ModelType]][[econ]]$Inputs$Y[1,]
        } else {
        idxJ <- (i-1)*J + 1
        avexpCS[ , 1]  <- ModelPara[[ModelType]]$Inputs$Y[idxJ,]
        }

    } else {

          for (t in 1:ncol(ZZ)) { # Per point in time

        # Initialization
        if (WishFP) {
          g <- matrix(NA, nrow(K0Z), matMAXAdj)
        } else {
          g <- matrix(NA, nrow(K0Z), ExpecCompLength[h])
        }
        rownames(g) <- rownames(ZZ)

        # Fitted P-dynamics
        g[, 1] <- ZZ[, t]
        if (WishFP) {
          loopLim <- matMAXAdj
        } else {
          loopLim <- ExpecCompLength[h]
        }
        for (j in 2:loopLim) {
          g[, j] <- K0Z + K1Z %*% g[, j - 1]
        }

        # Adjust `g_lim` dynamically
        g_lim <- if (WishFP) {
          round((matMIN / k) / dt):round((matMAX / k) / dt)
        } else {
          seq_len(ncol(g))
        }

        g <- g[IdxSpanned, g_lim] # extract relevant variables

        # Adjustment term
        if (ModelType %in% SepQ_Lab) {
          MaxExpec <- pmax(rhoList$rho0_PP[[i]] + (rhoList$rho1_PP[[i]] %*% g), 0)
        } else {
          MaxExpec <- pmax(rhoList$rho0_PP[i] + (rhoList$rho1_PP[i, ] %*% g), 0)
        }

        avexpCS[t, h] <- mean(MaxExpec)
      }
    }
    }
    avexp[[Economies[i]]] <- avexpCS * 100
  }

  return(avexp)
}

###################################################################################################
#' Compute the term premia
#'
#' @param ModelPara list of model parameter estimates
#' @param avexp list containing the country-specific pure expected component
#' @param ModelType desired model type
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

TermPremia <- function(ModelPara, avexp, ModelType, Economies) {
  # 0) Preliminary work
  YieldData <- list()
  TermPremium <- list()

  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")

  if (!ModelType %in% SepQ_Lab) {
    mat <- ModelPara[[ModelType]]$Inputs$mat
    J <- length(mat)
  }

  # 1) Compute the term premia
  for (i in 1:length(Economies)) {
    if (ModelType %in% SepQ_Lab) { # SepQ models
      Y <- ModelPara[[ModelType]][[Economies[i]]]$Inputs$Y
      YieldData[[Economies[i]]] <- t(Y) * 100
    } else { # jointQ models
      IdxRP <- (1:J) + J * (i - 1)
      Y <- ModelPara[[ModelType]]$Inputs$Y
      YieldData[[Economies[i]]] <- t(Y[IdxRP, ] * 100)
    }
    TermPremium[[Economies[i]]] <- YieldData[[Economies[i]]] - avexp[[Economies[i]]]
  }

  Output <- list(TermPremium, avexp)
  names(Output) <- c("Term Premia", "Expected Component")

  return(Output)
}
###############################################################################################
#' Compute the forward premia for all models
#'
#' @param ModelPara list of model parameter estimates
#' @param avexpFP list containing the country-specific expected component of the forward period
#' @param ModelType desired model type
#' @param FactorLabels List of factor labels
#' @param InputsForOutputs list containing the desired horizon of analysis for the model fit, IRFs, GIRFs, FEVDs, GFEVDs,
#'                        and risk premia decomposition
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

ForwardPremia <- function(ModelPara, avexpFP, ModelType, FactorLabels, InputsForOutputs, Economies) {
  # 0) Preliminary work: redefine necessary inputs
  SepQ_Lab <- c("JPS original", "JPS global", "GVAR single")

  if (ModelType %in% SepQ_Lab) {
    mat <- ModelPara[[ModelType]][[1]]$Inputs$mat
  } else {
    mat <- ModelPara[[ModelType]]$Inputs$mat
  }

  J <- length(mat)

  # Yield labels
  UnitYields <- InputsForOutputs$UnitMatYields
  if (UnitYields == "Month") {
    YLab <- "M"
  } else {
    YLab <- "Y"
  }

  # Inputs from the forward premia specification
  matAdjUnit <- MatAdjusted(mat, UnitYields)

  matMIN <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[1]
  matMAX <- InputsForOutputs[[ModelType]]$ForwardPremia$Limits[2]
  IDXMatMIN <- match(matMIN, matAdjUnit)
  IDXMatMAX <- match(matMAX, matAdjUnit)
  IDXMatBoth <- c(IDXMatMIN, IDXMatMAX)

  FR <- list()
  ForwardPremium <- list()

  # CASE 1: If any of the data of the specific maturity were not used in the estimation, then compute the fitted value
  if (anyNA(IDXMatBoth)) {
    MatBoth <- c(matMIN, matMAX)
    IdxNA <- which(is.na(IDXMatBoth))

    # a) Missing maturity (model-implied)
    MissingMat <- MatBoth[IdxNA] # Maturity not available
    YieldMissing <- YieldsFitAll(MissingMat, ModelPara, FactorLabels, ModelType, Economies, YLab)

    for (i in 1:length(Economies)) {
      if (ModelType %in% SepQ_Lab) {
        Y <- ModelPara[[ModelType]][[Economies[i]]]$Inputs$Y
      } else {
        Y <- ModelPara[[ModelType]]$Inputs$Y
      }

      # b) If both maturities are missing
      if (length(MissingMat) == 2) {
        YieldMIN <- t(t(YieldMissing[[Economies[i]]][1, ])) * 100
        YieldMAX <- t(t(YieldMissing[[Economies[i]]][2, ])) * 100
      } else {
        # c) If only one maturity is missing
        # Available maturity
        IdxNotMissing0 <- IDXMatBoth[!is.na(IDXMatBoth)]
        if (ModelType %in% SepQ_Lab) {
          IdxNotMissingCS <- IdxNotMissing0
        } else {
          IdxNotMissingCS <- IdxNotMissing0 + J * (i - 1)
        }

        YieldNotMissing <- t(t(Y[IdxNotMissingCS, ]))

        if (MissingMat == 1) {
          YieldMIN <- YieldNotMissing * 100
          YieldMAX <- t(YieldMissing[[Economies[i]]]) * 100
        } else {
          YieldMIN <- t(YieldMissing[[Economies[i]]]) * 100
          YieldMAX <- YieldNotMissing * 100
        }
      }

      FR[[Economies[i]]] <- (matMAX * YieldMAX - matMIN * YieldMIN) / (matMAX - matMIN) # Fitted forward rate
      ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
      colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep = "")
      colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep = "")
    }
  } else {
    # CASE 2: when all data is available
    YieldData <- list()

    for (i in 1:length(Economies)) {
      # SepQ models
      if (ModelType %in% SepQ_Lab) {
        Y <- ModelPara[[ModelType]][[Economies[i]]]$Inputs$Y
        YieldData[[Economies[i]]] <- t(Y) * 100

        YieldMIN <- t(t(YieldData[[Economies[i]]][, IDXMatMIN]))
        YieldMAX <- t(t(YieldData[[Economies[i]]][, IDXMatMAX]))
      } else {
        # jointQ models
        Y <- ModelPara[[ModelType]]$Inputs$Y

        IdxMinCS <- IDXMatMIN + J * (i - 1)
        IdxMaxCS <- IDXMatMAX + J * (i - 1)
        YieldMIN <- t(t(Y[IdxMinCS, ] * 100))
        YieldMAX <- t(t(Y[IdxMaxCS, ] * 100))
      }

      FR[[Economies[i]]] <- (matMAX * YieldMAX - matMIN * YieldMIN) / (matMAX - matMIN) # Fitted forward rate
      ForwardPremium[[Economies[i]]] <- FR[[Economies[i]]] - avexpFP[[Economies[i]]] # Forward Premia
      colnames(ForwardPremium[[Economies[i]]]) <- paste("FP_", matMIN, "-", matMAX, YLab, sep = "")
      colnames(FR[[Economies[i]]]) <- paste("Mat", matMIN, "-", matMAX, YLab, sep = "")
    }
  }

  OutputFP <- list(ForwardPremium, avexpFP, FR)
  names(OutputFP) <- c("Forward Premia", "Expected Component", "Forward Rate")

  return(OutputFP)
}

################################################################################################################
#' Find the indexes of the spanned factors
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param Economies  string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

IdxAllSpanned <- function(ModelType, FactorLabels, Economies) {
  G <- length(FactorLabels$Global)
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  IdxSpanned <- c()

  if (ModelType == "JPS original") {
    IdxSpanned <- (G + M + 1):(G + M + N)
  } else {
    idxSpa0 <- G + M
    for (j in 1:C) {
      idxSpa1 <- idxSpa0 + N

      if (j == 1) {
        IdxSpanned <- (idxSpa0 + 1):idxSpa1
      } else {
        IdxSpanned <- c(IdxSpanned, (idxSpa0 + 1):idxSpa1)
      }
      idxSpa0 <- idxSpa1 + M
    }
  }

  return(IdxSpanned)
}
