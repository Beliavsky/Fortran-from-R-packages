#' Numerical outputs (IRFs, GIRFs, FEVD, and GFEVD) for bootstrap
#'
#' @param ModelType A character vector indicating the model type to be estimated.
#' @param ModelParaBoot A list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#' @param InputsForOutputs A list containing the necessary inputs for generating IRFs, GIRFs, FEVDs, GFEVDs and Term Premia.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Economies A character vector containing the names of the economies included in the system.
#'
#' @keywords internal

NumOutputs_Bootstrap <- function(ModelType, ModelParaBoot, InputsForOutputs, FactorLabels, Economies) {
  # IRF and GIRF
  IRFout <- IRFandGIRF_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$IRF$horiz, FactorLabels, Economies)
  # FEVD and GFEVD
  FEVDout <- FEVDandGFEVD_BS(ModelType, ModelParaBoot, InputsForOutputs[[ModelType]]$FEVD$horiz, FactorLabels, Economies)

  # Output summary
  NumericalOutputs <- list(IRF = IRFout$IRFs, FEVD = FEVDout$FEVDs, GIRF = IRFout$GIRFs, GFEVD = FEVDout$GFEVDs)

  return(NumericalOutputs)
}

######################################################################################################
########################################### 1) IRFs and GIRFs ########################################
######################################################################################################
#' IRFs and GIRFs after bootstrap for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#' @param IRFhoriz single numerical vector containing the desired horizon of analysis for the IRFs
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

IRFandGIRF_BS <- function(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies) {
  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  # Pre-allocation
  GIRFoutputs <- list()
  IRFoutputs <- list()

  # 1) SINGLE COUNTRY MODELS
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (economy in Economies) {
      # a) Extract general parameters
      para_economy <- ModelParaBoot$ParaDraws[[ModelType]][[economy]]
      ndraws <- length(para_economy)
      K <- nrow(para_economy[[1]]$ModEst$P$K1Z)
      J <- length(ModelParaBoot$GeneralInputs$mat)
      YieldsLabel <- rownames(para_economy[[1]]$Inputs$Y)

      # b) Compute IRFs and GIRFs
      results <- lapply(seq_len(ndraws), function(tt) {
        # Extract parameters of each draw
        para_tt <- para_economy[[tt]]$ModEst$P
        SIGMA <- para_tt$SSZ # K x K (variance-covariance matrix)
        K1Z <- para_tt$K1Z # K x K (feedback matrix)
        G0.y <- para_tt$Gy.0

        B <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, economy, ModelType, tt)

        # Compute IRFs and GIRFs
        list(
          IRFs = ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, economy),
          GIRFs = ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, J, IRFhoriz, YieldsLabel, ModelType, economy)
        )
      })

      # c) Store results in respective lists
      IRFoutputs[[ModelType]][[economy]] <- lapply(results, `[[`, "IRFs")
      GIRFoutputs[[ModelType]][[economy]] <- lapply(results, `[[`, "GIRFs")
    }
  } else {
    # 2) JOINT COUNTRY MODELS
    JLL_label <- c("JLL original", "JLL No DomUnit", "JLL joint Sigma")
    ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])
    J <- length(ModelParaBoot$GeneralInputs$mat)
    K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[1]]$ModEst$P$K1Z)
    YieldsLabel <- rownames(ModelParaBoot$ParaDraws[[ModelType]][[1]]$Inputs$Y) # Yield labels

    for (tt in seq_len(ndraws)) {
      # a) Summarize inputs for the IRFs
      para_tt <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ModEst$P
      K1Z <- para_tt$K1Z # KxK (feedback matrix)
      G0.y <- para_tt$Gy.0
      SIGMA <- if (ModelType %in% JLL_label) {
        para_tt$JLLoutcomes$Sigmas$Sigma_Y
      } else {
        para_tt$SSZ
      }

      BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ModEst$Q$Load$P$B
      B <- BUnspannedAdapJoint(G, M, N, C, J, BSpanned)

      # b) Compute IRFs
      IRFoutputs[[ModelType]][[tt]] <- ComputeIRFs(SIGMA, K1Z, B, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType)

      # c) Compute GIRFs
      if (ModelType %in% JLL_label) {
        SIGMA <- para_tt$JLLoutcomes$Sigmas$VarCov_NonOrtho
      }
      GIRFoutputs[[ModelType]][[tt]] <- ComputeGIRFs(SIGMA, K1Z, B, G0.y, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType)

      # 3) JLL-BASED MODELS (orthogonalized outputs)
      if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
        # Summarize inputs for the IRFs
        K1Ze <- para_tt$JLLoutcomes$k1_e # KxK (feedback matrix)
        PI <- para_tt$JLLoutcomes$PI
        Se <- para_tt$JLLoutcomes$Sigmas$Sigma_Ye
        SIGMA_e <- para_tt$JLLoutcomes$Sigmas$VarCov_Ortho

        # Compute orthogonalized IRFs
        IRFOrtho <- ComputeIRFs(Se, K1Ze, B, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType,
          PI = PI, Mode = "Ortho"
        )
        # Store  IRFs
        IRFoutputs[[ModelType]][[tt]] <- list(
          Factors = list(NonOrtho = IRFoutputs[[ModelType]][[tt]]$Factors, Ortho = IRFOrtho$Factors),
          Yields  = list(NonOrtho = IRFoutputs[[ModelType]][[tt]]$Yields, Ortho = IRFOrtho$Yields)
        )

        # Compute orthogonalized GIRFs
        GIRFsOrtho <- ComputeGIRFs(SIGMA_e, K1Ze, B, G0.y, FactorLabels, K, C * J, IRFhoriz, YieldsLabel, ModelType,
          PI = PI, Mode = "Ortho"
        )
        # Store GIRFs
        GIRFoutputs[[ModelType]][[tt]] <- list(
          Factors = list(NonOrtho = GIRFoutputs[[ModelType]][[tt]]$Factors, Ortho = GIRFsOrtho$Factors),
          Yields  = list(NonOrtho = GIRFoutputs[[ModelType]][[tt]]$Yields, Ortho = GIRFsOrtho$Yields)
        )
      }
    }
  }

  Out <- list(IRFs = IRFoutputs, GIRFs = GIRFoutputs)
  return(Out)
}

#########################################################################################################
################################### 2) FEVD and GFEVDs ##################################################
#########################################################################################################
#' FEVDs and GFEVDs after bootstrap for all models
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#' @param FEVDhoriz single numerical vector containing the desired horizon of analysis for the FEVDs
#' @param FactorLabels string-list based which contains all the labels of all the variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#'
#' @keywords internal

FEVDandGFEVD_BS <- function(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies) {
  C <- length(Economies)
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N

  # Pre-allocation
  FEVDoutputs <- list()
  FEVDOrtho <- list()
  GFEVDoutputs <- list()
  GFEVDOrtho <- list()

  # 1) SINGLE COUNTRY MODELS
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    for (economy in Economies) {
      # a) Extract general parameters
      para_economy <- ModelParaBoot$ParaDraws[[ModelType]][[economy]]
      ndraws <- length(para_economy)
      K <- nrow(para_economy[[1]]$ModEst$P$K1Z)
      J <- length(ModelParaBoot$GeneralInputs$mat)
      YieldsLabel <- rownames(para_economy[[1]]$Inputs$Y)


      # b) Compute FEVDs and GFEVDs
      results <- lapply(seq_len(ndraws), function(tt) {
        # Extract parameters of each draw
        para_tt <- para_economy[[tt]]$ModEst$P
        SIGMA <- para_tt$SSZ # K x K (variance-covariance matrix)
        K1Z <- para_tt$K1Z # K x K (feedback matrix)
        G0 <- para_tt$Gy.0

        B <- BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, economy, ModelType, tt)

        # Compute FEVDs and GFEVDs
        list(
          FEVDs = ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, economy),
          GFEVDs = ComputeGFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, J, FEVDhoriz, YieldsLabel, ModelType, economy)
        )
      })

      # c) Store results in respective lists
      FEVDoutputs[[ModelType]][[economy]] <- lapply(results, `[[`, "FEVDs")
      GFEVDoutputs[[ModelType]][[economy]] <- lapply(results, `[[`, "GFEVDs")
    }
  } else {
    # 2) JOINT COUNTRY MODELS
    ndraws <- length(ModelParaBoot$ParaDraws[[ModelType]])
    J <- length(ModelParaBoot$GeneralInputs$mat)
    K <- nrow(ModelParaBoot$ParaDraws[[ModelType]][[1]]$ModEst$P$K1Z)
    YieldsLabel <- rownames(ModelParaBoot$ParaDraws[[ModelType]][[1]]$Inputs$Y) # Yield labels

    for (tt in 1:ndraws) {
      # a) Extract relevant inputs
      para_tt <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ModEst$P
      K1Z <- para_tt$K1Z
      SIGMA <- para_tt$SSZ
      G0 <- para_tt$Gy.0

      if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
        Chol_Fac_JLL <- para_tt$JLLoutcomes$Sigmas$Sigma_Y
      }

      BSpanned <- ModelParaBoot$ParaDraws[[ModelType]][[tt]]$ModEst$Q$Load$P$B
      B <- BUnspannedAdapJoint(G, M, N, C, J, BSpanned)

      # b) Compute FEVDs
      FEVDoutputs[[ModelType]][[tt]] <- ComputeFEVDs(SIGMA, K1Z, G0, B, FactorLabels, K, C * J, FEVDhoriz, YieldsLabel,
        ModelType,
        CholFac_JLL = Chol_Fac_JLL
      )

      # c) Compute GFEVDs
      GFEVDoutputs[[ModelType]][[tt]] <- ComputeGFEVDs(
        SIGMA, K1Z, G0, B, FactorLabels, K, C * J, FEVDhoriz,
        YieldsLabel, ModelType
      )

      # 3) JLL-BASED MODELS (orthogonalized outputs)
      if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
        # Extract relevant inputs
        K1Ze <- para_tt$JLLoutcomes$k1_e # KxK (feedback matrix)
        PI <- para_tt$JLLoutcomes$PI
        SIGMA_Ortho <- para_tt$JLLoutcomes$Sigmas$VarCov_Ortho
        Chol_Fac_JLL_Ortho <- para_tt$JLLoutcomes$Sigmas$Sigma_Ye

        # a) Compute FEVDs orthogonalized
        FEVDOrtho[[ModelType]][[tt]] <- ComputeFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C * J, FEVDhoriz,
          YieldsLabel, ModelType,
          CholFac_JLL = Chol_Fac_JLL_Ortho, PI = PI,
          Mode = "Ortho"
        )

        # Gather Outputs
        FEVDoutputs[[ModelType]][[tt]]$Factors <- list(
          NonOrtho = FEVDoutputs[[ModelType]][[tt]]$Factors,
          Ortho = FEVDOrtho[[ModelType]][[tt]]$Factors
        )
        FEVDoutputs[[ModelType]][[tt]]$Yields <- list(
          NonOrtho = FEVDoutputs[[ModelType]][[tt]]$Yields,
          Ortho = FEVDOrtho[[ModelType]][[tt]]$Yields
        )

        # b) Compute GFEVDs orthogonalized
        GFEVDOrtho[[ModelType]][[tt]] <- ComputeGFEVDs(SIGMA_Ortho, K1Ze, G0, B, FactorLabels, K, C * J, FEVDhoriz,
          YieldsLabel, ModelType,
          PI = PI, Mode = "Ortho"
        )

        # Gather Outputs
        GFEVDoutputs[[ModelType]][[tt]]$Factors <- list(
          NonOrtho = GFEVDoutputs[[ModelType]][[tt]]$Factors,
          Ortho = GFEVDOrtho[[ModelType]][[tt]]$Factors
        )
        GFEVDoutputs[[ModelType]][[tt]]$Yields <- list(
          NonOrtho = GFEVDoutputs[[ModelType]][[tt]]$Yields,
          Ortho = GFEVDOrtho[[ModelType]][[tt]]$Yields
        )
      }
    }
  }

  Out <- list(FEVDs = FEVDoutputs, GFEVDs = GFEVDoutputs)
  return(Out)
}

######################################################################################################
######################################## AUXILIARY FUNCTIONS #########################################
######################################################################################################
#' Obtain the full form of B unspanned for "sep Q" models within the bootstrap setting
#'
#' @param G number of global unspanned factors
#' @param M number of country-specific domestic unspanned factors
#' @param ModelParaBoot list of model parameter estimates (see the "Optimization" function) after a bootstrap draw
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param Economy string-vector containing the names of the economy under study
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param tt number of the bootstrap draw
#'
#' @keywords internal

BUnspannedAdapSep_BS <- function(G, M, ModelParaBoot, Economies, Economy, ModelType, tt) {
  C <- length(Economies)
  J <- length(ModelParaBoot$GeneralInputs$mat)
  i <- match(Economy, Economies)
  N <- ModelParaBoot$GeneralInputs$N

  # Extract required structure once
  Para_tt <- ModelParaBoot$ParaDraws[[ModelType]][[Economy]][[tt]]
  K <- nrow(Para_tt$ModEst$P$K1Z)
  BSpanned <- Para_tt$ModEst$Q$Load$P$B

  # Pre-allocate zero matrix
  BUnspanned <- matrix(0, nrow = J, ncol = K)

  if (ModelType == "JPS original") {
    BUnspanned[, (K - N + 1):K] <- BSpanned
  } else if (ModelType %in% c("JPS global", "GVAR single")) {
    idx0 <- G + M
    IDX <- lapply(seq_len(C), function(h) {
      idx1 <- idx0 + N
      result <- seq(idx0 + 1, idx1)
      idx0 <<- idx1 + M
      result
    })

    BUnspanned[, IDX[[i]]] <- BSpanned
  }

  return(BUnspanned)
}
