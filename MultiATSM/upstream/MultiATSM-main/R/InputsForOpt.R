#' Generates inputs necessary to build the likelihood function for the ATSM model
#'
#' @param InitialSampleDate Start date of the sample period in the format "dd-mm-yyyy"
#' @param FinalSampleDate End date of the sample period in the format "dd-mm-yyyy"
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param Yields numerical matrix with time series of yields (\code{J x Td} or \code{CJ x Td})
#' @param GlobalMacro numerical matrix with time series of the global risk factors (\code{G x Td})
#' @param DomMacro numerical matrix with time series of the country-specific risk factors for all \code{C} countries ( \code{C X Td} or  \code{CM x Td})
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param DataFrequency character. Data frequency. Permissible choices: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".
#' @param GVARlist list. Inputs for GVAR model estimation. See details below.
#' @param JLLlist list. Inputs for JLL model estimation. See details below.
#' @param WishBRW logical. Whether to estimate the physical parameter model with bias correction (see \code{\link{Bias_Correc_VAR}}). Default is FALSE.
#' @param BRWlist list. Inputs for bias-corrected estimation.
#' @param UnitYields character. Maturity unit of yields. Permissible choices: "Month" or "Year". Default is "Month".
#' @param CheckInputs logical. Whether to perform a prior check on the consistency of the provided input list. Default is TRUE.
#' @param BS_Adj logical. Whether to adjust the global series for the sepQ models in the Bootstrap setting. Default is FALSE.
#' @param verbose logical. Print progress messages. Default is TRUE.
#'
#' @section Permissible options for GVARlist:
#' \itemize{
#'   \item \code{VARXtype}: "unconstrained" or "constrained"
#'   \item \code{W_type}: "Time-varying" or "Sample Mean"
#'   \item \code{t_First_Wgvar}, \code{t_Last_Wgvar}: year as character
#' }
#'
#' @section Permissible options for JLLlist:
#' \itemize{
#'   \item \code{DomUnit}: name of the dominant economy or \code{None}
#'   \item \code{WishSigmas}: TRUE (estimate variance-covariance matrices) or FALSE
#'   \item \code{SigmaNonOrtho}: NULL or \code{K x K} matrix
#' }
#'
#' @section Permissible options for BRWlist:
#' \itemize{
#'   \item \code{BiasCorrection}: TRUE (bias-corrected) or FALSE
#'   \item \code{flag_mean}: TRUE (mean) or FALSE (median)
#'   \item \code{gamma}: numeric adjustment parameter
#'   \item \code{N_iter}: number of iterations
#'   \item \code{N_burn}: number of burn-in iterations
#'   \item \code{B}: number of bootstrap samples
#'   \item \code{checkBRW}: TRUE or FALSE
#'   \item \code{B_check}: number of bootstrap samples for closeness check
#' }
#'
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td} model time series dimension.
#'   \item \code{C} number of countries in the system.
#'   \item \code{G} number of global unspanned factors.
#'   \item \code{M} number of country-specific unspanned factors.
#'   \item \code{K} total number of risk factors.
#'   \item \code{J} number of bond yields per country used in estimation.
#' }
#'
#' @importFrom pracma null
#' @return An object of class 'ATSMModelInputs' containing the necessary inputs for performing the model optimization.
#'
#' @section Available Methods:
#' - `print(object)`
#' - `summary(object)`
#'
#' @examples
#' \donttest{
#' # Example 1:
#' data(GlobalMacro)
#' data(DomMacro)
#' data(Yields)
#'
#' ModelType <- "JPS original"
#' Economies <- "Mexico"
#' t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
#' tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
#' N <- 3
#' GlobalVar <- c("Gl_Eco_Act") # Global Variables
#' DomVar <- c("Eco_Act") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' DataFreq <- "Monthly"
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
#'   FactorLabels, Economies, DataFreq,
#'   CheckInputs = FALSE, verbose = FALSE
#' )
#'
#' # Example 2:
#' LoadData("CM_2024")
#'
#' ModelType <- "GVAR multi"
#'
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
#' tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
#' N <- 2
#' GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
#' DomVar <- c("Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' DataFreq <- "Monthly"
#' GVARlist <- list(
#'   VARXtype = "unconstrained", W_type = "Sample Mean",
#'   t_First_Wgvar = "2007", t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
#' )
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
#'   FactorLabels, Economies, DataFreq, GVARlist,
#'   CheckInputs = FALSE, verbose = FALSE
#' )
#'
#' # Example 3:
#' LoadData("CM_2024")
#'
#' ModelType <- "JLL original"
#'
#' Economies <- c("China", "Brazil", "Uruguay")
#' t0 <- "01-05-2007" # InitialSampleDate (Format: "dd-mm-yyyy")
#' tF <- "01-12-2019" # FinalSampleDate (Format: "dd-mm-yyyy")
#' N <- 2
#' GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation") # Global Variables
#' DomVar <- c("Eco_Act", "Inflation") # Domestic Variables
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#'
#' JLLinputs <- list(DomUnit = "China")
#'
#' DataFrequency <- "Monthly"
#'
#' ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
#'   FactorLabels, Economies, DataFreq,
#'   JLLlist = JLLinputs,
#'   CheckInputs = FALSE, verbose = FALSE
#' )
#' }
#' @export

InputsForOpt <- function(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
                         FactorLabels, Economies, DataFrequency, GVARlist = NULL, JLLlist = NULL,
                         WishBRW = FALSE, BRWlist = NULL, UnitYields = "Month", CheckInputs = TRUE,
                         BS_Adj = FALSE, verbose = TRUE) {
  # Print initial message with model type and sample period
  if (verbose) {
    message((paste(
      "1) PREPARING INPUTS FOR THE ESTIMATION OF THE MODEL:", ModelType, ". SAMPLE PERIOD:", InitialSampleDate,
      "-", FinalSampleDate
    )))
  }

  # Define labels for single and multi-country models
  Label_Single_Models <- c("JPS original", "JPS global", "GVAR single")
  Label_Multi_Models <- c("GVAR multi", "JPS multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma")

  # Check consistency of inputs across models if CheckInputs is TRUE
  if (isTRUE(CheckInputs)) {
    CheckInputsForMLE(
      InitialSampleDate, FinalSampleDate, Economies, DomMacro, GlobalMacro, UnitYields,
      DataFrequency, Label_Single_Models, Label_Multi_Models, FactorLabels, GVARlist, ModelType
    )
  }

  # Construct the time-series of the risk factors
  if (verbose) message(("1.1) Constructing the time-series of the risk factors"))
  RiskFactors <- BuildATSM_RiskFactors(
    InitialSampleDate, FinalSampleDate, Yields, GlobalMacro, DomMacro,
    Economies, FactorLabels, ModelType, BS_Adj
  )

  Yields <- AdjustYieldsDates(Yields, RiskFactors, Economies)

  # Build the vector of common maturities across countries
  mat <- Maturities(Yields, Economies, UnitYields)

  # Collect general input lists
  ModelInputsGen <- GeneralMLEInputs(
    Yields, RiskFactors, FactorLabels, mat, DataFrequency, Label_Multi_Models,
    Economies, ModelType, UnitYields
  )

  # Collect model-specific input lists (GVARinputs, JLLinputs, and BRWlist)
  ModelInputsSpe <- SpecificMLEInputs(
    ModelType, Economies, RiskFactors, FactorLabels, GVARlist, JLLlist,
    WishBRW, BRWlist
  )

  # Estimate P-dynamic parameters
  if (verbose) message(("1.2) Estimating model P-dynamics parameters:"))
  PdynPara <- GetPdynPara(
    RiskFactors, FactorLabels, Economies, ModelType, ModelInputsSpe$BRWinputs,
    ModelInputsSpe$GVARinputs, ModelInputsSpe$JLLinputs, CheckInputs, verbose
  )

  # Gather outputs to export
  Outputs <- Outputs2exportMLE(
    Label_Multi_Models, Economies, RiskFactors, Yields, mat, ModelInputsGen,
    ModelInputsSpe, PdynPara, ModelType
  )

  # Store metadata inside the class without explicitly exporting it
  attr(Outputs, "ModelInfo") <- list(
    ModelType = ModelType,
    Economies = Economies,
    InitialSampleDate = InitialSampleDate,
    FinalSampleDate = FinalSampleDate,
    RiskFactors = RiskFactors,
    Yields = Yields,
    DataFrequency = DataFrequency,
    Maturities = mat,
    G = length(FactorLabels$Global),
    NC = length(FactorLabels$Spanned) * length(Economies),
    MC = (length(FactorLabels$Domestic) - length(FactorLabels$Spanned)) * length(Economies)
  )

  # Return the structured Outputs object
  return(structure(Outputs, class = "ATSMModelInputs"))
}
##########################################################################################################################
#' Create a vector of numerical maturities in years
#'
#' @param DataYields matrix containing all yields of the system (JxT,if the model is single-country-based
#'                  or CJxT if the model is multy-country-based)
#' @param Economies vector containing names of all the economies of the system
#' @param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'
#' @return Vector containing all observed maturities expressed in years
#'
#' @keywords internal


Maturities <- function(DataYields, Economies, UnitYields) {
  if (UnitYields == "Month") {
    fac <- 12
  } else if (UnitYields == "Year") {
    fac <- 1
  }

  C <- length(Economies)
  s <- rownames(DataYields)
  AllMat <- as.numeric(gsub("[^0-9.-]", "", s))
  J <- length(AllMat) / C
  mat <- AllMat[1:J] / fac

  return(mat)
}

#########################################################################################################################################
#' Concatenate the model-specific inputs in a list
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param Economies string-vector containing the names of the economies of the system
#' @param RiskFactors time series of risk factors (F x T)
#' @param FactorLabels string-list based which contains the labels of all the variables present in the model
#' @param GVARlist A list of required inputs to estimate the GVAR-based setups:
#' \enumerate{
#'      \item VARXtype  string-vector containing the VARX feature (see "GVAR" function) (GVAR-based models)
#'      \item t_First_Wgvar Sample starting date (year) (GVAR-based models)
#'      \item t_Last_Wgvar  Sample last date (year) (GVAR-based models)
#'      \item W_type  Criterion used in the computation of the star variables (see "Transition_Matrix" function)
#'               (GVAR-based models)
#'  }
#' @param JLLlist A list of required inputs to estimate the JLL-based setups:
#' \enumerate{
#'      \item DomUnit name of the economy which is assigned as the dominant unit (JLL-based models)
#'      \item WishSigmas equal to "1" if one wishes the variance-covariance matrices and the Cholesky factorizations (JLL-based models)
#'      \item SigmaNonOrtho NULL or some F x F matrix from the non-orthogonalized dynamics (JLL-based models)
#' }
#' @param WishBRW Whether the user wishes to estimate the physical parameter model with the Bias correction model from BRW (2012) (see "Bias_Correc_VAR" function).\cr
#'              Default is set to 0.
#' @param BRWlist A list of required inputs to estimate the bias corrected setups of the type of BRW:
#' \enumerate{
#'      \item BiasCorrection binary variable. it takes value equal to 1 if the user whishes the estimates to be bias-corrected
#'                      and 0, otherwise. (BRW model)
#'      \item flag_mean flag whether mean- (TRUE) or median- (FALSE) unbiased estimation is desired
#'      \item gamma adjustment parameter (BRW model)
#'      \item N_iter number of iterations (BRW model)
#'      \item N_burn number of burn-in iterations (BRW model)
#'      \item B  number of bootstrap samples (BRW model)
#'      \item checkBRW flag whether the user wishes to perform the closeness check (BRW model)
#'      \item B_check number of bootstrap samples for closeness check
#' }
#'
#' @keywords internal


SpecificMLEInputs <- function(ModelType, Economies, RiskFactors, FactorLabels, GVARlist = NULL, JLLlist = NULL,
                              WishBRW = 0, BRWlist = NULL) {
  if (ModelType %in% c("GVAR single", "GVAR multi")) {
    # Compute the transition matrix and the
    # PRELIMINARY CHECK:  Makes sure that link matrices are correctly imputed in a model with time-varying interdependence
    if (GVARlist$W_type == "Time-varying" & GVARlist$t_First_Wgvar != GVARlist$t_Last_Wgvar) {
      stop("For estimating GVAR models with time-varying interdependence, the start and ending dates
            of the transition matrix must coincide!")
    }

    t0 <- GVARlist$t_First_Wgvar
    tF <- GVARlist$t_Last_Wgvar
    LinkageMeasure <- GVARlist$DataConnectedness
    W_type <- GVARlist$W_type
    Wgvar <- Transition_Matrix(t0, tF, Economies, W_type, LinkageMeasure)
    if (ModelType == "GVAR single") {
      RiskFactors <- RiskFactors[[1]]
    }

    GVARFactors <- DataSet_BS(ModelType, RiskFactors, Wgvar, Economies, FactorLabels)

    # Build the list of the necessary inputs to estimate a GVAR model:
    GVARinputs <- list(Economies = Economies, GVARFactors = GVARFactors, VARXtype = GVARlist$VARXtype)

    if (W_type == "Time-varying") {
      GVARinputs$Wgvar <- Wgvar[[t0]]
    } else {
      GVARinputs$Wgvar <- Wgvar
    } # constant interdependence
  } else {
    GVARinputs <- NULL
  }


  if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    JLLinputs <- list(
      Economies = Economies,
      DomUnit = JLLlist$DomUnit,
      WishSigmas = 1, # Sigma matrix is estimated within the "InputsForMLEdensity" function
      SigmaNonOrtho = NULL,
      JLLModelType = ModelType
    )
  } else {
    JLLinputs <- NULL
  }


  if (WishBRW == 1) {
    BRWinputs <- BRWlist
  } else {
    BRWinputs <- NULL
  }

  outputs <- list(GVARinputs = GVARinputs, JLLinputs = JLLinputs, BRWinputs = BRWinputs)
  return(outputs)
}

###########################################################################################################
#' Makes sure that the time series of yields and risk factors have coincident sample spans
#'
#' @param Yields time series of bond yields (CJ x T1 )
#' @param PdynamicsFactors time series of risk factors (F x T2)
#' @param Economies string-vector containing the names of the economies of the system.
#'
#' @keywords internal


AdjustYieldsDates <- function(Yields, PdynamicsFactors, Economies) {
  # Restrict the sample span
  if (is.list(PdynamicsFactors)) {
    t0 <- colnames(PdynamicsFactors[[1]])[1]
    tF <- colnames(PdynamicsFactors[[1]])[ncol(PdynamicsFactors[[1]])]
  } else {
    t0 <- colnames(PdynamicsFactors)[1]
    tF <- colnames(PdynamicsFactors)[ncol(PdynamicsFactors)]
  }

  Idx0 <- which(colnames(Yields) == t0)
  IdxF <- which(colnames(Yields) == tF)

  YieldsAdj <- Yields[, Idx0:IdxF]

  # Restrict the sample of countries
  IdxCountries <- sapply(Economies, function(economy) grep(economy, rownames(YieldsAdj)))
  YieldsAdj <- YieldsAdj[IdxCountries, ]

  return(YieldsAdj)
}



###############################################################################################################
#################################################################################################################
#' Loads data sets from several papers

#' @param DataPaper  Available options are \code{BR_2017} (Bauer and Rudebusch, 2017) , \code{CM_2023} (Candelon and Moura, 2023), \code{CM_2024} (Candelon and Moura, 2024)
#'
#'
#' @examples
#' # Example 1:
#' LoadData("BR_2017")
#'
#' # Example 2:
#' LoadData("CM_2023")
#'
#' # Example 3:
#' LoadData("CM_2024")
#'
#' @returns
#' Complete set of data from several papers.
#'
#' @references
#' \enumerate{
#' \item Bauer and Rudebusch (2017). "Resolving the Spanning Puzzle in Macro-Finance Term Structure Models" (Review of Finance)
#' \item Candelon and Moura (2023). "Sovereign yield curves and the COVID-19 in emerging markets" (Economic Modelling)
#' \item Candelon and Moura (2024). "A Multicountry Model of the Term Structures of Interest Rates with a GVAR" (Journal of Financial Econometrics)
#' }
#' @export


LoadData <- function(DataPaper) {
  switch(DataPaper,
    "BR_2017" = utils::data("BR_jps_out"),
    "CM_2024" = {
      utils::data("GlobalMacro")
      utils::data("DomMacro")
      utils::data("TradeFlows")
      utils::data("Yields")
    },
    "CM_2023" = {
      utils::data("GlobalMacro_covid")
      utils::data("DomMacro_covid")
      utils::data("TradeFlows_covid")
      utils::data("Yields_covid")
    },
    stop("Database unavailable.")
  )
}


#################################################################################################################
#' Builds the time series of the risk factors that are used in the estimation of the ATSM
#'
#' @param InitialSampleDate Sample starting date
#' @param FinalSampleDate Sample last date
#' @param Yields matrix  (J x T), where  J - the number of maturities and  T - time series length
#' @param GlobalMacroFactors time series of the global macroeconomic risk factors (G x T)
#' @param DomMacroFactors time series of the country-specific macroeconomic risk factors for all C countries (CM x T)
#' @param Economies      string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BS_Adj adjustment of global series for sepQ model in the Bootstrap setting. Default is set to FALSE.
#'
#' @return
#' Generates the complete set of risk factors that are used in the estimation of the ATSM
#'
#' @keywords internal


BuildATSM_RiskFactors <- function(InitialSampleDate, FinalSampleDate, Yields, GlobalMacroFactors,
                                  DomMacroFactors, Economies, FactorLabels, ModelType, BS_Adj = FALSE) {
  # Preliminary work
  N <- length(FactorLabels$Spanned)
  G <- length(FactorLabels$Global)
  M <- length(FactorLabels$Domestic) - N
  C <- length(Economies)

  AllFactorLabels <- c(FactorLabels$Global, FactorLabels$Tables$AllCountries)

  F_dim <- C * (M + N) + G
  T_dim <- ncol(DomMacroFactors)

  RiskFactors <- matrix(NA, nrow = F_dim, ncol = T_dim)
  colnames(RiskFactors) <- colnames(DomMacroFactors)
  rownames(RiskFactors) <- AllFactorLabels

  # Input the global factors
  if (!BS_Adj || (ModelType %in% c("JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
    IdxGlobal <- match(FactorLabels$Global, rownames(GlobalMacroFactors))
    IdxGlobalRF <- IdxGlobal[!is.na(IdxGlobal)]

    CommonLabelsGlob <- intersect(rownames(GlobalMacroFactors), AllFactorLabels)
    GlobalMacroFactorsClean <- GlobalMacroFactors[rownames(GlobalMacroFactors) %in% CommonLabelsGlob, ]

    RiskFactors[seq_len(G), ] <- GlobalMacroFactorsClean
  }
  # Input the spanned factors
  PP <- Spanned_Factors(Yields, Economies, N)
  IdxSpa <- match(rownames(PP), AllFactorLabels)
  RiskFactors[IdxSpa, ] <- PP

  # Input the unspanned factors
  IdxUnspa <- match(rownames(DomMacroFactors), AllFactorLabels)
  IdxUnspaRF <- IdxUnspa[!is.na(IdxUnspa)]

  CommonLabelsUnsp <- intersect(rownames(DomMacroFactors), AllFactorLabels)
  DomMacroFactorsClean <- DomMacroFactors[rownames(DomMacroFactors) %in% CommonLabelsUnsp, ]

  RiskFactors[IdxUnspaRF, ] <- DomMacroFactorsClean

  # Adjust dates to the desired sample sample
  Idx0 <- which(colnames(RiskFactors) == InitialSampleDate)
  IdxF <- which(colnames(RiskFactors) == FinalSampleDate)

  RiskFactors <- RiskFactors[, Idx0:IdxF]

  # Adapt to the specifies of the single-country-based models
  if (ModelType %in% c("JPS original", "JPS global", "GVAR single")) {
    RiskFactorsList <- list()

    for (i in 1:C) {
      if (BS_Adj) {
        RiskFactors[seq_len(G), ] <- GlobalMacroFactors[((1 + G * (i - 1)):(G + G * (i - 1))), ]
      }

      if (ModelType == "JPS original") {
        idxCountryFactors <- which(grepl(Economies[i], rownames(RiskFactors))) # index of the rows of the country-specific factors
        idxFactors <- c(seq_len(G), idxCountryFactors) # index global + country-specific factors
        RiskFactorsList[[Economies[i]]] <- RiskFactors[idxFactors, ]
      } else {
        RiskFactorsList[[Economies[i]]] <- RiskFactors
      }
    }
    RiskFactors <- RiskFactorsList
  }

  return(RiskFactors)
}

###################################################################################################################
#' Get delta t
#'
#' @param DataFrequency single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#'
#' @keywords internal

Getdt <- function(DataFrequency) {
  if (DataFrequency == "Daily All Days") {
    dt <- 1 / 365
  } else if (DataFrequency == "Daily Business Days") {
    dt <- 1 / 252
  } else if (DataFrequency == "Weekly") {
    dt <- 1 / 52
  } else if (DataFrequency == "Monthly") {
    dt <- 1 / 12
  } else if (DataFrequency == "Quarterly") {
    dt <- 1 / 4
  } else if (DataFrequency == "Annually") {
    dt <- 1
  }

  return(dt)
}

#############################################################################################################
#' Check consistence of inputs
#'
#' @param t0 Sample starting date
#' @param tF Sample last date
#' @param Economies string-vector containing the names of the economies of the system.
#' @param DomesticMacroFac time series of the country-specific macroeconomic risk factors for all C countries (CM x T)
#' @param GlobalMacroFac time series of the global macroeconomic risk factors (G x T)
#' @param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#' @param DataFreq single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#' @param Label_Single_Models string-vector containing the names of the single country setups
#' @param Label_Multi_Models string-vector containing the names of the multicountry setups
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param GVARlist list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal

CheckInputsForMLE <- function(t0, tF, Economies, DomesticMacroFac, GlobalMacroFac, UnitYields, DataFreq,
                              Label_Single_Models, Label_Multi_Models, FactorLabels, GVARlist, ModelType) {
  N <- length(FactorLabels$Spanned)
  DomVarLab <- utils::head(FactorLabels$Domestic, -N)
  GloVarLab <- FactorLabels$Global
  TimeSeriesLabels <- colnames(DomesticMacroFac)

  # CHECK 1: Check whether the model choice is compatible with the number of countries selected
  if (length(Economies) == 1 & (any(ModelType == Label_Multi_Models))) {
    stop("The models 'GVAR multi', 'JPS multi', 'JLL original', 'JLL No DomUnit', and 'JLL joint Sigma' are multicountry setups.
       Therefore, they require the estimation of several countries.")
  }

  # CHECK 2: Check whether the model type selected if available
  if (!(ModelType %in% c(Label_Single_Models, Label_Multi_Models))) {
    stop("Model type is invalid.")
  }

  # CHECK 3: Label consistency in terms of economies, domestic and global sets
  Check_label_consistency(Economies, DomesticMacroFac, GlobalMacroFac, DomVarLab, GloVarLab)

  # CHECK 4: consistency of the chosen data sample
  if (!(all(c(t0, tF) %in% TimeSeriesLabels))) {
    stop("Initial and/or final sample dates are not part of the data sample span.")
  }

  # CHECK 5: data frequency chosen
  if (!(DataFreq %in% c("Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"))) {
    stop(paste("Data frequency option", DataFreq, "is unavailable. Available options are: 'Daily All Days', 'Daily Business Days', 'Weekly', 'Monthly',  'Quarterly', 'Annually'."))
  }
  # CHECK 6: unit of bond yields chosen
  if (!(UnitYields %in% c("Month", "Year"))) {
    stop(paste("Bond yield time-unit option", UnitYields, "is unavailable. Available options are 'Month' or 'Year'."))
  }
}
############################################################################################################
#' Gathers the general inputs for model estimation
#'
#' @param Yields matrix (CJ x T) or a list containing those matrices, where C is the number of countries,
#'              J - the number of maturities and T - time series length \cr
#' @param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param mat vector of maturities (in years) used in the estimation
#' @param DataFrequency single element character-based vector. Available options are: "Daily All Days", \cr
#'                      "Daily Business Days", "Weekly", "Monthly",  "Quarterly", "Annually"
#' @param Label_Multi_Models string-vector containing the names of the multicountry setups
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#'
#' @keywords internal


GeneralMLEInputs <- function(Yields, RiskFactors, FactorLabels, mat, DataFrequency, Label_Multi_Models,
                             Economies, ModelType, UnitYields) {
  N <- length(FactorLabels$Spanned)
  dt <- Getdt(DataFrequency)
  J <- length(mat)

  # Compute the quantities of interest
  idxJ0 <- 0
  idxN0 <- 0

  ListInputs <- list()

  for (i in 1:length(Economies)) { # Country-specific inputs
    idxJ1 <- idxJ0 + J
    idxN1 <- idxN0 + N

    YCS <- Yields[(idxJ0 + 1):idxJ1, ] # Yields (JxT)
    WpcaCS <- pca_weights_one_country(YCS, Economy = Economies[i])[1:N, ] # matrix of weights for the portfolio without errors (N x J)
    WpcaCS <- 100 * matrix(WpcaCS, nrow = N) # useful command for the case N = 1
    WeCS <- t(pracma::null(matrix(WpcaCS, nrow = N))) # matrix of weights for the yield portfolios priced with errors
    WpcaFullCS <- rbind(WpcaCS, WeCS)
    PPCS <- Spanned_Factors(YCS, Economies = Economies[i], N) # Spanned Factors (N x T)
    K1XQCS <- FeedMat_Q(YCS, PPCS, Economies[i], UnitYields, dt)


    if (any(ModelType == Label_Multi_Models)) {
      if (i == 1) {
        K1XQ <- K1XQCS
        Wpca <- WpcaCS
        WpcaFull <- WpcaFullCS
        We <- WeCS
        PP <- PPCS
        Y <- YCS
      } else {
        K1XQ <- adiag(K1XQ, K1XQCS)
        Wpca <- adiag(Wpca, WpcaCS)
        We <- adiag(We, WeCS)
        WpcaFull <- adiag(WpcaFull, WpcaFullCS)
        PP <- rbind(PP, PPCS)
        Y <- rbind(Y, YCS)
      }
    } else {
      # Matrix of contemporaneous terms
      Gy.0 <- diag(nrow(RiskFactors[[1]]))
      # Store the outputs of interest for the single country model types
      ListInputs[[Economies[i]]] <- list(
        Wpca = WpcaCS, We = WeCS, WpcaFull = WpcaFullCS, Yields = YCS,
        SpaFact = PPCS, K1XQ = K1XQCS, Gy.0 = Gy.0
      )
    }

    idxJ0 <- idxJ1
    idxN0 <- idxN1
  }

  # Matrix of contemporaneous terms
  if (any(ModelType == Label_Multi_Models)) {
    colnames(Y) <- colnames(RiskFactors)
    Gy.0 <- diag(nrow(RiskFactors))
  }

  # Store the outputs of interest for the multicountry model types
  if (any(ModelType == Label_Multi_Models)) {
    ListInputs <- list(Wpca = Wpca, We = We, WpcaFull = WpcaFull, Yields = Y, SpaFact = PP, K1XQ = K1XQ, Gy.0 = Gy.0)
  }

  return(ListInputs)
}

####################################################################################################
#' Get an estimate for the risk-neutral (Q) feedback matrix
#'
#' @param Yields matrix of bond yields (J x T)
#' @param Spa_Fac matrix of spanned factors  (N x T)
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param UnitYields (i) "Month": if maturity of yields are expressed in months or
#'                  (ii) "Year": if maturity of yields are expressed in years
#'
#' @param time_step time unit of the model (scalar). For instance, if data is (i) monthly, dt <- 12; (ii) quarterly, dt <- 4; (iii) yearly, dt <- 1.
#' @param check_inputs Perform input validation. Default is TRUE.
#'
#' @references
#' Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
#' Estimation of Some Term Structure Models. EABCN Training School.\cr
#' This function offers an independent R implementation that is informed
#' by the conceptual framework outlined in Le and Singleton (2018), but adapted to the
#' present modeling context.
#'
#' @keywords internal

FeedMat_Q <- function(Yields, Spa_Fac, Economies, UnitYields, time_step, check_inputs = TRUE) {
  # 1) Input Validation
  mat_vec <- Maturities(Yields, Economies, UnitYields)

  if (check_inputs) {
    CheckInput_K1X(Yields, mat_vec, time_step)
  }

  # 2) Get interpolated yield set
  Y_Inta <- Intra_Yields(Yields, mat_vec, time_step)

  # 3) Regression: Y = Betas × P
  Betas <- Reg_demean(Y_Inta, Spa_Fac)

  # 4) Estimate K1Q at horizon h
  min_mat <- mat_vec[1]
  max_mat <- nrow(Y_Inta)
  Min_horiz <- max(1, round(min_mat / time_step)) # Ensure horiz >= 1
  K1_h <- Est_K1h(Min_horiz, min_mat, max_mat, time_step, Betas)

  # 5) Compute matrix h-th root via eigen-decomposition
  N <- nrow(Spa_Fac)
  K1_root <- RootEigen(K1_h, Min_horiz, N)
  CheckNumericalPrecision(K1_root) # Check for numerical issues

  # 6) Convert K1XQ to the Jordan form
  K1XQ <- JordanMat(K1_root)

  return(K1XQ)
}

###############################################################################
#' Input validation for the 'FeedMat_Q' function
#'
#' @param Yields matrix containing country-specific yields  (J x T)
#' @param mat_vec vector of maturities expressed in years (J x 1)
#' @param time_step time unit of the model (scalar).
#'
#' @keywords internal

CheckInput_K1X <- function(Yields, mat_vec, time_step) {
  if (!is.numeric(time_step) || length(time_step) != 1 || time_step <= 0) {
    stop("Argument 'time_step' must be a positive numeric scalar.")
  }

  if (nrow(Yields) != length(mat_vec)) {
    stop("Length of 'maturity_vec' must equal the number of rows in yield dataset.")
  }
}
###############################################################################
#' Fit the cross-section of yields using spline
#'
#' @param Yields matrix containing country-specific yields  (J x T)
#' @param mat_vec vector of maturities expressed in years (J x 1)
#' @param time_step time unit of the model (scalar)
#'
#' @keywords internal

Intra_Yields <- function(Yields, mat_vec, time_step) {
  max_mat <- max(mat_vec)
  min_mat <- min(mat_vec)
  grid_points <- seq(
    from = time_step,
    to = ceiling(max_mat / time_step) * time_step,
    by = time_step
  )
  n_grid <- length(grid_points)

  interp_fn <- function(mats, yields, target_mats) {
    spline_func <- stats::splinefun(mats, yields, method = "fmm")
    spline_func(target_mats)
  }
  interpolated_Yields <- apply(Yields, 2, function(y) interp_fn(mat_vec, y, grid_points))

  return(interpolated_Yields)
}
############################################################################
#' Perform a linear regression using demeaned variables
#'
#' @param LHS Left-hand side variables
#' @param RHS Right-hand side variables
#'
#' @return Beta coefficients
#'
#' @keywords internal

Reg_demean <- function(LHS, RHS) {
  lhs_deman <- scale(t(LHS), center = TRUE, scale = FALSE)
  rhs_deman <- scale(t(RHS), center = TRUE, scale = FALSE)
  B_coef <- stats::lm(lhs_deman ~ rhs_deman - 1)$coefficients
  B_coef <- t(B_coef)

  return(B_coef)
}
############################################################################
#' Estimate K1h
#'
#' @param horiz Minimum horizon observed
#' @param min_mat shortest maturity observed
#' @param max_mat longest maturity observed
#' @param time_step time unit of the model (scalar)
#' @param B_coef matrix of beta coefficients
#'
#' @keywords internal

Est_K1h <- function(horiz, min_mat, max_mat, time_step, B_coef) {
  n_indices <- (2 * horiz):max_mat
  if (length(n_indices) == 0) {
    stop("Not enough grid points for estimation. Increase time_step or data range.")
  }

  LHS <- RHS <- matrix(NA, nrow = length(n_indices), ncol = ncol(B_coef))
  for (i in seq_along(n_indices)) {
    n <- n_indices[i]
    LHS[i, ] <- B_coef[n, ] - (horiz / n) * B_coef[horiz, ]
    RHS[i, ] <- (1 - horiz / n) * B_coef[n - horiz, ]
  }

  # Solve for K1h
  K1h <- tryCatch(
    {
      qr.solve(RHS, LHS)
    },
    error = function(e) {
      stop("Failed to solve for K1Qh. System may be ill-conditioned: ", e$message)
    }
  )

  return(K1h)
}

############################################################################
#' Compute the root of the eigenvalue of K1h
#'
#' @param K1_h K1_h squared matrix
#' @param horiz Minimum horizon observed
#' @param N number of country-specific spanned factors
#'
#' @keywords internal

RootEigen <- function(K1_h, horiz, N) {
  eig <- eigen(K1_h)
  root_eigenvalues <- Mod(eig$values)^(1 / horiz) * exp(1i * Arg(eig$values) / horiz)
  if (N == 1) root_eigenvalues <- matrix(root_eigenvalues)

  K1Q_est <- Re(eig$vectors %*% diag(root_eigenvalues) %*% solve(eig$vectors))

  if (any(is.nan(K1Q_est)) || any(is.infinite(K1Q_est))) {
    stop("Numerical instability detected in the K1Q matrix")
  }

  return(K1Q_est)
}
############################################################################
#' Check Numerical Precision Issues of K1_root matrix
#'
#' @param matrix_in K1_root matrix
#' @param tolerance Numerical tolerance
#' @return List with precision indicators
#'
#' @keywords internal
CheckNumericalPrecision <- function(matrix_in, tolerance = 1e-10) {
  diagnostics <- list(
    is_precise = TRUE,
    has_nan = any(is.nan(matrix_in)),
    has_inf = any(is.infinite(matrix_in)),
    has_large_values = any(abs(matrix_in) > 1e6)
  )

  # Check for NaN/Inf
  if (diagnostics$has_nan || diagnostics$has_inf) {
    stop("Matrix contains NaN or Inf values")
    diagnostics$is_precise <- FALSE
  }

  # Check for extreme values that might indicate numerical issues
  if (diagnostics$has_large_values) {
    stop("Matrix contains unplausibly large values - possible numerical overflow")
    diagnostics$is_precise <- FALSE
  }

  return(diagnostics)
}

############################################################################
#' Convert a Matrix to Jordan-Like Form for Term Structure Models
#'
#' @param matrix_in squared matrix prior to Jordan form
#'
#' @return A matrix in a specialized block form used in term structure modeling
#'
#' @references
#' \itemize{
#'   \item The mathematical approach is based on methods described in:
#'   Dai, Q., & Singleton, K. J. (2000). Specification Analysis of Affine
#'   Term Structure Models. The Journal of Finance, 55(5), 1943-1978.
#'
#'   \item Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
#'   Estimation of Some Term Structure Models. EABCN Training School.
#'
#'   \item For theoretical background on Jordan forms in term structure models:
#'   Duffee, G. R. (2002). Term Premia and Interest Rate Forecasts in Affine
#'   Models. The Journal of Finance, 57(1), 405-443.
#' }
#' @keywords internal

JordanMat <- function(matrix_in) {
  # Eigenvalues of the input matrix
  eigvals <- eigen(matrix_in)$values

  # Separate real and complex eigenvalues
  is_real <- Im(eigvals) == 0
  real_vals <- Re(eigvals[is_real])
  complex_vals <- eigvals[!is_real]

  # Build the "lQ" vector from eigenvalue pairs
  lq <- numeric(0)

  # Handle odd real eigenvalue case
  if (length(real_vals) %% 2 == 1) {
    lq <- real_vals[1]
    real_vals <- real_vals[-1]
  }

  # Process real eigenvalue pairs
  if (length(real_vals) > 0) {
    for (h in seq_len(length(real_vals) / 2)) {
      avg <- 0.5 * (real_vals[2 * h - 1] + real_vals[2 * h])
      diff_sq <- (0.5 * (real_vals[2 * h - 1] - real_vals[2 * h]))^2
      lq <- c(lq, avg, diff_sq)
    }
  }

  # Process complex eigenvalue pairs (assumed conjugates)
  if (length(complex_vals) > 0) {
    for (h in seq_len(length(complex_vals) / 2)) {
      real_part <- Re(complex_vals[2 * h - 1])
      imag_sq <- -abs(Im(complex_vals[2 * h - 1]))^2
      lq <- c(lq, real_part, imag_sq)
    }
  }

  n <- length(lq)

  # Build Jordan block skeleton
  if (n == 1) {
    jordan_mat <- matrix(0, 1, 1)
  } else {
    jordan_mat <- rbind(rep(0, n), diag(1, n - 1, n))
  }

  # Fill in diagonal and superdiagonal blocks
  offset <- 0
  if (n %% 2 == 1) {
    jordan_mat[1, 1] <- lq[1]
    lq <- lq[-1]
    offset <- 1
  }

  if (length(lq) > 0) {
    for (h in seq_len(length(lq) / 2)) {
      idx <- offset + 2 * h - 1
      jordan_mat[idx, idx] <- lq[2 * h - 1]
      jordan_mat[idx + 1, idx + 1] <- lq[2 * h - 1]
      jordan_mat[idx, idx + 1] <- lq[2 * h]
    }
  }

  return(jordan_mat)
}

##############################################################################################################
#' Compute the parameters used in the P-dynamics of the model
#'
#' @param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BRWinputs list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#' @param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#' @param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#' @param CheckInputs Logical. Whether to perform a prior check on the consistency of the provided input list. Default is FALSE.
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal


GetPdynPara <- function(RiskFactors, FactorLabels, Economies, ModelType, BRWinputs, GVARinputs,
                        JLLinputs, CheckInputs = F, verbose) {
  N <- length(FactorLabels$Spanned)

  # Parameters of the of the P-dynamics
  # Get bias-corrected parameters, if selected
  if (!is.null(BRWinputs)) {
    if (verbose) message("- With the bias-correction procedure from BRW (2012)")

    if (CheckInputs) {
      if (any(ModelType == c("GVAR single", "GVAR multi"))) {
        CheckInputsGVAR(GVARinputs, N)
      } else if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
        CheckJLLinputs(RiskFactors, JLLinputs)
      }
    }

    tryCatch(
      {
        PdynPara <- GetPdynPara_BC(
          ModelType, BRWinputs, RiskFactors, Economies, FactorLabels, GVARinputs,
          JLLinputs, verbose
        )
      },
      error = function(err) {
        stop("BRW procedure leads to a highly collinear system. Please choose another combination of BRW parameters.")
      }
    )
  } else {
    if (verbose) message("- Without the bias-correction procedure \n")
    if (any(ModelType == c("JLL original", "JLL No DomUnit", "JLL joint Sigma"))) {
      if (verbose) {
        message("JLL-based setup in progress: Estimating the variance-covariance matrix numerically.
                             This may take some time. \n")
      }
    }
    # Otherwise compute the non-biased correction model parameters
    PdynPara <- GetPdynPara_NoBC(ModelType, RiskFactors, Economies, N, GVARinputs, JLLinputs, CheckInputs)
  }

  return(PdynPara)
}

#####################################################################################################
#' Prepares inputs to export
#'
#' @param Label_Multi_Models string-vector containing the names of the multicountry setups
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#' @param Yields matrix (CJ x T) or a list containing those matrices, where C is the number of countries,
#'              J - the number of maturities and T - time series length \cr
#' @param mat vector of maturities (in years) used in the estimation
#' @param ModelInputsGen List of generic inputs
#' @param ModelInputsSpec List of specific inputs
#' @param PdynPara Model parameters estimated in the P-dynamics the
#' @param ModelType string-vector containing the label of the model to be estimated
#'
#' @keywords internal


Outputs2exportMLE <- function(Label_Multi_Models, Economies, RiskFactors, Yields, mat, ModelInputsGen,
                              ModelInputsSpec, PdynPara, ModelType) {
  if (any(ModelType == Label_Multi_Models)) {
    Output <- list(
      mat = mat,
      Wpca = ModelInputsGen$Wpca,
      We = ModelInputsGen$We,
      WpcaFull = ModelInputsGen$WpcaFull,
      Yields = ModelInputsGen$Y,
      SpaFact = ModelInputsGen$SpaFact,
      RiskFactors = RiskFactors,
      Gy.0 = ModelInputsGen$Gy.0,
      K1XQ = ModelInputsGen$K1XQ,
      SSZ = PdynPara$SSZ,
      K0Z = PdynPara$K0Z,
      K1Z = PdynPara$K1Z,
      JLLinputs = ModelInputsSpec$JLLinputs,
      GVARinputs = ModelInputsSpec$GVARinputs
    )
  } else {
    Output <- stats::setNames(
      lapply(Economies, function(econ) {
        list(
          mat = mat,
          Wpca = ModelInputsGen[[econ]]$Wpca,
          We = ModelInputsGen[[econ]]$We,
          WpcaFull = ModelInputsGen[[econ]]$WpcaFull,
          Yields = ModelInputsGen[[econ]]$Y,
          SpaFact = ModelInputsGen[[econ]]$SpaFact,
          RiskFactors = RiskFactors[[econ]],
          Gy.0 = ModelInputsGen[[econ]]$Gy.0,
          K1XQ = ModelInputsGen[[econ]]$K1XQ,
          SSZ = PdynPara[[econ]]$SSZ,
          K0Z = PdynPara[[econ]]$K0Z,
          K1Z = PdynPara[[econ]]$K1Z,
          JLLinputs = ModelInputsSpec$JLLinputs,
          GVARinputs = ModelInputsSpec$GVARinputs
        )
      }), Economies
    )
  }
  return(Output)
}
##############################################################################################################
#' Compute P-dynamics parameters using the bias correction method from BRW (2012)
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param BRWinputs list of necessary inputs for performing the bias-corrected estimation (see "Bias_Correc_VAR" function)
#' @param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param FactorLabels string-list based which contains the labels of all variables present in the model
#' @param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#' @param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#' @param verbose Logical flag controlling function messaging.
#'
#' @keywords internal


GetPdynPara_BC <- function(ModelType, BRWinputs, RiskFactors, Economies, FactorLabels, GVARinputs, JLLinputs, verbose) {
  C <- length(Economies)

  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    PdynPara <- list()
    for (i in 1:C) {
      RiskFactors_CS <- RiskFactors[[Economies[i]]]

      if (verbose) message(paste("Estimation for country:", Economies[i]))

      BiasCorrec <- Bias_Correc_VAR(ModelType, BRWinputs, t(RiskFactors_CS), Economies, FactorLabels,
        GVARinputs, JLLinputs,
        verbose = verbose
      )

      PdynPara[[Economies[i]]] <- list(
        K0Z = BiasCorrec$K0Z_BC, K1Z = BiasCorrec$K1Z_BC,
        SSZ = BiasCorrec$SSZ_BC
      )
    }
  } else {
    BiasCorrec <- Bias_Correc_VAR(ModelType, BRWinputs, t(RiskFactors), Economies, FactorLabels,
      GVARinputs, JLLinputs,
      verbose = verbose
    )
    K0Z <- BiasCorrec$K0Z_BC
    K1Z <- BiasCorrec$K1Z_BC
    SSZ <- BiasCorrec$SSZ_BC
  }

  # Export outputs
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Out <- PdynPara
  } else {
    Out <- list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ)
  }

  return(Out)
}


####################################################################################################################
#' Compute P-dynamics parameters without using the bias correction method from BRW (2012)
#'
#' @param ModelType string-vector containing the label of the model to be estimated
#' @param RiskFactors time series of risk factors (F x T). Could be stored in a list depending on the model
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param N number of country-specific spanned factors
#' @param GVARinputs list of necessary inputs for the estimation of GVAR-based models (see "GVAR" function)
#' @param JLLinputs list of necessary inputs for the estimation of JLL-based models (see "JLL" function)
#'
#' @keywords internal

GetPdynPara_NoBC <- function(ModelType, RiskFactors, Economies, N, GVARinputs, JLLinputs, CheckInpts = F) {
  # JPS-related models
  if (any(ModelType == c("JPS original", "JPS global"))) {
    PdynPara <- list()
    for (i in 1:length(Economies)) {
      RiskFactorsMat <- RiskFactors[[Economies[i]]]
      VARpara <- VAR(RiskFactorsMat, VARtype = "unconstrained")
      PdynPara[[Economies[i]]] <- VARpara
    }
  } else if (ModelType == "JPS multi") {
    VARpara <- VAR(RiskFactors, VARtype = "unconstrained")
    K0Z <- VARpara$K0Z
    K1Z <- VARpara$K1Z
    SSZ <- VARpara$SSZ
  }
  # GVAR-related models
  else if (any(ModelType == c("GVAR single", "GVAR multi"))) {
    GVARpara <- GVAR(GVARinputs, N, CheckInpts)
    K0Z <- GVARpara$F0
    K1Z <- GVARpara$F1
    SSZ <- GVARpara$Sigma_y

    if (ModelType == "GVAR single") {
      PdynPara <- list()
      for (i in 1:length(Economies)) {
        PdynPara[[Economies[i]]] <- list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ)
      }
    }
  }

  # JLL-related models
  else if (ModelType %in% c("JLL original", "JLL No DomUnit", "JLL joint Sigma")) {
    JLLinputs$WishSigmas <- 1
    JLLPara <- JLL(RiskFactors, N, JLLinputs, CheckInpts)
    K0Z <- JLLPara$k0
    K1Z <- JLLPara$k1
    JLLinputs$WishSigmas <- 0 # Ensures that the variance-covariance matrix will no longer be estimated within the JLL function

    if (any(ModelType == c("JLL original", "JLL No DomUnit"))) {
      SSZ <- JLLPara$Sigmas$VarCov_NonOrtho
    } else {
      SSZ <- JLLPara$Sigmas$VarCov_Ortho
      # NOTE: the required vectorization that precedes the numerical optimization of the variance-covariance matrix
      # is done on the orthogonalized variance-covariance matrix because we want to preserve the restrictions imposed in this matrix.
      # However, to compute the loadings A and B from Y= A + B*P we have to use the non-orthogonalized variance-covariance
      # matrix. We make this adjustment in the function 'MLEdensity" by redefining SSZ before computing A and B.
      # (this avoids overcomplicating the overall structure of the code")
    }
  }

  # Export outputs
  if (any(ModelType == c("JPS original", "JPS global", "GVAR single"))) {
    Out <- PdynPara
  } else {
    Out <- list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ)
  }

  return(Out)
}
###################################################################################################################
#' Check consistency of labels (economies, domestic and global variables)
#'
#' @param Economies string-vector containing the names of the economies which are part of the economic system
#' @param DomesticMacroFac time series of the country-specific domestic risk factors (C(M+N) x T)
#' @param GlobalMacroFac time series of the global risk factors (G x T)
#' @param DomVarLab string-vector containing the names of the desired domestic risk factors
#' @param GloVarLab string-vector containing the names of the desired global risk factors
#'
#' @keywords internal


Check_label_consistency <- function(Economies, DomesticMacroFac, GlobalMacroFac, DomVarLab, GloVarLab) {
  DomVarLabels_ALL <- row.names(DomesticMacroFac)
  GlobalVarLabels_ALL <- row.names(GlobalMacroFac)

  # a) Check consistence in the economy set
  missing_economies <- Economies[!sapply(Economies, function(economy) any(grepl(economy, DomVarLabels_ALL)))]
  if (length(missing_economies) > 0) {
    print(paste(
      "The following economy names are misspecified:", paste(missing_economies, collapse = ", "),
      ". Check whether (i) the label of the domestic variables is followed by the economy's name or (ii) the economy name is correctly spelled."
    ))
  }

  # b) Check consistence in the domestic variable set
  missing_DomVariables <- DomVarLab[!sapply(DomVarLab, function(DomVarLab) any(grepl(DomVarLab, DomVarLabels_ALL)))]
  if (length(missing_DomVariables) > 0) {
    print(paste("Inconsitent domestic variable name(s):", paste(missing_DomVariables, collapse = ", ")))
  }

  # c) Check consistence in the global variable set
  missing_GlobalVariables <- GloVarLab[!sapply(GloVarLab, function(GloVarLab) any(grepl(GloVarLab, GlobalVarLabels_ALL)))]
  if (length(missing_GlobalVariables) > 0) {
    print(paste("Inconsitent global variable name(s):", paste(missing_GlobalVariables, collapse = ", ")))
  }

  if (any(lengths(list(missing_economies, missing_DomVariables, missing_GlobalVariables)) > 0)) {
    stop("Inconsistent model input specifications. Check messages above.")
  }
}
