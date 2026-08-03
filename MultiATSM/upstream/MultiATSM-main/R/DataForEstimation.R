#' Retrieves data from Excel and builds the database used in the model estimation
#'
#' @param t0 character. Start date of the sample period in the format yyyy-mm-dd.
#' @param tF character. End date of the sample period in the format yyyy-mm-dd.
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param N positive integer. Number of country-specific spanned factors.
#' @param FactorLabels list. Labels for all variables present in the model, as returned by \code{\link{LabFac}}.
#' @param ModelType character. Model type to be estimated. Permissible choices: "JPS original", "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
#' @param DataFrequency character. Data frequency. Permissible choices: "Daily All Days", "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".
#' @param Macro_FullData list. Full set of macroeconomic data.
#' @param Yields_FullData list. Full set of bond yield data.
#' @param DataConnect list. Data for computing bilateral connectedness measures. Default is NULL. Required for GVAR-based models.
#' @param W_type character. Weight matrix type. Permissible choices: "Full Sample" (all years), "Sample Mean" (average over sample), or a specific year (e.g. "1998", "2005"). Default is NULL.
#' @param t_First_Wgvar character. First year for weight matrix computation. Default is NULL.
#' @param t_Last_Wgvar character. Last year for weight matrix computation. Default is NULL.
#'
#' @return A list containing:
#'   \enumerate{
#'     \item Yields: matrix (\code{J x Td} or \code{CJ x Td}) of bond yields for all countries.
#'     \item RiskFactors: matrix (\code{K x Td}) of risk factors for all countries.
#'     \item GVARFactors: list of variables used in VARX estimation (see \code{GVARFactors} data file). NULL if not GVAR-based.
#'   }
#'
#' @section General Notation:
#'   \itemize{
#'     \item \code{Td}: model time series dimension.
#'     \item \code{C}: number of countries in the system.
#'     \item \code{N}: number of country-specific spanned factors.
#'     \item \code{K}: total number of risk factors.
#'     \item \code{J}: number of bond yields per country used in estimation.
#'   }
#'
#' @seealso \code{\link{Load_Excel_Data}}
#'
#' @examples
#' DomVar <- c("Eco_Act", "Inflation")
#' GlobalVar <- c("GBC", "CPI_OECD")
#' t0 <- "2006-09-01"
#' tF <- "2019-01-01"
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay", "Russia")
#' N <- 2
#' ModelType <- "JPS original"
#' FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
#' DataFrequency <- "Monthly"
#' MacroData <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
#' YieldData <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
#' DataModel <- DataForEstimation(
#'   t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency,
#'   MacroData, YieldData
#' )
#'
#' @export

DataForEstimation <- function(t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency, Macro_FullData,
                              Yields_FullData, DataConnect = NULL, W_type = NULL, t_First_Wgvar = NULL,
                              t_Last_Wgvar = NULL) {
  C <- length(Economies)

  # Compute the transition matrix if the model is GVAR-based
  if (any(ModelType == c("GVAR single", "GVAR multi"))) {
    Wgvar <- Transition_Matrix(t_First_Wgvar, t_Last_Wgvar, Economies, W_type, DataConnect)
  }

  # Prepare the database
  FactorSet <- DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelType, Macro_FullData, Yields_FullData, Wgvar)

  # Gather all bond yields of all countries
  Yields <- c()
  for (i in 1:C) {
    if (i == 1) {
      Yields <- FactorSet[[Economies[i]]]$Yields
    } else {
      Yields <- rbind(Yields, FactorSet[[Economies[i]]]$Yields)
    }
  }

  # Gather the spanned and unspanned factors
  RiskFactors <- RiskFactorsPrep(FactorSet, Economies, FactorLabels, t0, tF, DataFrequency)

  # Gather the factors used in the GVAR estimation
  if (any(ModelType == c("GVAR single", "GVAR multi"))) {
    GVARFactors <- DataSet_BS(ModelType, RiskFactors, Wgvar, Economies, FactorLabels)
  } else {
    GVARFactors <- NULL
  }

  Outputs <- list(Yields, RiskFactors, GVARFactors)
  names(Outputs) <- c("Yields", "RiskFactors", "GVARFactors")

  return(Outputs)
}

##########################################################################################################
#' Read data from Excel files and return a named list of data frames
#'
#' @param ExcelFilePath character. Path to the Excel file (.xlsx) to load. Must be a valid file path. The file can contain multiple sheets; each sheet will be loaded as a separate data frame in the output list.
#'
#' @return Named list of data frames, one for each sheet in the Excel file. The names of the list elements correspond to the sheet names.
#'
#' @details
#' Uses the readxl package to read all sheets from the specified Excel file. Each sheet is returned as a data frame. The output is a named list, with names matching the sheet names in the Excel file.
#'
#' @examples
#' if (!requireNamespace("readxl", quietly = TRUE)) {
#'   stop(
#'     "Please install package \"readxl\" to use this feature.",
#'     call. = FALSE
#'   )
#'
#'   Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
#'   Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
#' }
#' @export

Load_Excel_Data <- function(ExcelFilePath) {
  sheets <- readxl::excel_sheets(ExcelFilePath)
  out <- lapply(sheets, function(x) readxl::read_excel(ExcelFilePath, sheet = x))
  names(out) <- sheets
  out
}
