#' Computes the transition matrix required in the estimation of the GVAR model
#'
#' @param t_First character. Sample starting date (format: yyyy).
#' @param t_Last character. Sample ending date (format: yyyy).
#' @param Economies character vector. Names of the \code{C} economies included in the system.
#' @param type character. Method for computing interdependence. Possible options:
#' \itemize{
#'   \item \code{"Time-varying"}: Computes time-varying interdependence and returns weight matrices for each year.
#'   \item \code{"Sample Mean"}: Returns a single weight matrix with average weights over the sample period.
#'   \item Specific year (e.g., "1998", "2005"): Computes time-invariant interdependence for the specified year.
#' }
#' @param DataConnectedness list or data frame. Data used to compute the transition matrix (e.g., trade flows).
#'
#' @return matrix or list of matrices. Time-varying or time-invariant transition matrix depending on 'type'.
#'
#' @details
#' If there is missing data for any country in a particular year, the transition matrix will include only NAs.
#'
#' @examples
#' t_First <- "2006"
#' t_Last <- "2019"
#' Economies <- c("China", "Brazil", "Mexico", "Uruguay")
#' type <- "Sample Mean"
#' # Load data if Connectedness data from excel, otherwise use pre-saved data
#' GetExcelData <- FALSE
#'
#' if (GetExcelData) {
#'   if (!requireNamespace("readxl", quietly = TRUE)) {
#'     stop(
#'       "Please install package \"readxl\" to use this feature.",
#'       call. = FALSE
#'     )
#'     DataPath <- system.file("extdata", "TradeData.xlsx", package = "MultiATSM")
#'     tab_names_Trade <- readxl::excel_sheets(DataPath)
#'     list_all_Trade <- suppressMessages(lapply(tab_names_Trade, function(x) {
#'       readxl::read_excel(path = DataPath, sheet = x)
#'     }))
#'     names(list_all_Trade) <- tab_names_Trade
#'
#'     L <- length(list_all_Trade)
#'
#'     for (i in 1:L) {
#'       Countries <- list_all_Trade[[i]][[1]]
#'       list_all_Trade[[i]] <- as.data.frame(list_all_Trade[[i]][, -1])
#'       rownames(list_all_Trade[[i]]) <- Countries
#'     }
#'
#'     DataConnectedness <- list_all_Trade
#'   }
#' } else {
#'   data(TradeFlows)
#'   DataConnectedness <- TradeFlows
#' }
#'
#' W_mat <- Transition_Matrix(t_First, t_Last, Economies, type, DataConnectedness)
#'
#' @export

Transition_Matrix <- function(t_First, t_Last, Economies, type, DataConnectedness) {
  # 1) Pre-allocation of variables
  DataAdj <- lapply(DataConnectedness, function(x) as.data.frame(t(x)))
  C <- length(Economies)
  T_dim <- ncol(DataConnectedness[[1]])
  WgvarAllYears <- vector(mode = "list", length = T_dim)
  names(WgvarAllYears) <- colnames(DataConnectedness[[1]])
  num <- matrix(NA, nrow = C, ncol = C)
  dem <- c()
  Wyear <- matrix(NA, nrow = C, ncol = C)

  # 2) Generate the matrix of weights year-by-year
  for (k in 1:T_dim) {
    for (h in 1:C) {
      for (j in 1:C) {
        num[h, j] <- DataAdj[[Economies[h]]][[Economies[j]]][k]
      }
      dem[h] <- sum(num[h, ])
      Wyear[h, ] <- num[h, ] / dem[h]
    }
    if (any(is.na(Wyear))) { # If there is missing data for any country of the system for that particularly year, then return a matrix of NAs.
      Wyear <- matrix(NA, nrow = C, ncol = C)
    }
    WgvarAllYears[[k]] <- Wyear
  }

  # Desired output
  # a) Time-varying: Gathers the interdependence matrices from the full Sample
  if (type == "Time-varying") {
    Wgvar <- Filter(function(x) all(stats::complete.cases(x)), WgvarAllYears)
  }
  # b) Sample Mean
  else if (type == "Sample Mean") {
    Y_First <- substring(t_First, 1, 4)
    Y_Last <- substring(t_Last, 1, 4)
    TimeLable <- colnames(DataConnectedness[[1]])
    idx0 <- which(TimeLable == Y_First)
    idx1 <- which(TimeLable == Y_Last)

    WgvarSubSample <- WgvarAllYears[idx0:idx1]
    Wgvar <- apply(simplify2array(WgvarSubSample), 1:2, mean)

    colnames(Wgvar) <- Economies
    rownames(Wgvar) <- Economies
  }
  # c) Specific Year
  else {
    Wgvar <- WgvarAllYears[[type]]
    colnames(Wgvar) <- Economies
    rownames(Wgvar) <- Economies
  }

  return(Wgvar)
}
