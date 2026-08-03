#' Builds the complete set of time series of the risk factors (spanned and unspanned)
#'
#' @param FactorSet Factor set list (see e.g. "GVARFactors" data file)
#' @param Economies A character vector containing the names of the economies included in the system.
#' @param FactorLabels A list of character vectors with labels for all variables in the model.
#' @param Initial_Date Start date of the sample period in the format yyyy-mm-dd
#' @param Final_Date   End date of the sample period in the format yyyy-mm-dd
#' @param DataFrequency A character vector specifying the data frequency. Available options: "Daily All Days", "Daily Business Days",
#'                      "Weekly", "Monthly", "Quarterly", "Annually".
#'
#'
#' @keywords internal
#'
#' @return
#' Risk factors used in the estimation of the desired ATSM
#'


RiskFactorsPrep <- function(FactorSet, Economies, FactorLabels, Initial_Date, Final_Date, DataFrequency) {
  Initial_Date <- as.Date(Initial_Date)
  Final_Date <- as.Date(Final_Date)

  G <- length(FactorSet$Global)
  C <- length(Economies) # We subtract 1 because we want to exclude the list of global factors
  N <- length(FactorLabels$Spanned)
  M <- length(FactorLabels$Domestic) - N
  T_dim <- length(FactorSet[[1]][[1]][[1]])


  ZZfull <- matrix(NA, nrow = C * (N + M) + G, ncol = T_dim)
  ZZfull[seq_len(G), ] <- do.call(rbind, lapply(FactorSet$Global, matrix, ncol = T_dim))

  # 1) Assign variables:
  idx0 <- G
  for (i in 1:C) {
    idx1 <- idx0 + N + M
    ZZfull[(idx0 + 1):idx1, ] <- do.call(rbind, lapply(FactorSet[[Economies[i]]]$Factors[1:(N + M)], matrix, ncol = T_dim))
    idx0 <- idx1
  }

  # 2) Assign labels:
  # a) Factor names
  idx0 <- 0
  LabCountries <- rep(NA, times = C * (N + M))
  for (i in 1:C) {
    idx1 <- idx0 + N + M
    LabCountries[(idx0 + 1):idx1] <- FactorLabels$Tables[[Economies[i]]]
    idx0 <- idx1
  }

  rownames(ZZfull) <- c(FactorLabels$Global, LabCountries)

  # b) Time series labels
  # (i) Daily All Days data
  if (DataFrequency == "Daily All Days") {
    DAD <- seq(Initial_Date, Final_Date, by = "1 day")
    DateLabels <- format(DAD, "%d-%m-%Y")
  }
  # (ii) Daily Business Days data
  if (DataFrequency == "Daily Business Days") {
    DateLabels <- 1:T_dim
  }
  # (iii) Weekly data
  if (DataFrequency == "Weekly") {
    week <- seq(Initial_Date, Final_Date, by = "1 week")
    DateLabels <- format(week, "%d-%m-%Y")
  }
  # (iv) Monthly data
  if (DataFrequency == "Monthly") {
    month <- seq(Initial_Date, Final_Date, by = "1 month")
    DateLabels <- format(month, "%d-%m-%Y")
  }
  # (v) Quarterly data
  if (DataFrequency == "Quarterly") {
    quarter <- seq(Initial_Date, Final_Date, by = "1 quarter")
    DateLabels <- format(quarter, "%m-%Y")
  }
  # (vi) Annually data
  if (DataFrequency == "Annually") {
    year <- seq(Initial_Date, Final_Date, by = "1 year")
    DateLabels <- format(year, "%Y")
  }


  colnames(ZZfull) <- DateLabels


  return(ZZfull)
}
