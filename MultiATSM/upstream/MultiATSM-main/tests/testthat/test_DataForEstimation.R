library(testthat)
library(MultiATSM)


#  Retrieve data from excel files
MacroAllData <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
YieldsAllData <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))

# Define common inputs of interest
DomVar <- c("Eco_Act", "Inflation")
GlobalVar <- c("GBC", "CPI_OECD")
t0 <- "2006-09-01"
tF <- "2011-09-01"
Economies <- c("China", "Brazil")
N <- 1

# 1) Test output structure for non-GVAR model
test_that("DataForEstimation returns correct structure for non-GVAR model", {
  ModelType <- "JPS original"
  K <- (N + length(DomVar)) * length(Economies) + length(GlobalVar)
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  DataFrequency <- "Monthly"
  res <- DataForEstimation(
    t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency,
    MacroAllData, YieldsAllData
  )
  expect_type(res, "list")
  expect_true(all(c("Yields", "RiskFactors", "GVARFactors") %in% names(res)))
  expect_null(res$GVARFactors)
  expect_equal(nrow(res$RiskFactors), K)
})


# 2) Test output structure for GVAR-based models
test_that("DataForEstimation returns correct structure for GVAR model", {
  ModelType <- "GVAR multi"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  DataFrequency <- "Monthly"
  #data("CM_Trade")
  res <- DataForEstimation(t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency, MacroAllData,
    YieldsAllData, TradeFlows,
    W_type = "Sample Mean",
    t_First_Wgvar = "2006", t_Last_Wgvar = "2006"
  )
  expect_type(res, "list")
  expect_true(all(c("Yields", "RiskFactors", "GVARFactors") %in% names(res)))
  expect_type(res$GVARFactors, "list")
  expect_equal(names(res$GVARFactors), c(Economies, "Global"))
})


# 3) Test error for invalid input
# a) Type 1
test_that("DataForEstimation throws error for invalid input", {
  Economies <- character(0)
  ModelType <- "JPS original"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  DataFrequency <- "Monthly"
  expect_error(DataForEstimation(
    t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency, MacroAllData,
    YieldsAllData
  ))
})

# b) Type 2
test_that("DataForEstimation throws error for invalid input", {
  Economies <- character(0)
  ModelType <- "JPS original"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  DataFrequency <- "Quarterly"
  expect_error(DataForEstimation(t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency))
})
