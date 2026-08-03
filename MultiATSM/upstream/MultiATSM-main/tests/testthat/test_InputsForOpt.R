library(testthat)
library(MultiATSM)

# Define inputs
LoadData("CM_2024")
InitialSampleDate <- "01-05-2007"
FinalSampleDate <- "01-12-2018"
DataFrequency <- "Monthly"
GlobalVar <- c("Gl_Eco_Act")
DomVar <- c("Eco_Act", "Inflation")
N <- 2


# 1) Test output structure (Model type "JPS original")
# a) No bias correction
test_that("InputsForOpt returns correct output structure (JPS model)", {
  ModelType <- "JPS original"
  Economies <- "Mexico"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ",
    "SSZ", "K0Z", "K1Z"
  ) %in% names(res$Mexico)))

  # Test methods
  out_print <- capture.output(print(res))
  expect_true(length(out_print) > 0)
  expect_type(out_print, "character")

  out_summary <- capture.output(summary(res))
  expect_true(length(out_summary) > 0)
  expect_type(out_summary, "character")
})

# b) With bias correction
test_that("InputsForOpt returns correct output structure (JPS model with bias correction)", {
  ModelType <- "JPS original"
  Economies <- "Mexico"
  WishBC <- TRUE
  BRWlist <- within(list(
    Cent_Measure = "Mean", gamma = 0.1, N_iter = 10, B = 20, checkBRW = TRUE,
    B_check = 1000, Eigen_rest = 1
  ), N_burn <- round(N_iter * 0.15))

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency,
    WishBRW = TRUE, BRWlist = BRWlist,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ",
    "SSZ", "K0Z", "K1Z"
  ) %in% names(res$Mexico)))
})


# 2) Test error for inconsistent inputs
# a) Example 1
test_that("InputsForOpt throws error for inconsistent input (ModelType mispelled)", {
  ModelType <- "JPS Original" # Typo, capital o in "original"
  Economies <- "Mexico"
  expect_error(InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency,
    CheckInputs = TRUE, verbose = FALSE
  ))
})

# Example 2
test_that("InputsForOpt throws error for inconsistent input (date mispecified)", {
  ModelType <- "JPS original"
  InitialSampleDate <- "15-01-2007" # wrong format
  Economies <- "Mexico"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  expect_error(InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency,
    CheckInputs = TRUE, verbose = FALSE
  ))
})

# Example 3
test_that("InputsForOpt throws error for inconsistent input (frequency mispecified)", {
  ModelType <- "JPS original"
  InitialSampleDate <- "15-01-2007"
  Economies <- "Mexico"
  DataFrequency <- "Bi-annual" # unavailable option
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  expect_error(InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency,
    CheckInputs = TRUE, verbose = FALSE
  ))
})


# 3) Test output structure (Model type "GVAR multi")
# a) version 1
test_that("InputsForOpt returns correct output structure (GVAR single model)", {
  ModelType <- "GVAR single"
  Economies <- c("China", "Mexico")
  GVARlist <- list(
    VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
    t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
  )
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency, GVARlist,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res$China, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ", "SSZ",
    "K0Z", "K1Z"
  ) %in% names(res$China)))
  expect_type(res$China$GVARinputs, "list")
})

# b) version 2
test_that("InputsForOpt returns correct output structure (GVAR multi model)", {
  ModelType <- "GVAR multi"
  Economies <- c("China", "Mexico")
  GVARlist <- list(
    VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
    t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
  )
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency, GVARlist,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ", "SSZ",
    "K0Z", "K1Z"
  ) %in% names(res)))
  expect_type(res$GVARinputs, "list")
})

# c) With bias correction
test_that("InputsForOpt returns correct output structure (GVAR multi model with bias correction)", {
  ModelType <- "GVAR multi"
  Economies <- c("China", "Mexico")
  GVARlist <- list(
    VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
    t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
  )
  WishBC <- TRUE
  BRWlist <- within(list(
    Cent_Measure = "Mean", gamma = 0.1, N_iter = 10, B = 15, checkBRW = TRUE,
    B_check = 100, Eigen_rest = 1
  ), N_burn <- round(N_iter * 0.15))

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency, GVARlist,
    WishBRW = TRUE, BRWlist = BRWlist,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ", "SSZ",
    "K0Z", "K1Z"
  ) %in% names(res)))
  expect_type(res$GVARinputs, "list")
})


# 4) Test output structure (Model type "JLL original")
test_that("InputsForOpt returns correct output structure (JLL model)", {
  ModelType <- "JLL original"
  Economies <- c("China", "Mexico")
  JLLlist <- list(DomUnit = "China")
  GVARlist <- NULL

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

  res <- InputsForOpt(InitialSampleDate, FinalSampleDate, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFrequency, GVARlist, JLLlist,
    CheckInputs = FALSE, verbose = FALSE
  )
  expect_type(res, "list")
  expect_s3_class(res, "ATSMModelInputs")
  expect_true(all(c(
    "mat", "Wpca", "We", "WpcaFull", "Yields", "SpaFact", "RiskFactors", "Gy.0", "K1XQ", "SSZ",
    "K0Z", "K1Z"
  ) %in% names(res)))
  expect_type(res$JLLinputs, "list")
})
