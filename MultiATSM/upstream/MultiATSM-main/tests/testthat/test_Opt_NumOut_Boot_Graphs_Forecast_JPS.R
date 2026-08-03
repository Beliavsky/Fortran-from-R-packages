library(testthat)
library(MultiATSM)


# 0) Load inputs
extract_ggplots <- function(x) {
  if (inherits(x, "ggplot")) {
    list(x)
  } else if (is.list(x)) {
    unlist(lapply(x, extract_ggplots), recursive = FALSE)
  } else {
    list()
  }
}

# Inputs required for the optimization
LoadData("CM_2024")
Economy <- "Brazil"
ModelType <- "JPS original"
t0 <- "01-05-2007"
tF <- "01-12-2018"
N <- 1
GlobalVar <- "Gl_Eco_Act"
DomVar <- "Eco_Act"
DataFreq <- "Monthly"
StatQ <- TRUE
WishFP <- FALSE
FPLim <- c(24, 36)
WishBC <- FALSE

# Inputs required for the computation of the numerical outputs
Horiz <- 4
OutputLabel <- "Model_demo"
DesiredGraphs <- c()
WGYields <- FALSE
WGFac <- FALSE
WGJLL <- FALSE

# Bootstrap settings
WishBoot <- TRUE
BootList <- list(methodBS = "bs", BlockLength = 4, ndraws = 2, pctg = 95)

# Forecasting setting
WishFor <- TRUE
ForList <- list(ForHoriz = 6, t0Sample = 1, t0Forecast = 133, ForType = "Expanding")


# 1) Set of tests
test_that("Optimization + Outputs + Graphs return correct structure (JPS model)", {
  # Run optimization
  FacLab <- LabFac(N, DomVar, GlobalVar, Economy, ModelType)
  ATSMInputs <- InputsForOpt(
    t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
    FacLab, Economy, DataFreq,
    CheckInputs = FALSE, verbose = FALSE
  )

  res_Opt <- Optimization(ATSMInputs, StatQ, DataFreq, FacLab, Economy, ModelType, verbose = FALSE)

  # --- A) Optimization structure ---
  expect_type(res_Opt, "list")
  expect_s3_class(res_Opt, "ATSMModelOutputs")
  out_summary <- capture.output(summary(res_Opt))
  expect_true(length(out_summary) > 0)
  expect_type(out_summary, "character")

  expect_true(all(c("Inputs", "ModEst") %in% names(res_Opt[[ModelType]][[Economy]])))
  expect_true(all(c("Y", "AllFactors", "mat", "N", "dt", "Wpca") %in%
    names(res_Opt[[ModelType]][[Economy]]$Inputs)))
  expect_true(all(c("Max_llk", "Q", "P") %in% names(res_Opt[[ModelType]][[Economy]]$ModEst)))
  expect_true(all(c("K1XQ", "r0", "se", "VarYields") %in% names(res_Opt[[ModelType]][[Economy]]$ModEst$Q)))
  expect_true(all(c("SSZ", "K0Z", "K1Z", "Gy.0") %in% names(res_Opt[[ModelType]][[Economy]]$ModEst$P)))
  expect_type(res_Opt[[ModelType]][[Economy]]$ModEst$Max_llk, "double")

  # --- B) Numerical outputs ---
  InputsForOutputs <- InputsForOutputs(
    ModelType, Horiz, DesiredGraphs, OutputLabel, StatQ, DataFreq, WGYields,
    WGFac, WGJLL, WishFP, FPLim, WishBoot, BootList, WishFor, ForList
  )
  res_NumOut <- NumOutputs(ModelType, res_Opt, InputsForOutputs, FacLab, Economy, verbose = FALSE)

  expect_type(res_NumOut, "list")
  expect_s3_class(res_NumOut, "ATSMNumOutputs")

  expect_true(all(c("PC var explained", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremiaDecomp") %in% names(res_NumOut)))
  expect_equal(length(res_NumOut$`PC var explained`[[Economy]]), N)

  # compact checks of nested names
  checks <- list(
    "res_NumOut$Fit[[ModelType]][[Economy]]"         = c("Yield Fit", "Yield Model Implied"),
    "res_NumOut$IRF[[ModelType]][[Economy]]"         = c("Factors", "Yields"),
    "res_NumOut$FEVD[[ModelType]][[Economy]]"        = c("Factors", "Yields"),
    "res_NumOut$GIRF[[ModelType]][[Economy]]"        = c("Factors", "Yields"),
    "res_NumOut$GFEVD[[ModelType]][[Economy]]"       = c("Factors", "Yields"),
    "res_NumOut$TermPremiaDecomp$RiskPremia"         = c("Term Premia", "Expected Component")
  )

  for (nm in names(checks)) {
    expect_true(all(checks[[nm]] %in% names(eval(parse(text = nm)))))
  }

  # --- C) Graphs ---
  plot_types <- c(
    "RiskFactors", "Fit",
    #"IRF_Factors", "IRF_Yields", # To save time (will be tested in the bootstrap setting)
    #"GIRF_Factors", "GIRF_Yields", "FEVD_Factors", # To save time (will be tested in the bootstrap setting)
    #"FEVD_Yields", "GFEVD_Factors", "GFEVD_Yields", # To save time (will be tested in the bootstrap setting)
    "TermPremia"
  )

  plot_list <- lapply(plot_types, function(tp) autoplot(res_NumOut, type = tp))

  all_plots <- extract_ggplots(plot_list)

  for (p in all_plots) {
    expect_s3_class(p, "ggplot")
  }

  # --- D) Bootstrap analysis ---
  res_Boot <- Bootstrap(ModelType, res_Opt, res_NumOut, Economy, InputsForOutputs, FacLab,
    JLLlist = NULL, GVARlist = NULL, WishBC, BRWlist = NULL, verbose = FALSE
  )
  expect_type(res_Boot, "list")
  expect_s3_class(res_Boot, "ATSMModelBoot")

  plot_types_Boot <- c(
    #"IRF_Factors_Boot",
    "IRF_Yields_Boot",
    #"GIRF_Factors_Boot",
    "GIRF_Yields_Boot",
    #"FEVD_Factors_Boot",
    "FEVD_Yields_Boot",
    #"GFEVD_Factors_Boot",
    "GFEVD_Yields_Boot"
  )

  plot_list_Boot <- lapply(plot_types_Boot, function(tp) autoplot(res_Boot, res_NumOut, type = tp))

  all_plots <- extract_ggplots(plot_list_Boot)

  for (p in all_plots) {
    expect_s3_class(p, "ggplot")
  }

  # --- E) Out-of-sample forecasting analysis ---
  res_For <- ForecastYields(ModelType, res_Opt, InputsForOutputs, FacLab, Economy, WishBRW = WishBC, verbose = FALSE)
  expect_type(res_For, "list")
  expect_s3_class(res_For, "ATSMModelForecast")

  out_plot <- capture.output(print(res_For))
  expect_true(length(out_plot) > 0)
  expect_type(out_plot, "character")
})
