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
Economies <- c("China", "Mexico")
ModelType <- "JLL original"
t0 <- "01-05-2007"
tF <- "01-12-2018"
N <- 1
GlobalVar <- "Gl_Eco_Act"
DomVar <- "Eco_Act"
DataFreq <- "Monthly"
StatQ <- FALSE
WishFP <- TRUE
FPLim <- c(24, 36)

WishBC <- FALSE
JLLlist <- list(DomUnit = "China")
GVARlist <- NULL

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


test_that("Optimization + Outputs + Graphs return correct structure (JLL model)", {
  # Run optimization
  FacLab <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  ATSMInputs <- InputsForOpt(
    t0, tF, ModelType, Yields, GlobalMacro, DomMacro, FacLab, Economies,
    DataFreq, GVARlist, JLLlist,
    CheckInputs = FALSE, verbose = FALSE
  )

  res_Opt <- Optimization(ATSMInputs, StatQ, DataFreq, FacLab, Economies, ModelType, verbose = FALSE)

  # --- A) Optimization structure ---
  expect_type(res_Opt, "list")
  expect_s3_class(res_Opt, "ATSMModelOutputs")
  out_summary <- capture.output(summary(res_Opt))
  expect_true(length(out_summary) > 0)
  expect_type(out_summary, "character")

  expect_true(all(c("Inputs", "ModEst") %in% names(res_Opt[[ModelType]])))
  expect_true(all(c("Y", "AllFactors", "mat", "N", "dt", "Wpca") %in%
    names(res_Opt[[ModelType]]$Inputs)))
  expect_true(all(c("Max_llk", "Q", "P") %in% names(res_Opt[[ModelType]]$ModEst)))
  expect_true(all(c("K1XQ", "r0", "se", "VarYields") %in% names(res_Opt[[ModelType]]$ModEst$Q)))
  expect_true(all(c("SSZ", "K0Z", "K1Z", "Gy.0") %in% names(res_Opt[[ModelType]]$ModEst$P)))
  expect_type(res_Opt[[ModelType]]$ModEst$Max_llk, "double")

  # --- B) Numerical outputs ---
  InputsForOutputs <- InputsForOutputs(
    ModelType, Horiz, DesiredGraphs, OutputLabel, StatQ, DataFreq, WGYields,
    WGFac, WGJLL, WishFP, FPLim, WishBoot, BootList, WishFor, ForList
  )
  res_NumOut <- NumOutputs(ModelType, res_Opt, InputsForOutputs, FacLab, Economies, verbose = FALSE)

  expect_type(res_NumOut, "list")
  expect_s3_class(res_NumOut, "ATSMNumOutputs")

  expect_true(all(c("PC var explained", "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremiaDecomp") %in% names(res_NumOut)))
  expect_equal(length(res_NumOut$`PC var explained`), length(Economies)*N)

  # compact checks of nested names
  checks <- list(
    # IRF
    "res_NumOut$IRF[[ModelType]]" = c("Factors", "Yields"),
    "res_NumOut$IRF[[ModelType]]$Factors" = c("NonOrtho", "Ortho"),
    "res_NumOut$IRF[[ModelType]]$Yields" = c("NonOrtho", "Ortho"),
    # FEVD
    "res_NumOut$FEVD[[ModelType]]" = c("Factors", "Yields"),
    "res_NumOut$FEVD[[ModelType]]$Factors" = c("NonOrtho", "Ortho"),
    "res_NumOut$FEVD[[ModelType]]$Yields" = c("NonOrtho", "Ortho"),
    # GIRF
    "res_NumOut$GIRF[[ModelType]]" = c("Factors", "Yields"),
    "res_NumOut$GIRF[[ModelType]]$Factors" = c("NonOrtho", "Ortho"),
    "res_NumOut$GIRF[[ModelType]]$Yields" = c("NonOrtho", "Ortho"),
    # GFEVD
    "res_NumOut$GFEVD[[ModelType]]" = c("Factors", "Yields"),
    "res_NumOut$GFEVD[[ModelType]]$Factors" = c("NonOrtho", "Ortho"),
    "res_NumOut$GFEVD[[ModelType]]$Yields" = c("NonOrtho", "Ortho"),
    # TP
    "res_NumOut$TermPremiaDecomp$RiskPremia" = c("Term Premia", "Expected Component")
  )

  for (nm in names(checks)) {
    expect_true(all(checks[[nm]] %in% names(eval(parse(text = nm)))))
  }

  # Handle "Fit" separately
  check_Fit <- list(Fit = c("Yield Fit", "Yield Model Implied"))
  for (eco in Economies) {
    expect_true(all(check_Fit$check_Fit %in% names(res_NumOut$Fit[[ModelType]][[eco]])))
  }

  # --- C) Graphs ---
  plot_types <- c(
    "RiskFactors", "Fit",
    #"IRF_Factors", "IRF_Yields", "IRF_Factors_Ortho", # To save time
    "IRF_Yields_Ortho",
    #"GIRF_Factors", "GIRF_Yields", "GIRF_Factors_Ortho", # To save time
     "GIRF_Yields_Ortho",
    #"FEVD_Factors", "FEVD_Yields", "FEVD_Factors_Ortho", # To save time
     "FEVD_Yields_Ortho",
    #"GFEVD_Factors", "GFEVD_Yields", "GFEVD_Factors_Ortho", # To save time
     "GFEVD_Yields_Ortho",
    "TermPremia"
  )

  plot_list <- lapply(plot_types, function(tp) autoplot(res_NumOut, type = tp))

  all_plots <- extract_ggplots(plot_list)

  for (p in all_plots) {
    expect_s3_class(p, "ggplot")
  }

  # --- D) Bootstrap analysis ---
  res_Boot <- Bootstrap(ModelType, res_Opt, res_NumOut, Economies, InputsForOutputs, FacLab,
    JLLlist,
    GVARlist = NULL, WishBC = FALSE, BRWlist = NULL, verbose = FALSE
  )
  expect_type(res_Boot, "list")
  expect_s3_class(res_Boot, "ATSMModelBoot")

  plot_types_Boot <- c(
    #"IRF_Factors_Boot", "IRF_Yields_Boot", "IRF_Factors_Ortho_Boot",   # to save time
     "IRF_Yields_Ortho_Boot",
    #"GIRF_Factors_Boot", "GIRF_Yields_Boot", "GIRF_Factors_Ortho_Boot", # to save time
     "GIRF_Yields_Ortho_Boot",
    #"FEVD_Factors_Boot", "FEVD_Yields_Boot", "FEVD_Factors_Ortho_Boot", # to save time
     "FEVD_Yields_Ortho_Boot",
    #"GFEVD_Factors_Boot", "GFEVD_Yields_Boot", "GFEVD_Factors_Ortho_Boot", # to save time
     "GFEVD_Yields_Ortho_Boot"
  )

  plot_list_Boot <- lapply(plot_types_Boot, function(tp) autoplot(res_Boot, res_NumOut, type = tp))

  all_plots <- extract_ggplots(plot_list_Boot)

  for (p in all_plots) {
    expect_s3_class(p, "ggplot")
  }


  # --- E) Out-of-sample forecasting analysis ---
  res_For <- ForecastYields(ModelType, res_Opt, InputsForOutputs, FacLab, Economies, JLLlist,
    GVARlist = NULL,
    WishBRW = WishBC, verbose = FALSE
  )
  expect_type(res_For, "list")
  expect_s3_class(res_For, "ATSMModelForecast")

  out_plot <- capture.output(print(res_For))
  expect_true(length(out_plot) > 0)
  expect_type(out_plot, "character")
})
