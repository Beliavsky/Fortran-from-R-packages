## code to prepare pre-loaded outputs

########################################################################################################
#################################### USER INPUTS #######################################################
########################################################################################################
rm(list = ls())
cat("\014")
library(MultiATSM)
# A) Load database data
LoadData("CM_2024")

# B) GENERAL model inputs
ModelType <- "JPS original"

Economies <- c("Brazil")
GlobalVar <- c("Gl_Eco_Act")
DomVar <- c("Eco_Act")
N <- 1

t0_sample <- "01-07-2006"
tF_sample <- "01-07-2019"

OutputLabel <- "Test"
DataFreq <- "Monthly"

Folder2Save <- NULL
StatQ <- FALSE

# Model-specific list
GVARlist <- list(
  VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
  t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
)

JLLlist <- list(DomUnit = "China")

WishBC <- FALSE
BRWlist <- within(list(
  Cent_Measure = "Mean", gamma = 0.05, N_iter = 250, B = 50, checkBRW = TRUE,
  B_check = 1000, Eigen_rest = 1
), N_burn <- round(N_iter * 0.15))


# C) Decide on Settings for numerical outputs
WishFPremia <- TRUE
FPmatLim <- c(60, 120)
Horiz <- 30
DesiredGraphs <- c()
WishGraphRiskFac <- FALSE
WishGraphYields <- TRUE
WishOrthoJLLgraphs <- FALSE

# D) Bootstrap settings
WishBootstrap <- TRUE
BootList <- list(methodBS = "bs", BlockLength = 4, ndraws = 5, pctg = 95)

# E) Out-of-sample forecast
WishForecast <- FALSE
ForecastList <- list(ForHoriz = 12, t0Sample = 1, t0Forecast = 131, ForType = "Rolling")

# 2) Minor preliminary work: get the sets of factor labels and  a vector of common maturities
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

# 3) Prepare the inputs of the likelihood function
ATSMInputs <- InputsForOpt(
  t0_sample, tF_sample, ModelType, Yields, GlobalMacro, DomMacro,
  FactorLabels, Economies, DataFreq, GVARlist, JLLlist, WishBC, BRWlist
)

# 4) Optimization of the ATSM (Point Estimates)
ParaSetEx <- Optimization(ATSMInputs, StatQ, DataFreq, FactorLabels, Economies, ModelType)

# 5) Numerical and graphical outputs
# a) Prepare list of inputs for graphs and numerical outputs
InpForOutEx <- InputsForOutputs(
  ModelType, Horiz, DesiredGraphs, OutputLabel, StatQ, DataFreq,
  WishGraphYields, WishGraphRiskFac, WishOrthoJLLgraphs, WishFPremia,
  FPmatLim, WishBootstrap, BootList, WishForecast, ForecastList
)

# b) Fit, IRF, FEVD, GIRF, GFEVD, and Term Premia
NumOutEx <- NumOutputs(
  ModelType, ParaSetEx, InpForOutEx, FactorLabels, Economies,
  Folder2Save
)

# Output to export
usethis::use_data(ParaSetEx, overwrite = TRUE)
usethis::use_data(InpForOutEx, overwrite = TRUE)
usethis::use_data(NumOutEx, overwrite = TRUE)
