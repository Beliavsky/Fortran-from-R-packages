## -----------------------------------------------------------------------------
library(MultiATSM)

## ----ModFea, message=FALSE, echo=FALSE----------------------------------------
ModelLabels <- c(
  "JPS original", "JPS global", "JPS multi", "GVAR single", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma"
)

# Rows
Tab <- data.frame(matrix(nrow = length(ModelLabels), ncol = 0))
rownames(Tab) <- ModelLabels

# Empty columns
EmptyCol <- c("", "", "", "", "", "", "", "")
Tab$EmptyCol0 <- EmptyCol
# P-dynamics + 2 empty spaces
Tab$PdynIndUnco <- c("x", "", "", "", "", "", "", "")
Tab$PdynIndCo <- c("", "", "", "", "", "", "", "")
Tab$PdynJointUnco <- c("", "x", "x", "", "", "", "", "")
Tab$PdynJointJLL <- c("", "", "", "", "", "x", "x", "x")
Tab$PdynJointGVAR <- c("", "", "", "x", "x", "", "", "")
Tab$EmptyCol1 <- EmptyCol
Tab$EmptyCol2 <- EmptyCol
# Q-dynamics + 2 empty spaces
Tab$QdynInd <- c("x", "x", "", "x", "", "", "", "")
Tab$QdynJoint <- c("", "", "x", "", "x", "x", "x", "x")
Tab$EmptyCol3 <- EmptyCol
Tab$EmptyCol4 <- EmptyCol
# Sigma + 2 empty spaces
Tab$Ponly <- c("", "", "", "", "", "x", "x", "")
Tab$PandQ <- c("x", "x", "x", "x", "x", "", "", "x")
Tab$EmptyCol5 <- EmptyCol
Tab$EmptyCol6 <- EmptyCol
# Dominant Unit
Tab$DomUnit <- c("", "", "", "", "", "x", "", "x")

# Adjust column names
ColNames <- c("", "", "", "", "JLL", "GVAR", "", "", "", "", "", "", "", "", "", "", "")
colnames(Tab) <- ColNames


# Generate the table
suppressWarnings(library(magrittr))


kableExtra::kbl(
  Tab,
  align = "c",
  caption = "Model Features"
) %>%
  kableExtra::kable_classic("striped", full_width = F) %>%
  kableExtra::row_spec(0, font_size = 10) %>%
  kableExtra::add_header_above(c(" " = 2, "Unrestricted" = 1, "Restricted" = 1, "Unrestricted" = 1, "Restricted" = 2, " " = 11)) %>%
  kableExtra::add_header_above(c(" " = 2, "Individual" = 2, "Joint" = 3, " " = 2, "Individual" = 1, "Joint" = 1, " " = 2, "P only" = 1, "P and Q" = 1, " " = 3)) %>%
  kableExtra::add_header_above(c(" " = 2, "P-dynamics" = 5, " " = 2, "Q-dynamics" = 2, " " = 2, "Sigma matrix estimation" = 2, " " = 2, "Dominant Country" = 1), bold = T) %>%
  kableExtra::pack_rows("Unrestricted VAR", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  kableExtra::pack_rows("Restricted VAR (GVAR)", 4, 5, label_row_css = "background-color: #666; color: #fff;") %>%
  kableExtra::pack_rows("Restricted VAR (JLL)", 6, 8, label_row_css = "background-color: #666; color: #fff;")

## -----------------------------------------------------------------------------
LoadData("CM_2024")

## -----------------------------------------------------------------------------
data("Yields")

## -----------------------------------------------------------------------------
data("GlobalMacro")
data("DomMacro")

## -----------------------------------------------------------------------------
data("TradeFlows")

## -----------------------------------------------------------------------------
MacroData <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
YieldData <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))

## -----------------------------------------------------------------------------
Initial_Date <- "2006-09-01" # Format "yyyy-mm-dd"
Final_Date <- "2019-01-01" # Format "yyyy-mm-dd"
DataFrequency <- "Monthly"
GlobalVar <- c("GBC", "VIX") # Global Variables
DomVar <- c("Eco_Act", "Inflation", "Com_Prices", "Exc_Rates") #  Domestic variables
N <- 3 # Number of spanned factors per country
Economies <- c("China", "Mexico", "Uruguay", "Brazil", "Russia")
ModelType <- "JPS original"

## -----------------------------------------------------------------------------
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
RiskFac_TS <- DataForEstimation(
  Initial_Date, Final_Date, Economies, N, FactorLabels, ModelType,
  DataFrequency, MacroData, YieldData
)

## -----------------------------------------------------------------------------
# 1) Model type
ModelType <- "JPS original"

# 2) Risk factor set
Economies <- c("China", "Brazil", "Mexico", "Uruguay")
GlobalVar <- c("GBC", "CPI_OECD") # Global variables
DomVar <- c("Eco_Act", "Inflation") # Domestic variables
N <- 3 # number of spanned factors per country

# 3) Sample span
Initial_Date <- "01-05-2005" # Format: "dd-mm-yyyy"
Final_Date <- "01-12-2019" # Format: "dd-mm-yyyy"

# 4) Frequency of the data
DataFrequency <- "Monthly"

# 5) Risk-neutral stationary constraint
StationarityUnderQ <- FALSE

# 6) Path to save outputs
Folder2Save <- NULL

# 7) Output label
OutputLabel <- "Model_demo"

## -----------------------------------------------------------------------------
VARXtype <- "unconstrained"

## -----------------------------------------------------------------------------
W_type <- "Sample Mean" # Method to compute the transition matrix
t_First_Wgvar <- "2000" # First year of the sample
t_Last_Wgvar <- "2015" # Last year of the sample
DataConnectedness <- TradeFlows # Measure of connectedness across countries

## -----------------------------------------------------------------------------
GVARlist <- list(
  VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2000",
  t_Last_Wgvar = "2015", DataConnectedness = TradeFlows
)

## -----------------------------------------------------------------------------
JLLlist <- list(DomUnit = "China")

## -----------------------------------------------------------------------------
BRWlist <- within(list(
  Cent_Measure = "Mean", gamma = 0.2, N_iter = 500, B = 50,
  checkBRW = TRUE, B_check = 1000, Eigen_rest = 1
), N_burn <- round(N_iter * 0.15))

## ----eval=FALSE---------------------------------------------------------------
# LoadData("CM_2024")
# 
# ModelType <- "JPS original"
# Economies <- "Mexico"
# t0 <- "01-05-2007" # Initial Sample Date (Format: "dd-mm-yyyy")
# tF <- "01-12-2018" # Final Sample Date (Format: "dd-mm-yyyy")
# N <- 3
# GlobalVar <- c("Gl_Eco_Act") # Global Variables
# DomVar <- c("Eco_Act") # Domestic Variables
# FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
# 
# DataFreq <- "Monthly"
# 
# ATSMInputs <- InputsForOpt(t0, tF, ModelType, Yields, GlobalMacro, DomMacro,
#   FactorLabels, Economies, DataFreq,
#   CheckInputs = FALSE
# )

## -----------------------------------------------------------------------------
Horiz <- 100

## -----------------------------------------------------------------------------
DesiredGraphs <- c("Fit", "GIRF", "GFEVD") # Available options are: "Fit", "IRF", "FEVD", "GIRF",
# "GFEVD", "TermPremia".

## -----------------------------------------------------------------------------
WishGraphRiskFac <- FALSE
WishGraphYields <- TRUE
WishOrthoJLLgraphs <- FALSE

## -----------------------------------------------------------------------------
WishFPremia <- TRUE
FPmatLim <- c(60, 120)

## -----------------------------------------------------------------------------
Bootlist <- list(methodBS = "block", BlockLength = 4, ndraws = 50, pctg = 95)

## -----------------------------------------------------------------------------
ForecastList <- list(ForHoriz = 12, t0Sample = 1, t0Forecast = 70, ForType = "Rolling")

## -----------------------------------------------------------------------------
w <- pca_weights_one_country(Yields, Economy = "Uruguay")

## ----fig.cap = "Yield loadings on the spanned factors", echo=FALSE------------
LabSpaFac <- c("Level", "Slope", "Curvature")
N <- length(LabSpaFac)

w_pca <- data.frame(t(w[1:N, ]))
colnames(w_pca) <- LabSpaFac
w_pca$mat <- c(0.25, 0.5, 1, 3, 5, 10) # vector of maturitie

# Prepare plots
colors <- c("Level" = "blue", "Slope" = "green", "Curvature" = "red")

g <- ggplot2::ggplot(data = w_pca, ggplot2::aes(x = mat)) +
  ggplot2::geom_line(ggplot2::aes(y = Level, color = "Level"), linewidth = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = Slope, color = "Slope"), linewidth = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = Curvature, color = "Curvature"), linewidth = 0.7) +
  ggplot2::labs(color = "Legend") +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position = "top", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 8)) +
  ggplot2::xlab("Maturity (Years)") +
  ggplot2::geom_hline(yintercept = 0)

print(g)

## -----------------------------------------------------------------------------
data("Yields")
Economies <- c("China", "Brazil", "Mexico", "Uruguay")
N <- 3
SpaFact <- Spanned_Factors(Yields, Economies, N)

## -----------------------------------------------------------------------------
data(RiskFacFull)
PdynPara <- VAR(RiskFacFull, VARtype = "unconstrained")

## -----------------------------------------------------------------------------
FactorsChina <- RiskFacFull[1:7, ]
PdynPara <- VAR(FactorsChina, VARtype = "unconstrained")

## -----------------------------------------------------------------------------
data(GVARFactors)
GVARinputs <- list(Economies = Economies, GVARFactors = GVARFactors, VARXtype = "constrained: Inflation")

## -----------------------------------------------------------------------------
data("TradeFlows")
t_First <- "2006"
t_Last <- "2019"
Economies <- c("China", "Brazil", "Mexico", "Uruguay")
type <- "Sample Mean"
W_gvar <- Transition_Matrix(t_First, t_Last, Economies, type, TradeFlows)
print(W_gvar)

## -----------------------------------------------------------------------------
data("GVARFactors")
GVARinputs <- list(
  Economies = Economies, GVARFactors = GVARFactors, VARXtype = "unconstrained",
  Wgvar = W_gvar
)
N <- 3

GVARpara <- GVAR(GVARinputs, N, CheckInputs = TRUE)

## -----------------------------------------------------------------------------
ModelType <- "JLL original"
JLLinputs <- list(
  Economies = Economies, DomUnit = "China", WishSigmas = TRUE, SigmaNonOrtho = NULL,
  JLLModelType = ModelType
)

## ----eval=FALSE---------------------------------------------------------------
# data("RiskFacFull")
# N <- 3
# JLLpara <- JLL(RiskFacFull, N, JLLinputs, CheckInputs = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# ########################################################################################################
# #################################### USER INPUTS #######################################################
# ########################################################################################################
# library(MultiATSM)
# # A) Load database data
# LoadData("CM_2024")
# 
# # B) GENERAL model inputs
# ModelType <- "JPS original" # available options: "JPS original", "JPS global", "GVAR single", "JPS multi",
# # "GVAR multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma".
# 
# Economies <- c("China", "Brazil") # Names of the economies from the economic system
# GlobalVar <- c("Gl_Eco_Act") # Global Variables
# DomVar <- c("Eco_Act") # Country-specific variables
# N <- 2 # Number of spanned factors per country
# 
# t0_sample <- "01-12-2007" # Format: "dd-mm-yyyy"
# tF_sample <- "01-12-2019" # Format: "dd-mm-yyyy"
# 
# OutputLabel <- "Test" # label of the model for saving the file
# DataFreq <- "Monthly" # Frequency of the data
# 
# Folder2Save <- NULL # Folder to save the required graphical outputs. NULL -> temporary directory
# StatQ <- FALSE # Wish to impose stationary condition for the eigenvalues of each country
# 
# # B.1) SPECIFIC model inputs
# #################################### GVAR-based models ##################################################
# GVARlist <- list(
#   VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
#   t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
# )
# # VARXtype: Available options "unconstrained" or "constrained" (VARX)
# # W_type: Method to compute the transition matrix. Options:"Time-varying" or "Sample Mean"
# # t_First_Wgvar: First year of the sample (transition matrix)
# # t_Last_Wgvar:  Last year of the sample (transition matrix)
# # DataConnectedness: measure of connectedness across countries
# #################################### JLL-based models ###################################################
# JLLlist <- list(DomUnit = "China")
# # DomUnit: name of the economy of the economic system, or "None" for the model "JLL No DomUnit"
# ###################################### BRW inputs  ######################################################
# WishBC <- FALSE # Wish to estimate the model with the bias correction method of BRW (2012)
# BRWlist <- within(list(
#   Cent_Measure = "Mean", gamma = 0.05, N_iter = 250, B = 50, checkBRW = TRUE,
#   B_check = 1000, Eigen_rest = 1
# ), N_burn <- round(N_iter * 0.15))
# # Cent_Measure: "Mean" = compute the mean; "Median" = compute the median
# # gamma: Adjustment parameter
# # N_iter:  Number of iteration to be conserved
# # N_burn:  Number of iteration to be discarded
# # B: Number of bootstrap samples
# # checkBRW: wishes to perform closeness check: TRUE or FALSE
# # B_check: If checkBRW is chosen as TRUE, then choose number of bootstrap samples used in the check
# 
# # C) Decide on Settings for numerical outputs
# WishFPremia <- TRUE # Wish to estimate the forward premia
# FPmatLim <- c(60, 120) #  If the forward premia is desired, then choose the Min and max maturities of the
# # forward premia. Otherwise set NA
# Horiz <- 30
# DesiredGraphs <- c("Fit", "IRF", "TermPremia") # "Fit", "IRF", "FEVD", "GIRF", "GFEVD", "TermPremia"
# WishGraphRiskFac <- FALSE
# WishGraphYields <- TRUE
# WishOrthoJLLgraphs <- FALSE
# 
# # D) Bootstrap settings
# WishBootstrap <- TRUE
# BootList <- list(methodBS = "bs", BlockLength = 4, ndraws = 4, pctg = 95)
# # methodBS: bootstrap method. Available options: (i) 'bs' ; (ii) 'wild'; (iii) 'block'
# # BlockLength: Block length, necessary input for the block bootstrap method
# # ndraws: number of bootstrap draws
# # pctg: confidence level
# 
# # E) Out-of-sample forecast
# WishForecast <- TRUE
# ForecastList <- list(ForHoriz = 12, t0Sample = 1, t0Forecast = 131, ForType = "Rolling")
# # ForHoriz: forecast horizon
# # t0Sample:   initial sample date
# # t0Forecast:  last sample date for the first forecast
# # ForType: Available options are "Rolling" or "Expanding"
# 
# #########################################################################################################
# ## NO NEED TO MAKE CHANGES FROM HERE:
# ## The sections below automatically process the inputs provided above, run the model
# ## estimation, generate the numerical and graphical outputs, and save results.
# #########################################################################################################
# 
# # 2) Minor preliminary work: get the sets of factor labels and  a vector of common maturities
# FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
# 
# # 3) Prepare the inputs of the likelihood function
# ATSMInputs <- InputsForOpt(
#   t0_sample, tF_sample, ModelType, Yields, GlobalMacro, DomMacro,
#   FactorLabels, Economies, DataFreq, GVARlist, JLLlist, WishBC, BRWlist
# )
# 
# # 4) Optimization of the ATSM (Point Estimates)
# ModelParaList <- Optimization(ATSMInputs, StatQ, DataFreq, FactorLabels, Economies, ModelType)
# 
# # 5) Numerical and graphical outputs
# # a) Prepare list of inputs for graphs and numerical outputs
# InputsForOutputs <- InputsForOutputs(
#   ModelType, Horiz, DesiredGraphs, OutputLabel, StatQ, DataFreq,
#   WishGraphYields, WishGraphRiskFac, WishOrthoJLLgraphs, WishFPremia,
#   FPmatLim, WishBootstrap, BootList, WishForecast, ForecastList
# )
# 
# # b) Fit, IRF, FEVD, GIRF, GFEVD, and Term Premia
# NumericalOutputs <- NumOutputs(
#   ModelType, ModelParaList, InputsForOutputs, FactorLabels, Economies,
#   Folder2Save
# )
# 
# # c) Confidence intervals (bootstrap analysis)
# BootstrapAnalysis <- Bootstrap(
#   ModelType, ModelParaList, NumericalOutputs, Economies, InputsForOutputs,
#   FactorLabels, JLLlist, GVARlist, WishBC, BRWlist, Folder2Save
# )
# 
# # 6) Out-of-sample forecasting
# Forecasts <- ForecastYields(
#   ModelType, ModelParaList, InputsForOutputs, FactorLabels, Economies,
#   JLLlist, GVARlist, WishBC, BRWlist, Folder2Save
# )

## -----------------------------------------------------------------------------
data("Out_Example")
print(Out_Example$ATSMInputs)

## -----------------------------------------------------------------------------
summary(Out_Example$ATSMInputs)

## -----------------------------------------------------------------------------
summary(Out_Example$ModelPara)

## ----fig.alt="GIRFs - Responses of yields curves in China and Brazil to a global economic shocks"----
GIRF_all <- autoplot(Out_Example$NumOut, type = "GIRF_Yields")
GIRF_all$Brazil$Gl_Eco_Act
GIRF_all$China$Gl_Eco_Act

## ----fig.alt="GIRFs - Responses of 12-month maturity bonds in China and Brazil to several shocks"----
GIRF_all_Boot <- autoplot(Out_Example$Bootstrap, Out_Example$NumOut, type = "GIRF_Yields_Boot")
GIRF_all_Boot$China$Y12M_China
GIRF_all_Boot$Brazil$Y12M_Brazil

## ----fig.alt="Root mean square error (RMSE)"----------------------------------
plot(Out_Example$Forecasts)

