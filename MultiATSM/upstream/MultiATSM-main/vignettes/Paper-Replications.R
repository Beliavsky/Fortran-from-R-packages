## -----------------------------------------------------------------------------
library(MultiATSM)

## -----------------------------------------------------------------------------
# A) Load database data
LoadData("BR_2017")

# B) GENERAL model inputs
ModelType <- "JPS original"

Economies <- c("US") # Names of the economies from the economic system
GlobalVar <- c() # Global Variables
DomVar <- c("GRO", "INF") # Country-specific variables
N <- 3 # Number of spanned factors per country

t0_sample <- "January-1985"
tF_sample <- "December-2007"

DataFreq <- "Monthly" # Frequency of the data

StatQ <- FALSE # Stationary condition

#########################################################################################################
############################### NO NEED TO MAKE CHANGES FROM HERE #######################################
#########################################################################################################
# 2) Minor preliminary work
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

Yields <- t(BR_jps_out$Y)
DomesticMacroVar <- t(BR_jps_out$M.o)
GlobalMacroVar <- c()

# 3) Prepare the inputs of the likelihood function
ATSMInputs <- InputsForOpt(t0_sample, tF_sample, ModelType, Yields, GlobalMacroVar, DomesticMacroVar,
  FactorLabels, Economies, DataFreq,
  verbose = FALSE
)

# 4) Optimization of the model
ModelPara <- Optimization(ATSMInputs, StatQ, DataFreq, FactorLabels, Economies, ModelType,
  verbose = FALSE
)

## ----echo= FALSE--------------------------------------------------------------
options(scipen = 100) # eliminate the scientific notation

RowsQ <- c("$r0$", "$\\lambda_1$", "$\\lambda_2$", "$\\lambda_3$")
TableQ <- data.frame(matrix(NA, ncol = 0, nrow = length(RowsQ)))
row.names(TableQ) <- RowsQ

PackageQ <- c(ModelPara$`JPS original`$US$ModEst$Q$r0, diag(ModelPara$`JPS original`$US$ModEst$Q$K1XQ))
BRq <- c(BR_jps_out$est.llk$rho0.cP, diag(BR_jps_out$est.llk$KQ.XX))
TableQ$MultiATSM <- PackageQ
TableQ$"BR (2017)" <- BRq

TableQ <- round(TableQ, digits = 5)

suppressWarnings(library(magrittr))

kableExtra::kbl(TableQ, align = "c", caption = "$Q$-dynamics parameters") %>%
  kableExtra::kable_classic("striped", full_width = F) %>%
  kableExtra::row_spec(0, font_size = 14) %>%
  kableExtra::footnote(general = " $\\lambda$'s are the eigenvalues from the risk-neutral feedback matrix and $r0$ is the long-run mean of the short rate under Q.")

## ----PdynTab, echo= FALSE-----------------------------------------------------
RowsP <- c("PC1", "PC2", "PC3", "GRO", "INF")
ColP <- c(" ", RowsP)

# 1) K0Z and K1Z
# a) Bauer and Rudebusch coefficients
TablePbr <- data.frame(matrix(NA, ncol = length(ColP), nrow = length(RowsP)))
row.names(TablePbr) <- RowsP
colnames(TablePbr) <- ColP

TablePbr[[ColP[1]]] <- BR_jps_out$est.llk$KP.0Z
for (j in 1:length(RowsP)) {
  TablePbr[[RowsP[j]]] <- BR_jps_out$est.llk$KP.ZZ[, j]
}

TablePbr <- round(TablePbr, digits = 5)

# b) MultiATSM coefficients
TablePMultiATSM <- data.frame(matrix(NA, ncol = length(ColP), nrow = length(RowsP)))
row.names(TablePMultiATSM) <- RowsP
colnames(TablePMultiATSM) <- ColP

# Estimate the P-dynamics using the weights provided by BR necessary for replication)
PP <- BR_jps_out$W[1:N, ] %*% Yields
ZZ <- rbind(PP, DomesticMacroVar)
Pdyncoef <- VAR(ZZ, "unconstrained")

TablePMultiATSM[[ColP[1]]] <- Pdyncoef$K0Z
for (j in 1:length(RowsP)) {
  TablePMultiATSM[[RowsP[j]]] <- Pdyncoef$K1Z[, j]
}


TablePMultiATSM <- round(TablePMultiATSM, digits = 5)
TableP <- rbind(TablePbr, TablePMultiATSM)
row.names(TableP) <- c(RowsP, paste(RowsP, " ", sep = "")) # trick to avoid rows to have the same name

kableExtra::kbl(TableP, align = "c", caption = "$P$-dynamics parameters") %>%
  kableExtra::kable_classic("striped", full_width = F) %>%
  kableExtra::row_spec(0, font_size = 14) %>%
  kableExtra::add_header_above(c(" " = 1, "K0Z" = 1, "K1Z" = 5), bold = T) %>%
  kableExtra::pack_rows("BR (2017)", 1, 5) %>%
  kableExtra::pack_rows("MultiATSM", 6, 10) %>%
  kableExtra::footnote(general = " $K0Z$ is the intercept and $K1Z$ is feedback matrix from the $P$-dynamics.")

## ----echo= FALSE--------------------------------------------------------------
se <- data.frame(BR_jps_out$est.llk$sigma.e, ModelPara$`JPS original`$US$ModEst$Q$se)
rownames(se) <- "se"
colnames(se) <- c("MultiATSM", "BR (2017)")
se <- round(se, digits = 7)

kableExtra::kbl(se, align = "c", caption = "Portfolio of yields with errors") %>%
  kableExtra::kable_classic("striped", full_width = F) %>%
  kableExtra::row_spec(0, font_size = 14) %>%
  kableExtra::footnote(general = " $se$ is the standard deviation of the portfolio of yields observed with errors.")

## ----eval=FALSE---------------------------------------------------------------
# # A) Load database data
# LoadData("CM_2024")
# 
# # B) GENERAL model inputs
# ModelType <- "GVAR multi" # Options: "GVAR multi" or "JLL original".
# 
# Economies <- c("China", "Brazil", "Mexico", "Uruguay")
# GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation")
# DomVar <- c("Eco_Act", "Inflation")
# N <- 3
# 
# t0_sample <- "01-06-2004"
# tF_sample <- "01-01-2020"
# 
# OutputLabel <- "CM_jfec"
# DataFreq <- "Monthly"
# Folder2Save <- NULL
# 
# StatQ <- FALSE
# 
# # B.1) SPECIFIC model inputs
# #################################### GVAR-based models ##################################################
# GVARlist <- list(
#   VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2004",
#   t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
# )
# #################################### JLL-based models ###################################################
# JLLlist <- list(DomUnit = "China")
# ###################################### BRW inputs  ######################################################
# WishBC <- TRUE
# BRWlist <- within(list(
#   Cent_Measure = "Mean", gamma = 0.001, N_iter = 200, B = 50, checkBRW = TRUE,
#   B_check = 1000, Eigen_rest = 1
# ), N_burn <- round(N_iter * 0.15))
# 
# # C) Decide on Settings for numerical outputs
# WishFPremia <- TRUE
# FPmatLim <- c(24, 36)
# 
# Horiz <- 25
# DesiredGraphs <- c("GIRF", "GFEVD", "TermPremia")
# WishGraphRiskFac <- FALSE
# WishGraphYields <- TRUE
# WishOrthoJLLgraphs <- TRUE
# 
# # D) Bootstrap settings
# WishBootstrap <- FALSE #  Set it to TRUE, if bootstrap is desired
# BootList <- list(methodBS = "bs", BlockLength = 4, ndraws = 1000, pctg = 95)
# 
# # E) Out-of-sample forecast
# WishForecast <- TRUE
# ForecastList <- list(ForHoriz = 12, t0Sample = 1, t0Forecast = 100, ForType = "Rolling")
# 
# #########################################################################################################
# ############################### NO NEED TO MAKE CHANGES FROM HERE #######################################
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
#   ModelType, ModelParaList, InputsForOutputs, FactorLabels,
#   Economies, Folder2Save
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

## ----eval=FALSE---------------------------------------------------------------
# # A) Load database data
# LoadData("CM_2023")
# 
# # B) GENERAL model inputs
# ModelType <- "GVAR multi"
# 
# Economies <- c("Brazil", "India", "Russia", "Mexico")
# GlobalVar <- c("US_Output_growth", "China_Output_growth", "SP500")
# DomVar <- c("Inflation", "Output_growth", "CDS", "COVID")
# N <- 2
# 
# t0_sample <- "22-03-2020"
# tF_sample <- "26-09-2021"
# 
# OutputLabel <- "CM_EM"
# DataFreq <- "Weekly"
# F2S <- NULL
# 
# StatQ <- FALSE
# 
# # B.1) SPECIFIC model inputs
# #################################### GVAR-based models ##################################################
# GVARlist <- list(
#   VARXtype = "constrained: COVID", W_type = "Sample Mean", t_First_Wgvar = "2015",
#   t_Last_Wgvar = "2020", DataConnectedness = TradeFlows_covid
# )
# ###################################### BRW inputs  ######################################################
# WishBC <- FALSE
# 
# # C) Decide on Settings for numerical outputs
# WishFPremia <- TRUE
# FPmatLim <- c(47, 48)
# 
# Horiz <- 12
# DesiredGraphs <- c("GIRF", "GFEVD", "TermPremia")
# WishGraphRiskFac <- FALSE
# WishGraphYields <- TRUE
# WishOrthoJLLgraphs <- FALSE
# 
# # D) Bootstrap settings
# WishBootstrap <- TRUE #  YES: 1; No = 0.
# BootList <- list(methodBS = "bs", BlockLength = 4, ndraws = 100, pctg = 95)
# 
# #########################################################################################################
# ############################### NO NEED TO MAKE CHANGES FROM HERE #######################################
# #########################################################################################################
# 
# # 2) Minor preliminary work: get the sets of factor labels and  a vector of common maturities
# FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
# 
# # 3) Prepare the inputs of the likelihood function
# ATSMInputs <- InputsForOpt(
#   t0_sample, tF_sample, ModelType, Yields_covid, GlobalMacro_covid, DomMacro_covid, FactorLabels,
#   Economies, DataFreq, GVARlist
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
#   FPmatLim, WishBootstrap, BootList
# )
# 
# # b) Fit, IRF, FEVD, GIRF, GFEVD, and Term Premia
# NumericalOutputs <- NumOutputs(
#   ModelType, ModelParaList, InputsForOutputs, FactorLabels,
#   Economies, F2S
# )
# 
# # c) Confidence intervals (bootstrap analysis)
# BootstrapAnalysis <- Bootstrap(ModelType, ModelParaList, NumericalOutputs, Economies, InputsForOutputs,
#   FactorLabels,
#   JLLlist = NULL, GVARlist, WishBC, BRWlist, F2S
# )

