library(testthat)
library(MultiATSM)

# Load inputs
LoadData("CM_2024")

Economies <- c("China", "Brazil")
GlobalVar <- c("Gl_Eco_Act")
DomVar <- c("Eco_Act")
N <- 2
t0_sample <- "01-01-2005"
tF_sample <- "01-12-2019"
DataFreq <- "Monthly"

# Inputs for a BRW model
WishBC <- TRUE # useful to get the risk factors, bias-corrected features will be tested below
BRWlist <- list(Cent_Measure = "Mean", gamma = 0.2, N_iter = 5, N_burn = 2, B = 2, check = FALSE)

GVARlist <- list(
  VARXtype = "unconstrained", W_type = "Sample Mean", t_First_Wgvar = "2005",
  t_Last_Wgvar = "2019", DataConnectedness = TradeFlows
)
JLLlist <- list(DomUnit = "China")


# 1) Test output structure - "JPS original"
test_that("Bias_Correc_VAR returns correct output structure (JPS model)", {
  ModelType <- "JPS original"

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  ATSMInputs <- InputsForOpt(t0_sample, tF_sample, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFreq, GVARlist, JLLlist, WishBC, BRWlist,
    verbose = FALSE
  )

  for (i in 1:length(Economies)) {
    RiskFactors_CS <- ATSMInputs[[Economies[i]]]$RiskFactors
    res <- Bias_Correc_VAR(ModelType, BRWlist, t(RiskFactors_CS), Economies, FactorLabels, verbose = FALSE)
    expect_type(res, "list")
    expect_true(all(c("K0Z_BC", "K1Z_BC", "SSZ_BC", "dist") %in% names(res)))
    K <- nrow(RiskFactors_CS)
    expect_equal(dim(res$K0Z_BC), c(K, 1))
    expect_equal(dim(res$K1Z_BC), c(K, K))
    expect_equal(dim(res$SSZ_BC), c(K, K))
  }
})



# 2) Test output structure - "GVAR multi"
test_that("Bias_Correc_VAR returns correct output structure (GVAR model)", {
  set.seed(1)
  ModelType <- "GVAR multi"

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  ATSMInputs <- InputsForOpt(t0_sample, tF_sample, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFreq, GVARlist, JLLlist, WishBC, BRWlist,
    verbose = FALSE
  )

  RiskFactors <- ATSMInputs$RiskFactors
  GVARinputs <- ATSMInputs$GVARinputs

  BRWlist$check <- TRUE
  BRWlist$B_check <- 10

  res <- Bias_Correc_VAR(ModelType, BRWlist, t(RiskFactors), Economies, FactorLabels, GVARinputs, verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("K0Z_BC", "K1Z_BC", "SSZ_BC", "dist") %in% names(res)))
  K <- nrow(RiskFactors)
  expect_equal(dim(res$K0Z_BC), c(K, 1))
  expect_equal(dim(res$K1Z_BC), c(K, K))
  expect_equal(dim(res$SSZ_BC), c(K, K))
  expect_type(res$dist, "double")
})


# 3) Test output structure - "JLL original"
test_that("Bias_Correc_VAR returns correct output structure (JLL model)", {
  set.seed(20)
  ModelType <- "JLL original"

  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  ATSMInputs <- InputsForOpt(t0_sample, tF_sample, ModelType, Yields, GlobalMacro, DomMacro,
    FactorLabels, Economies, DataFreq, GVARlist, JLLlist, WishBC, BRWlist,
    verbose = FALSE
  )

  RiskFactors <- ATSMInputs$RiskFactors
  GVARinputs <- NULL
  JLLinputs <- ATSMInputs$JLLinputs

  res <- Bias_Correc_VAR(ModelType, BRWlist, t(RiskFactors), Economies, FactorLabels, GVARinputs,
    JLLinputs,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_true(all(c("K0Z_BC", "K1Z_BC", "SSZ_BC", "dist") %in% names(res)))
  K <- nrow(RiskFactors)
  expect_equal(dim(res$K0Z_BC), c(K, 1))
  expect_equal(dim(res$K1Z_BC), c(K, K))
  expect_equal(dim(res$SSZ_BC), c(K, K))
})
