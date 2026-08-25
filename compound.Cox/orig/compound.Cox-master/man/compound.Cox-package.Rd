\name{compound.Cox-package}
\alias{compound.Cox-package}
\alias{compound.Cox}
\docType{package}
\title{
Univariate Feature Selection and Compound Covariate for Predicting Survival, Including Copula-Based Analyses for Dependent Censoring
}
\description{
Univariate feature selection and compound covariate methods under the Cox model with high-dimensional features (e.g., gene expressions).
 Available are survival data for non-small-cell lung cancer patients with gene expressions (Chen et al 2007 New Engl J Med),
 statistical methods in Emura et al (2012 PLoS ONE),
 Emura & Chen (2016 Stat Methods Med Res), and Emura et al. (2019 Comput Methods Programs Biomed).
 Algorithms for generating correlated gene expressions are also available.
 Estimation of survival functions via copula-graphic (CG) estimators is also implemented, which is useful for
 sensitivity analyses under dependent censoring (Yeh et al 2023 Biomedicines)
 and factorial survival analyses (Emura et al 2024 Stat Methods Med Res).
}
\details{
\tabular{ll}{
Package: \tab compound.Cox\cr
Type: \tab Package\cr
Version: \tab 3.33\cr
Date: \tab 2025-6-16\cr
License: \tab GPL-2\cr
}
}

\author{
Takeshi Emura, Hsuan-Yu Chen, Shigeyuki Matsui, Yi-Hau Chen;
Maintainer: Takeshi Emura <takeshiemura@gmail.com>
}
\references{
Chen HY, Yu SL, Chen CH, et al (2007). A Five-gene Signature and Clinical Outcome in Non-small-cell Lung Cancer, N Engl J Med 356: 11-20.

Emura T, Chen YH, Chen HY (2012). Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models.
 PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring:
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57

Emura T, Matsui S, Chen HY (2019). compound.Cox: Univariate Feature Selection and Compound Covariate for Predicting Survival,
Computer Methods and Programs in Biomedicine 168: 21-37

Matsui S (2006). Predicting Survival Outcomes Using Subsets of Significant Genes in Prognostic Marker Studies with Microarrays. BMC Bioinformatics: 7:156.

Yeh CT, Liao GY, Emura T (2023). Sensitivity analysis for survival prognostic prediction
with gene selection: a copula method for dependent censoring, Biomedicines 11(3):797.

Emura T, Ditzhaus M, Dobler D, Murotani K (2024). Factorial survival analysis for treatment effects under dependent censoring, Stat Methods Med Res 33(1):61-79.
}
\keyword{ package }
