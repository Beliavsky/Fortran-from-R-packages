\name{uni.selection}
\alias{uni.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Univariate feature selection based on univariate significance tests}
\description{
This function performs univariate feature selection using significance tests (Wald tests or score tests) based on 
association between individual features and survival. 
Features are selected if their P-values are less than a given threshold (P.value). 
}
\usage{
uni.selection(t.vec, d.vec, X.mat, P.value=0.001,K=10,score=TRUE,d0=0,
                       randomize=FALSE,CC.plot=FALSE,permutation=FALSE,M=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators (1=death, 0=censoring)}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{P.value}{A threshold for selecting features}
  \item{K}{The number of cross-validation folds}
  \item{score}{If TRUE, the score tests are used; if not, the Wald tests are used}
  \item{d0}{A positive constant to stabilize the variance of score statistics (Witten & Tibshirani 2010)}
  \item{randomize}{If TRUE, randomize patient ID's before cross-validation}
  \item{CC.plot}{If TRUE, the compound covariate (CC) predictors are plotted}
  \item{permutation}{If TRUE, the FDR is computed by a permutation method (Witten & Tibshirani 2010; Emura et al. 2019).}
  \item{M}{The number of permutations to calculate the FDR}
}
\details{
The cross-validated likelihood (CVL) value is computed for selected features (Matsui 2006; Emura et al. 2019). 
A high CVL value corresponds to a better predictive ability of selected features. 
Hence, the CVL value can be used to find the optimal set of features. 
The CVL value is computed by a K-fold cross-validation, where the number K can be chosen by user. 
The false discovery rate (FDR) is also computed by a formula and a permutation test (if "permutation=TRUE"). 
The RCVL1 and RCVL2 are "re-substitution" CVL values and provide upper control limits for the CVL value.
If the CVL value is less than RCVL1 and RCVL2 values, the CVL value would be in-control. 
On the other hand, if the CVL value exceeds either RCVL1 or RCVL2 value, then 
the CVL may be computed again after changing the sample allocation.  
}
\value{
  \item{gene }{Gene symbols}
  \item{beta }{Estimated regression coefficients}
  \item{Z }{Z-values for significance tests}
  \item{P }{P-values for significance tests}
  \item{CVL }{The value of CVL, RCVL1, and RCVL2 (Emura et al. 2019)}
  \item{Genes }{The number of genes, the number of selected genes, and the number of falsely selected genes}
  \item{FDR }{False discovery rate (by a formula or a permutation method)}
}
\references{
Emura T, Matsui S, Chen HY (2019). compound.Cox: Univariate Feature Selection and Compound Covariate for Predicting Survival, 
Computer Methods and Programs in Biomedicine 168: 21-37.

Matsui S (2006). Predicting Survival Outcomes Using Subsets of Significant Genes in Prognostic Marker Studies with Microarrays. BMC Bioinformatics: 7:156.

Witten DM, Tibshirani R (2010) Survival analysis with high-dimensional covariates. Stat Method Med Res 19:29-51
}
\author{Takeshi Emura}
\examples{
data(Lung)
t.vec=Lung$t.vec[Lung$train==TRUE]
d.vec=Lung$d.vec[Lung$train==TRUE]
X.mat=Lung[Lung$train==TRUE,-c(1,2,3)]
uni.selection(t.vec, d.vec, X.mat, P.value=0.05,K=5,score=FALSE)
## the outputs reproduce Table 3 of Emura and Chen (2016) ## 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univariate Cox regression }
\keyword{ compound covariate }
\keyword{ cross-validated partial likelihood }
\keyword{ false discovery rate }
\keyword{ feature selection }% __ONLY ONE__ keyword per line
