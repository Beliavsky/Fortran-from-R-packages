\name{compound.reg}
\alias{compound.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compound shrinkage estimation under the Cox model 
}
\description{
This function implements the "compound shrinkage estimator" to calculate the regression coefficients
 of the Cox model, which was proposed by Emura, Chen & Chen (2012).
 The method is a variant of the Cox partial likelihood estimator such that
 the regression coefficients are mixed with the univariate Cox regression estimators.
 The resultant estimator is applicable even when the number of covariates is greater than the number of
 samples (the p>n setting).
 The standard errors (SEs) are calculated based on the asymptotic theory (see Emura et al., 2012).    
}
\usage{
compound.reg(t.vec, d.vec, X.mat, K = 5, delta_a = 0.025, a_0 = 0, var = FALSE,
plot=TRUE, randomize = FALSE, var.detail = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{K}{The number of cross validation folds, K=n corresponds to a leave-one-out cross validation (default=5)}
  \item{delta_a}{The step size for a grid search for the maximum of the cross-validated likelihood (default=0.025)}
  \item{a_0}{The starting value of a grid search for the maximum of the cross-validated likelihood (default=0)}
  \item{var}{If TRUE, the standard deviations and confidence intervals are given (default=FALSE, to reduce the computational cost)}
  \item{plot}{If TRUE, the cross validated likelihood curve and its maximized point are drawn}
  \item{randomize}{If TRUE, randomize the subject ID's so that the subjects in the cross validation folds
are randomly chosen. Otherwise, the cross validation folds are constructed in the ascending sequence}
  \item{var.detail}{Detailed information about the covariance matrix, which is mainly used for theoretical purposes. 
Please consult Takeshi Emura for more details (default=FALSE)}
}
\details{
K=5 cross validation is recommended for computational efficiency, though the results appear to be robust
against the choice of the number K. If the number of covariates is greater than 200, the computational time becomes very long.
 In such a case, the univariate pre-selection is recommended to reduce the number of covariates.
}
\value{
\item{a }{An optimized value of the shrinkage parameter (0<=a<=1)}
\item{beta }{Estimated regression coefficients}
\item{SE }{Standard errors for estimated regression coefficients}
\item{Lower95CI }{Lower ends of 95 percent confidence intervals (beta_hat-1.96*SE)}
\item{Upper95CI }{Upper ends of 95 percent confidence intervals (beta_hat+1.96*SE)}
\item{Sigma }{Covariance matrix for estimated regression coefficients}
\item{V}{Estimates of the information matrix (-[Hessian of the loglikelihood]/n)}
\item{Hessian_CV}{Second derivative of the cross-validated likelihood. Normally negative
since the cross-validated curve is concave}
\item{h_dot}{Derivative of Equation (8) of Emura et al. (2012) with respect to a shrinkage parameter "a"}
}
\references{
Emura T, Chen Y-H, Chen H-Y (2012) Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models.
 PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627
}
\author{Takeshi Emura & Yi-Hau Chen}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
### A simulation study ###
n=50 ### sample size
beta_true=c(1,1,0,0,0)
p=length(beta_true) 
t.vec=d.vec=numeric(n)
X.mat=matrix(0,n,p)

set.seed(1)
for(i in 1:n){
  X.mat[i,]=rnorm(p,mean=0,sd=1)
  eta=sum( as.vector(X.mat[i,])*beta_true )
  T=rexp(1,rate=exp(eta))
  C=runif(1,min=0,max=5)
  t.vec[i]=min(T,C)
  d.vec[i]=(T<=C)
}
compound.reg(t.vec,d.vec,X.mat,delta_a=0.1) 
### compare the estimates (beta) with the true value ###
beta_true

### Lung cancer data analysis (Emura et al. 2012 PLoS ONE) ###
data(Lung)
temp=Lung[,"train"]==TRUE
t.vec=Lung[temp,"t.vec"]
d.vec=Lung[temp,"d.vec"]
X.mat=as.matrix( Lung[temp,-c(1,2,3)] )
#compound.reg(t.vec=t.vec,d.vec=d.vec,X.mat=X.mat,delta_a=0.025) # time-consuming process
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compound covariate }
\keyword{ shrinkage estimation }% __ONLY ONE__ keyword per line
