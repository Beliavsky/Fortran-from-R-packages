\name{X.pathway}
\alias{X.pathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a matrix of gene expressions in the presence of gene pathways}
\description{Generate a matrix of gene expressions in the presence of gene pathways (Scenario 2 of Emura et al. (2012)).}
\usage{X.pathway(n,p,q1,q2,rho1=0.5,rho2=0.5)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{the number of individuals (sample size)}
  \item{p}{the number of genes}
  \item{q1}{the number of genes in the first pathway}
  \item{q2}{the number of genes in the second pathway}
  \item{rho1}{the correlation coefficient for the first pathway}
  \item{rho2}{the correlation coefficient for the second pathway}
}
\details{
n by p matrix of gene expressions are generated. Correlation between columns (genes) is introduced to reflect the presence of gene pathways. The distribution of each column (gene) is standardized to have mean=0 and SD=1.
Two blocks of correlated genes (i.e., two gene pathways) are generated, where the first pathway include "q1" genes and the second pathway includes "q2" genes. The last p-q1-q2 genes are independent genes.
If two genes are correlated, the correlation is 0.5 (can be changed by option "rho1" and "rho2").
Details are referred to p.4 of Emura et al. (2012).
This deta generation scheme was used in the simulations of Emura et al. (2012), Emura and Chen (2016) and
Emura et al. (2017).
}
\value{
 \item{X}{n by p matrix of gene expressions}
}
\references{
Emura T, Chen YH, Chen HY (2012). Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models.
 PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring:
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57

Emura T, Nakatochi M, Matsui S, Michimae H, Rondeau V (2017) Personalized dynamic prediction of death according to
tumour progression and high-dimensional genetic factors: meta-analysis with a joint model,
Stat Methods Med Res, doi:10.1177/0962280216688032
}

\author{Takeshi Emura & Yi-Hau Chen}

\examples{
## generate 6 gene expressions from 10 individuals
X.pathway(n=10,p=6,q1=2,q2=2)

## generate 200 gene expressions and check the mean and SD
X.mat=X.pathway(n=200,p=100,q1=10,q2=10)
round( colMeans(X.mat),3 ) ## mean ~ 0 ##
round( apply(X.mat, MARGIN=2, FUN=sd),3) ## SD ~ 1 ##

## Change correlation coefficients by option "rho1" and "rho2" ##
X.mat=X.pathway(n=10000,p=6,q1=2,q2=2,rho1=0.2,rho2=0.8)
round( colMeans(X.mat),2 ) ## mean ~ 0 ##
round( apply(X.mat, MARGIN=2, FUN=sd),2) ## SD ~ 1 ##
round( cor(X.mat),1) ## Correlation matrix ##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene expression }
\keyword{ gene pathway }
