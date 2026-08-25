\name{X.tag}
\alias{X.tag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a matrix of gene expressions in the presence of tag genes}
\description{Generate a matrix of gene expressions in the presence of tag genes (Scenario 1 of Emura et al. (2012)).
}
\usage{
X.tag(n, p, q, s = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{the number of individuals (sample size)}
  \item{p}{the number of genes}
  \item{q}{the number of non-null genes}
  \item{s}{the number of null genes correlated with a non-null gene (tag)}
}
\details{n by p matrix of gene expressions are generated. 
Correlation between columns is introduced to reflect the presence of tag genes. 
The distribution of each column is standardized to have mean=0 and SD=1.
If two genes are correlated, the correlation is 0.5. Otherwise, the correlation is 0. 
Details are referred to p.4 of Emura et al. (2012).
This deta generation scheme was used in the simulations of Emura et al. (2012) and Emura and Chen (2016).
}
\value{
 \item{X}{n by p matrix of gene expressions}
}
\references{
Emura T, Chen YH, Chen HY (2012). Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models.
 PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.
}
\author{Takeshi Emura & Yi-Hau Chen}

\examples{
X.mat=X.tag(n=200,p=100,q=10,s=4)
round( colMeans(X.mat),3 ) ## mean ~ 0 ##
round( apply(X.mat, MARGIN=2, FUN=sd),3) ## SD ~ 1 ##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene expression }
\keyword{ tag gene }
