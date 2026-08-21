\name{princomp.aplus}
\alias{princomp.aplus}
\alias{print.princomp.aplus}
\alias{plot.princomp.aplus}
\alias{plot.princomp.aplus}
\alias{predict.princomp.aplus}
\title{Principal component analysis for amounts in log geometry}
\description{
  A principal component analysis is done in the Aitchison geometry
  (i.e. ilt-transform). Some gimics simplify the
  interpretation of the computed components as perturbations of amounts.
}
\usage{
\method{princomp}{aplus}(x,\dots,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),
                           robust=getOption("robust"))
\method{print}{princomp.aplus}(x,\dots)
\method{plot}{princomp.aplus}(x,y=NULL,\dots, npcs=min(10,length(x$sdev)),
               type=c("screeplot","variance","biplot","loadings","relative"),
               main=NULL,scale.sdev=1)
\method{predict}{princomp.aplus}(object,newdata,\dots)
}
\arguments{
  \item{x}{an aplus dataset (for princomp) or a result from
    princomp.aplus}
  \item{y}{not used}
  \item{scores}{a logical indicating whether scores should be computed or not}
  \item{npcs}{the number of components to be drawn in the scree plot}
  \item{type}{type of the plot: \code{"screeplot"} is a lined screeplot,
    \code{"variance"} is a boxplot-like screeplot, \code{"biplot"} is a
    biplot, \code{"loadings"} displays the loadings as a 
    \code{\link{barplot.acomp}}}
  \item{scale.sdev}{the multiple of sigma to use when plotting the loadings}
  \item{main}{title of the plot}
  \item{object}{a fitted princomp.aplus object}
  \item{newdata}{another amount dataset of class aplus}
  \item{\dots}{further arguments to pass to internally-called functions}
  \item{covmat}{provides the covariance matrix to be used for the
    principle component analysis}
  \item{center}{provides the  be used for the computation of scores}
  \item{robust}{Gives the robustness type for the calculation of the
    covariance matrix. See \code{\link{var.rmult}} for details.}
}
\value{
  \code{princomp} gives an object of type
  \code{c("princomp.acomp","princomp")} with the following content:
  \item{sdev}{the standard deviation of the principal components}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix which
    columns contain the eigenvectors). This is of class
    \code{"loadings"}. }
  \item{center}{the ilt-transformed vector of means used to center the dataset}
  \item{Center}{the \code{\link{aplus}} vector of means used to center the dataset}
  \item{scale}{the scaling applied to each variable}
  \item{n.obs}{number of observations}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied data
    on the principal components. Scores are coordinates in a basis given by the principal
    components and thus not compositions}
  \item{call}{the matched call}
  \item{na.action}{not clearly understood}
  \item{Loadings}{vectors of amounts that represent a perturbation with the
    vectors represented by the loadings of each of the factors}
  \item{DownLoadings}{vectors of amounts that represent a perturbation with the
    inverses of the vectors represented by the loadings of each of the
    factors}
  \code{predict} returns a matrix of scores of the observations in the
  \code{newdata} dataset
  \cr.
  The other routines are mainly called for their side effect of plotting or
  printing and return the object \code{x}.   
}

\details{
  As a metric euclidean space, the positive real space described in
  \code{\link{aplus}} has its own
  principal component analysis, that can be performed either in terms of the
  covariance matrix or the correlation matrix. However, since all parts in a composition 
  or in an amount vector share a natural scaling, they do not need the 
  standardization (which in fact would produce a loss of important information). 
  For this reason, \code{princomp.aplus} works on the covariance matrix.
  \cr
  To aid the interpretation we added some extra functionality to a
  normal \code{princomp(ilt(x))}. First of all the result contains as
  additional information the amount representation of
  returned vectors in the space of the data: the center as an amount
  \code{Center}, and the loadings in terms of amounts to perturbe
  with, either positively
  (\code{Loadings}) or negatively (\code{DownLoadings}). The Up- and
  DownLoadings are normalized to the number of parts  
  and not to one to simplify the interpretation. A value of about one
  means no change in the specific component. 
  \cr
  The plot routine provides screeplots (\code{type = "s"},\code{type=
    "v"}), biplots (\code{type = "b"}), plots of the effect of
  loadings (\code{type = "b"}) in \code{scale.sdev*sdev}-spread, and
  loadings of pairwise (log-)ratios (\code{type = "r"}).
  \cr
  The interpretation of a screeplot does not differ from ordinary
  screeplots. It shows the eigenvalues of the covariance matrix, which
  represent the portions of variance explained by the principal
  components. 
  \cr
  The interpretation of the the biplot uses, additionally to the
  classical one, a compositional concept: The
  differences between two arrowheads can be interpreted as log-ratios
  between the two components represented by the arrows. 
  \cr
  The amount loading plot is introduced with this
  package. The loadings of all component can be seen as an orthogonal basis
  in the space of \code{\link{ilt}}-transformed data. These vectors are displayed by a barplot with
  their corresponding amounts. A portion of one means no change of this
  part. This is equivalent to a zero loading in a real principal component analysis. 
  \cr
  The loadings plot can work in two different modes. If
  \code{scale.sdev} is set to \code{NA} it displays the amount vector
  being represented by the unit vector of loadings in the ilt-transformed space. If
  \code{scale.sdev} is numeric we use this amount vector scaled by the
  standard deviation of the respective component. 
  \cr
  The relative plot displays the \code{\link{relativeLoadings}} as a
  barplot. The deviation from a unit bar shows the effect of each principal component 
  on the respective ratio. The
  interpretation of the ratios plot may only be done in an Aitchison-compositional framework 
  (see \code{\link{princomp.acomp}}).  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilt}},\code{\link{aplus}}, \code{\link{relativeLoadings}} 
  \code{\link{princomp.acomp}}, \code{\link{princomp.rplus}},
  \code{\link{barplot.aplus}}, \code{\link{mean.aplus}},
}

\examples{
data(SimulatedAmounts)
pc <- princomp(aplus(sa.lognormals5))
pc
summary(pc)
plot(pc)           #plot(pc,type="screeplot")
plot(pc,type="v")
plot(pc,type="biplot")
plot(pc,choice=c(1,3),type="biplot")
plot(pc,type="loadings")
plot(pc,type="loadings",scale.sdev=-1) # Downward
plot(pc,type="relative",scale.sdev=NA) # The directions
plot(pc,type="relative",scale.sdev=1) # one sigma Upward 
plot(pc,type="relative",scale.sdev=-1) # one sigma Downward
biplot(pc)
screeplot(pc)
loadings(pc)
relativeLoadings(pc,mult=FALSE)
relativeLoadings(pc)
relativeLoadings(pc,scale.sdev=1)
relativeLoadings(pc,scale.sdev=2)

pc$Loadings
pc$DownLoadings
barplot(pc$Loadings)
pc$sdev^2
cov(predict(pc,sa.lognormals5))
}
\keyword{multivariate}
