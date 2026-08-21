\name{princomp.rcomp}
\alias{princomp.rcomp}
\alias{print.princomp.rcomp}
\alias{plot.princomp.rcomp}
\alias{predict.princomp.rcomp}
\title{Principal component analysis for real compositions}
\description{
  A principal component analysis is done in real geometry
  (i.e. cpt-transform) of the simplex. Some gimics simplify the
  interpretation of the obtained components.  
}
\usage{
 \method{princomp}{rcomp}(x,\dots,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),
                           robust=getOption("robust"))
 \method{print}{princomp.rcomp}(x,\dots)
 \method{plot}{princomp.rcomp}(x,y=NULL,\dots,npcs=min(10,length(x$sdev)),
        type=c("screeplot","variance","biplot","loadings","relative"),
        main=NULL,scale.sdev=1)
 \method{predict}{princomp.rcomp}(object,newdata,\dots)
}
\arguments{
  \item{x}{an rcomp dataset (for princomp) or a result from
    princomp.rcomp}
  \item{y}{not used}
  \item{scores}{a logical indicating whether scores should be computed or not}
  \item{npcs}{the number of components to be drawn in the scree plot}
  \item{type}{type of the plot: \code{"screeplot"} is a lined screeplot,
    \code{"variance"} is a boxplot-like screeplot, \code{"biplot"} is a
    biplot, \code{"loadings"} displays the loadings as a 
    \code{\link[graphics]{barplot}}}
  \item{scale.sdev}{the multiple of sigma to use when plotting the loadings}
  \item{main}{title of the plot}
  \item{object}{a fitted princomp.rcomp object}
  \item{newdata}{another compositional dataset of class rcomp}
  \item{\dots}{further arguments to pass to internally-called functions}
  \item{covmat}{provides the covariance matrix to be used for the
    principle component analysis}
  \item{center}{provides the  be used for the computation of scores}
  \item{robust}{Gives the robustness type for the calculation of the
    covariance matrix. See \code{\link{var.rmult}} for details.}
}
\value{
  \code{princomp} gives an object of type
  \code{c("princomp.rcomp","princomp")} with the following content:
%  \itemize{
    \item{sdev}{the standard deviation of the principal components.}
    \item{loadings}{the matrix of variable loadings (i.e., a matrix which
      columns contain the eigenvectors). This is of class
      \code{"loadings"}. The last eigenvalue is removed since it should
      contain the irrelevant scaling.}
    \item{Loadings}{the loadings as an rmult-object}
    \item{center}{the cpt-transformed vector of means used to center the dataset}
    \item{Center}{the \code{\link{rcomp}} vector of means used to center the dataset}
    \item{scale}{the scaling applied to each variable}
    \item{n.obs}{number of observations}
    \item{scores}{if \code{scores = TRUE}, the scores of the supplied data
      on the principal components. Scores are coordinates in a basis given by the principal
      components and thus not compositions}
    \item{call}{the matched call}
    \item{na.action}{not clearly understood}
%  }
  \code{predict} returns a matrix of scores of the observations in the
  \code{newdata} dataset.
  The other routines are mainly called for their side effect of plotting or
  printing and return the object \code{x}.   
}
\details{
  Mainly a \code{princomp(cpt(x))} is performed.  To avoid confusion, the
  meaningless last principal component is removed.
  \cr
  The plot routine provides screeplots (\code{type = "s"},\code{type=
    "v"}), biplots (\code{type = "b"}), plots of the effect of
  loadings (\code{type = "b"}) in \code{scale.sdev*sdev}-spread, and
  loadings of pairwise differences (\code{type = "r"}).
  \cr
  The interpretation of a screeplot does not differ from ordinary
  screeplots. It shows the eigenvalues of the covariance matrix, which
  represent the portions of variance explained by the principal
  components. 
  \cr
  The interpretation of the biplot strongly differs from a classical one.
  The relevant variables are not the arrows drawn (one for each component), 
   but rather the links (i.e., the differences) between two
  arrow heads, which represents the difference between the two
  components represented by the arrows, or the transfer of mass from 
  one to the other.
  \cr
  The compositional loading plot is more or less a standard
  one. The loadings are displayed by a barplot as positve and
  negative changes of amounts.
  \cr
  The loading plot can work in
  two different modes: If
  \code{scale.sdev} is set to \code{NA} it displays the composition
  being represented by the unit vector of loadings in cpt-transformed space. If
  \code{scale.sdev} is numeric we use this composition scaled by the
  standard deviation of the respective component. 
  \cr
  The relative plot displays the \code{\link{relativeLoadings}} as a
  barplot. The deviation from a unit bar shows the effect of each
  principal component on the respective differences. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cpt}},\code{\link{rcomp}}, \code{\link{relativeLoadings}} 
  \code{\link{princomp.acomp}}, \code{\link{princomp.rplus}},
}
\examples{
data(SimulatedAmounts)
pc <- princomp(rcomp(sa.lognormals5))
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

pc$sdev^2
cov(predict(pc,sa.lognormals5))
}
\keyword{multivariate}
