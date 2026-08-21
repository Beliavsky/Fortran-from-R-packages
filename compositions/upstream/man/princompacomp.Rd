\name{princomp.acomp}
\alias{princomp.acomp}
\alias{print.princomp.acomp}
\alias{plot.princomp.acomp}
\alias{plot.princomp.acomp}
\alias{predict.princomp.acomp}
\title{Principal component analysis for Aitchison compositions}
\description{
  A principal component analysis is done in the Aitchison geometry
  (i.e. clr-transform) of the simplex. Some gimics simplify the
  interpretation of the computed components as compositional perturbations.  
}
\usage{
 \method{princomp}{acomp}(x,\dots,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),
                           robust=getOption("robust"))
 \method{print}{princomp.acomp}(x,\dots)
 \method{plot}{princomp.acomp}(x,y=NULL,\dots, npcs=min(10,length(x$sdev)),
        type=c("screeplot","variance","biplot","loadings","relative"),
        main=NULL,scale.sdev=1)
 \method{predict}{princomp.acomp}(object,newdata,\dots)

}
\arguments{
  \item{x}{a acomp-dataset (in princomp) or a result from
    princomp.acomp}
  \item{y}{not used}
  \item{scores}{a logical indicating whether scores should be computed or not}
  \item{npcs}{the number of components to be drawn in the scree plot}
  \item{type}{type of the plot: \code{"screeplot"} is a lined screeplot,
    \code{"variance"} is a boxplot-like screeplot, \code{"biplot"} is a
    biplot, \code{"loadings"} displays the loadings as a 
    \code{\link{barplot.acomp}}}
  \item{scale.sdev}{the multiple of sigma to use plotting the loadings}
  \item{main}{title of the plot}
  \item{object}{a fitted princomp.acomp object}
  \item{newdata}{another compositional dataset of class acomp}
  \item{\dots}{further arguments to pass to internally-called functions}
  \item{covmat}{provides the covariance matrix to be used for the
    principle component analysis}
  \item{center}{provides the  be used for the computation of scores}
  \item{robust}{Gives the robustness type for the calculation of the
    covariance matrix. See \code{\link{robustnessInCompositions}} for details.}
}
\value{
  \code{princomp} gives an object of type
  \code{c("princomp.acomp","princomp")} with the following content:
  \item{sdev}{the standard deviation of the principal components}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix which
    columns contain the eigenvectors). This is of class
    \code{"loadings"}. The last eigenvector is removed since it should
    contain the irrelevant scaling.}
  \item{center}{the clr-transformed vector of means used to center the dataset}
  \item{Center}{the \code{\link{acomp}} vector of means used to center the dataset}
  \item{scale}{the scaling applied to each variable}
  \item{n.obs}{number of observations}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied data
    on the principal components. Scores are coordinates in a basis given by the principal
    components and thus not compositions}
  \item{call}{the matched call}
  \item{na.action}{not clearly understood}
  \item{Loadings}{compositions that represent a perturbation with the
    vectors represented by the loadings of each of the factors }
  \item{DownLoadings}{compositions that represent a perturbation with the
    inverse of the vectors represented by the loadings of each of the
    factors}
  \code{predict} returns a matrix of scores of the observations in the
  \code{newdata} dataset
  \cr.
  The other routines are mainly called for their side effect of plotting or
  printing and return the object \code{x}.   
}

\details{
  As a metric euclidean space the Aitchison simplex has its own
  principal component analysis, that should be performed in terms of the
  covariance matrix and not in terms of the meaningless correlation
  matrix. 
  \cr
  To aid the interpretation we added some extra functionality to a
  normal \code{princomp(clr(x))}. First of all the result contains as
  additional information the compositional representation of the
  returned vectors in the space of the data: the center as a composition
  \code{Center}, and the loadings in terms of a composition to perturbe
  with, either positively
  (\code{Loadings}) or negatively (\code{DownLoadings}). The Up- and
  DownLoadings are normalized to the number of parts in the simplex
  and not to one to simplify the interpretation. A value of about one
  means no change in the specific component. To avoid confusion the
  meaningless last principal component is removed.
  \cr
  The \code{plot} routine provides screeplots (\code{type = "s"},\code{type=
    "v"}), biplots (\code{type = "b"}), plots of the effect of
  loadings (\code{type = "b"}) in \code{scale.sdev*sdev}-spread, and
  loadings of pairwise (log-)ratios (\code{type = "r"}).
  \cr
  The interpretation of a screeplot does not differ from ordinary
  screeplots. It shows the eigenvalues of the covariance matrix, which
  represent the portions of variance explained by the principal
  components. 
  \cr
  The interpretation of the biplot strongly differs from a classical one.
  The relevant variables are not the arrows drawn (one for each component), 
   but rather the links (i.e., the differences) between two
  arrow heads, which represents the log-ratio between the two
  components represented by the arrows. 
  \cr
  The compositional loading plot is introduced with this
  package. The loadings of all component can be seen as an orthogonal basis
  in the space of clr-transformed data. These vectors are displayed by a barplot with
  their corresponding composition. For a better
  interpretation the total of these compositons is set to the number of
  parts in the composition, such that a portion of one means no
  effect. This is similar to (but not exactly the same as) a zero loading in a real
  principal component analysis. 
  \cr
  The loadings plot can work in two different modes: if
  \code{scale.sdev} is set to \code{NA} it displays the composition
  beeing represented by the unit vector of loadings in the clr-transformed space. If
  \code{scale.sdev} is numeric we use this composition scaled by the
  standard deviation of the respective component. 
  \cr
  The relative plot displays the \code{\link{relativeLoadings}} as a
  barplot. The deviation from a unit bar shows the effect of each
  principal component on the respective ratio. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{clr}},\code{\link{acomp}}, \code{\link{relativeLoadings}} 
  \code{\link{princomp.aplus}}, \code{\link{princomp.rcomp}},
  \code{\link{barplot.acomp}}, \code{\link{mean.acomp}},
  \code{\link{var.acomp}} 
}

\references{
  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Aitchison, J. and M. Greenacre (2002) \emph{Biplots for Compositional
    Data} Journal of the Royal Statistical Society, Series C (Applied Statistics)
    \bold{51} (4) 375-392\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

}
\examples{
data(SimulatedAmounts)
pc <- princomp(acomp(sa.lognormals5))
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
p = predict(pc,sa.lognormals5)
cov(p)
}
\keyword{multivariate}
