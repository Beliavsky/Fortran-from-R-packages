\name{acomp}
\alias{acomp}
\title{Aitchison compositions}
\description{
  A class providing the means to analyse compositions in the
  philosophical framework of the Aitchison Simplex.
}
\usage{
  acomp(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,
          detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{parts}{vector containing the indices xor names of the columns to be used}
  \item{total}{the total amount to be used, typically 1 or 100}
  \item{warn.na}{should the user be warned in case of NA,NaN or 0
    coding different types of missing values?}
  \item{detectionlimit}{a number, vector or matrix of positive
    numbers giving the detection limit of all values, all columns or
    each value, respectively}
  \item{BDL}{the code for 'Below Detection Limit' in X}
  \item{SZ}{the code for 'Structural Zero' in X}
  \item{MAR}{the code for 'Missing At Random' in X}
  \item{MNAR}{the code for 'Missing Not At Random' in X}
}
\value{
  a vector of class \code{"acomp"} representing one closed composition
  or a matrix of class \code{"acomp"} representing
  multiple closed compositions each in one row.

}

\details{
  Many multivariate datasets essentially describe amounts of D different
  parts in a whole. This has some important implications justifying to
  regard them as a scale for its own, called a
  composition. This scale  was in-depth analysed by Aitchison
  (1986) and the functions around the class \code{"acomp"} follow his
  approach.\cr
  Compositions have some important properties: Amounts are always
  positive. The amount of every part is limited to the whole. The
  absolute amount of the whole is noninformative since it is typically due
  to artifacts on the measurement procedure. Thus only relative changes
  are relevant. If the relative amount of one part
  increases, the amounts of other parts must decrease, introducing
  spurious anticorrelation (Chayes 1960), when analysed directly. Often
  parts (e.g H2O, Si) are missing in the dataset leaving the total
  amount unreported and longing for analysis procedures avoiding
  spurious effects when applied to such subcompositions. Furthermore, 
  the result of an analysis should be indepent of the units (ppm, g/l, vol.\%, mass.\%, molar
  fraction) of the dataset.
  \cr
  From these properties Aitchison showed that the
  analysis should be based on ratios or log-ratios only. He introduced 
  several transformations (e.g. \code{\link{clr}},\code{\link{alr}}),
  operations (e.g. \code{\link{perturbe}}, \code{\link{power.acomp}}),
  and a distance (\code{\link{dist}}) which are compatible
with these
  properties. Later it was found that the set of compostions equipped with
  perturbation as addition and power-transform as scalar multiplication
  and the \code{\link{dist}} as distance form a D-1 dimensional
  euclidean vector space (Billheimer, Fagan and Guttorp, 2001), which 
  can be mapped isometrically to a usual real vector space by \code{\link{ilr}} 
  (Pawlowsky-Glahn and Egozcue, 2001).
  \cr
  The general approach in analysing acomp objects is thus to perform
  classical multivariate analysis on clr/alr/ilr-transformed coordinates
  and to backtransform or display the results in such a way that they
  can be interpreted in terms of the original compositional parts.   
  \cr
  A side effect of the procedure is to force the compositions to sum up to a
  \var{total}, which is done by the closure operation \code{\link{clo}} . 
}
\section{Missing Policy}{
  The policy of treatment of zeroes, missing values and values 
  below detecion limit is explained in depth in \link{compositions.missing}. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{clr}},\code{\link{rcomp}}, \code{\link{aplus}},
  \code{\link{princomp.acomp}}, 
  \code{\link{plot.acomp}}, \code{\link{boxplot.acomp}},
  \code{\link{barplot.acomp}}, \code{\link{mean.acomp}},
  \code{\link{var.acomp}}, \code{\link{variation.acomp}},
  \code{\link{cov.acomp}}, \code{\link{msd}}
}

\references{
Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Billheimer, D., P. Guttorp, W.F. and Fagan (2001) Statistical interpretation of species composition,
  \emph{Journal of the American Statistical Association}, \bold{96} (456), 1205-1214\cr

  Chayes, F. (1960). On correlation between variables of constant
  sum. Journal of Geophysical Research 65~(12), 4185--4193.

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr

  Pawlowsky-Glahn, V. (2003) Statistical modelling on coordinates. In: 
  Thi\'o-Henestrosa, S. and Mart\'in-Fern\'andez, J.A. (Eds.)
  \emph{Proceedings of the 1st International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-111-X, \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr

  Mateu-Figueras, G. and Barcel\'o-Vidal, C. (Eds.)
  \emph{Proceedings of the 2nd International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-222-1, \url{https://ima.udg.edu/Activitats/CoDaWork05/}\cr

    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
}
\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
}
\keyword{classes}
