\name{otrimlesimg}

\alias{otrimlesimg}
\alias{summary.otrimlesimgdens}
\alias{print.summary.otrimlesimgdens}
\alias{plot.summary.otrimlesimgdens}


\title{Adequacy approach for number of clusters for OTRIMLE}

\description{
  \code{otrimlesimg} computes Optimally Tuned Robust Improper Maximum
  Likelihood Clustering
  (OTRIMLE), see \code{\link{otrimle}} for a range of values of the
  number of clusters, and also for artificial datasets simulated from
  the model parameters estimated on the original data. The
  \code{summary}-methods present and evaluate the results so that a
  smallest adequate number of clusters can be found as the smallest one
  for which the value of the density-based cluster quality statistics Q
  on the original data 
  is compatible with its distribution on the artificial datasets with
  the same number of clusters, see Hennig and Coretto 2021 for details.
}

\usage{
otrimlesimg(dataset, G=1:6, multicore=TRUE,
ncores=detectCores(logical=FALSE)-1, erc=20, beta0=0, simruns=20,
sim.est.logicd=FALSE, 
monitor=1)

\method{summary}{otrimlesimgdens}(object, noisepenalty=0.05 , sdcutoff=2
, ...)

\method{print}{summary.otrimlesimgdens}(x, ...)

\method{plot}{summary.otrimlesimgdens}(x , plot="criterion", penx=NULL,
peny=NULL, pencex=1, cutoff=TRUE, ylim=NULL, ...)

}



\arguments{
  \item{dataset}{something that can be coerced into an
    observations times variables matrix. The dataset.}
  \item{G}{vector of integers (normally starting from 1). Numbers of
    clusters to be considered.}
  \item{multicore}{logical. If \code{TRUE}, parallel computing is used
    through the function \code{\link{mclapply}} from package
    \code{parallel}; read warnings there if you intend to use this; it
    won't work on Windows.}
  \item{ncores}{integer. Number of cores for parallelisation.}
  \item{erc}{A number larger or equal than one specifying the maximum
    allowed ratio between within-cluster covariance matrix
    eigenvalues. See \code{\link{otrimle}}.}
  \item{beta0}{A non-negative constant, penalty term for noise, to be
    passed as \code{beta} to  \code{\link{otrimle}}, see documentation
    there.}
  \item{simruns}{integer. Number of replicate artificial datasets drawn
    from each model.} 
  \item{sim.est.logicd}{logical. If \code{TRUE}, the logarithm of the improper
    constant density \code{logicd}, see \code{\link{otrimle}}, is
    re-estimated when running \code{\link{otrimle}} on the artificial
    datasets. Otherwise the value estimated on the original data is
    taken as fixed. \code{TRUE} requires much longer computation time,
    but can be seen as generating more realistic variation.}
  \item{monitor}{0 or 1. If 1, progress messages are printed on screen.}
  \item{noisepenalty}{number between 0 and 1. \code{p_0} in Hennig and
    Coretto (2021); normally small. The method prefers to treat a
    proportion of \code{<=noisepenalty} of points as outliers to adding a
    cluster.}
  \item{sdcutoff}{numerical. \code{c} in formula (7) in Hennig and
    Coretto (2021). A clustering is treated as adequate if its value of
    the density-based cluster quality measure Q calibrated (i.e.,
    mean/sd-standardised) by the values on the artificial datasets
    generated from the estimated model is \code{<=sdcutoff}.}
  \item{plot}{\code{"criterion"} or \code{"noise"}, see details.}
  \item{penx}{\code{FALSE, NULL}, or numerical. x-coordinate from where
    the simplicity ordering of clustering is given (as test in the
    plot). If \code{FALSE}, this is not added to the plot. If
    \code{NULL} a default guess is made for a good position (which
    doesn't always work well).}
  \item{peny}{\code{NULL}, or numerical. x-coordinate from where
    the simplicity ordering of clustering is given (as test in the
    plot). If
    \code{NULL}, a default guess is made for a good position (which
    doesn't always work well).}
  \item{pencex}{numeric. Magnification factor (parameter \code{cex} to
    be passed on to \code{\link{legend}}) for simplicity ordering, see
    parameter \code{penx}.}
  \item{cutoff}{logical. If \code{TRUE}, the \code{"criterion"}-plot
    shows the cutoff value below which numbers of clusters are adequate,
    see details.}
  \item{ylim}{vector of two numericals, range of the y-axis to be passed
    on to \code{\link{plot}}. If \code{NULL}, the range is chosen
    automatically (but can be different from the \code{\link{plot}}
    default).}
  \item{object}{an object of class \code{'otrimlesimgdens'} obtained
    from calling \code{\link{otrimlesimg}}
  }
  \item{x}{an object of class \code{'summary.otrimlesimgdens'} obtained
    from calling \code{\link{summary}} function over an object of class 
    \code{'otrimlesimgdens'} obtained from calling \code{\link{otrimlesimg}}.
  }
  \item{...}{optional parameters to be passed on to \code{\link{plot}}.}
}



\details{
  The method is fully described in Hennig and Coretto
  (2021). The required tuning constants for choosing an optimal number
  of clusters, the smallest percentage of additional noise that the user
  is willing to trade in for adding another cluster (\code{p_0} in the
  paper, \code{noisepenalty} here) and the critical value (\code{c} in
  the paper, \code{sdcutoff} here) for adequacy of the standardised
  density based quality measure Q are provided to the summary function,
  which is required to choose the best (simplest adequate) number of
  clusters.

  The plot function \code{plot.summary.otrimlesimgdens} can produce two
  plots. If \code{plot="criterion"}, the standardised density-based
  cluster quality
  measure Q is plotted against the number of clusters. The values for
  the simulated artificial datasets are points, the values for the
  original dataset are given as line type. If \code{cutoff="TRUE"}, the
  critical values (see above) are added as red crosses; a number of
  clusters is adequate if the value of the original data is below the
  critical value, i.e., Q is not significantly larger than for the
  artificial datasets generated from the fitted model. Using
  \code{penx}, the ordered numbers of clusters from the simplest to the
  least simple can also be indicated in the plot, where simplicitly is
  defined as the number of clusters plus the estimated noise proportion
  divided by \code{noisepenalty}, see above. The chosen number of
  clusters is the simplest adequate one, meaning that a low number of
  clusters and a low noise proportion are preferred.

  If \code{plot="noise"}, the noise proportion (black) and the
  simplicity (red) are plotted against the numnber of clusters.
}

\value{
  \code{otrimlesimg} returns a list of type \code{"otrimlesimgdens"}
  containing the components \code{result, simresult, simruns}.   

  \item{result}{output object of \code{\link{otrimleg}} (list of results on
    original data) run with the parameters provided to
    \code{\link{otrimlesimg}}.}
  \item{simresult}{list of length \code{simruns} of output objects of
    \code{\link{otrimleg}} for all the simulated artificial datasets.}
  \item{simruns}{input parameter \code{simruns}.}

  \code{summary.otrimlesimgdens} returns a list of type
  \code{"summary.otrimlesimgdens"} with components \code{G, simeval,
    ssimruns, npr, nprdiff, logicd, denscrit, peng,
    penorder, bestG, sdcutoff, bestresult,
    cluster}. \code{simruns} 

  \item{G}{\code{otrmlesimg} input parameter \code{G} (numbers of
    clusters).}
  \item{simeval}{list with components \code{denscrit,
      meandens, sddens, standens, errors}, defined below.}
  \item{ssimruns}{\code{otrmlesimg} input parameter \code{simruns}.}
  \item{npr}{vector of estimated noise proportions on the original data
    for all numbers of clusters, \code{exproportion[1]} from
    \code{\link{otrimle}}.}
  \item{nprdiff}{vector for all numbers of clusters of differences
    between estimated smallest cluster proportion and noise proportion
    on the original data.}
  \item{logicd}{vector of logs of improper constant density values on
    the original data for all numbers of clusters.}
  \item{denscrit}{vector over all numbers of clusters of density-based
    cluster quality statistics Q
    on original data as provided by the \code{measure}-component of
    \code{\link{kerndensmeasure}}.}
  \item{peng}{vector of simplicity values (see Details) over all numbers
    of clusters.}
  \item{penorder}{simplicity order of number of clusters.}
  \item{bestG}{best (i.e., most simple adequate) number of clusters.}
  \item{sdcutoff}{input parameter \code{sdcutoff}.}
  \item{result}{output of \code{\link{otrimle}} for the best number of
    clusters \code{bestG}.}
  \item{cluster}{clustering vector for the best number of
    clusters \code{bestG}. \code{0} corresponds to noise/outliers.}

  Components of \code{summary.otrimlesimgdens} output component
  \code{simeval}:

  \item{denscritmatrix}{maximum number of clusters times \code{simruns} matrix
    of \code{denscrit}-vectors for all clusterings on simulated data.}
  \item{meandens}{vector over numbers of clusters of robust estimator of
    the mean of \code{denscrit} over simulated datasets, computed by
    \code{\link[robustbase]{scaleTau2}}.}    
  \item{sddens}{vector over numbers of clusters of robust estimator of
    the standard deviation of \code{denscrit} over simulated datasets,
    computed by \code{\link[robustbase]{scaleTau2}}.}    
  \item{standens}{vector over numbers of clusters of \code{denscrit} of
    the original data standardised by \code{meandens} and
    \code{sddens}.}

  \item{errors}{vector over numbers of clusters of numbers of times that
    \code{\link{otrimle} led to an error.}
    
    \code{plot.summary.otrimlesimgdens} will return the output
    of \code{\link{par}()} before anything was changed by the plot
    function. 
  }
}


  \references{
    Coretto, P.  and C. Hennig (2016).  Robust improper maximum
    likelihood: tuning, computation, and a comparison with other
    methods for robust Gaussian clustering.  \emph{Journal of the American
    Statistical Association}, Vol. 111(516), pp. 1648-1659.
    \doi{10.1080/01621459.2015.1100996}

    P. Coretto and C. Hennig (2017).  Consistency, breakdown
    robustness, and algorithms for robust improper maximum likelihood
    clustering. \emph{Journal of Machine Learning Research}, Vol. 18(142),
    pp. 1-39.  \url{https://jmlr.org/papers/v18/16-382.html}

    Hennig, C. and P.Coretto (2021). An adequacy approach for deciding the
    number of clusters for OTRIMLE robust Gaussian mixture based
    clustering. To appear in \emph{Australian and New Zealand Journal of
    Statistics}, \url{https://arxiv.org/abs/2009.00921}.
  }

  
  \author{Christian Hennig
    \email{christian.hennig@unibo.it}
    \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
   }


  \seealso{
    \code{\link{otrimle}}, \code{\link{rimle}}, \code{\link{otrimleg}},
    \code{\link{kerndensmeasure}} 
  }

\examples{
## otrimlesimg is computer intensive, so only a small data subset
## is used for speed.
data(banknote)
selectdata <- c(1:30,101:110,117:136,160:161)
set.seed(555566)
x <- banknote[selectdata,5:7]
   
## simruns=2 chosen for speed. This is not recommended in practice. 
obanknote <- otrimlesimg(x,G=1:2,multicore=FALSE,simruns=2,monitor=0)
sobanknote <- summary(obanknote)
print(sobanknote)
plot(sobanknote,plot="criterion",penx=1.4)
plot(sobanknote,plot="noise",penx=1.4)
plot(x,col=sobanknote$cluster+1,pch=c("N","1","2")[sobanknote$cluster+1])
}

\keyword{cluster}
\keyword{robust}

