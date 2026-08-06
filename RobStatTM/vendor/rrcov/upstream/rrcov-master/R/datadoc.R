######
##  VT::03.08.2019
##
##
##  roxygen2::roxygenise("c:/Users/valen/OneDrive/MyRepo/R/rrcov", load_code=roxygen2:::load_installed, clean=TRUE)
##
#'
#'
#'
#' Computer Hardware
#'
#' A data set containing relative CPU performance data of 209 machines on 8 variables.
#'  The \code{rownames} are the vendor and model descriptions. Six of the variables
#'  are predictive, one (\code{PRP}) is the goal field and one (\code{ERP}) is the
#'  linear regression's guess. The estimated relative performance values were
#'  estimated by the authors using a linear regression method.  See their article
#'  (Ein-Dor and Feldmesser, CACM 4/87, pp 308-317) for more details on how the
#'  relative performance values were set.
#'
#' @name machines
#' @docType data
#' @usage data(machines)
#' @format A data frame with 209 rows and 8 variables
#' The variables are as follows:
#'
#' \itemize{
#'   \item MMIN: minimum main memory in kilobytes (integer)
#'   \item MMAX: maximum main memory in kilobytes (integer)
#'   \item CACH: cache memory in kilobytes (integer)
#'   \item CHMIN: minimum channels in units (integer)
#'   \item CHMAX: maximum channels in units (integer)
#'   \item PRP: published relative performance (integer)
#'   \item ERP: estimated relative performance from the original article (integer)
#' }
#'
#' @source \href{http://archive.ics.uci.edu/ml/datasets/Computer+Hardware?ref=datanews.io}{UCI Archive}
#'
#' @references
#'  Phillip Ein-Dor and Jacob Feldmesser (1987), Attributes of the performance
#'      of central processing units: A relative performance prediction model,
#'      \emph{Communications of the ACM}, \bold{30}, 4, pp 308-317.
#'
#'    M. Hubert, P. J. Rousseeuw and T. Verdonck (2009), Robust PCA for skewed data and
#'    its outlier map, \emph{Computational Statistics & Data Analysis}, \bold{53}, 2264--2274.
#'
#'
#' @examples
#'
#'  data(machines)
#'
#'  ## Compute the medcouple of each variable of the Computer hardware data
#'      data.frame(MC=round(apply(machines, 2, mc),2))
#'
#'  ## Plot a pairwise scaterplot matrix
#'      pairs(machines[,1:6])
#'
#'      mcd <- CovMcd(machines[,1:6])
#'      plot(mcd, which="pairs")
#'
#'  ##  Remove the rownames (too long)
#'      rownames(machines) <- NULL
#'
#'  ## Start with robust PCA based on MCD (P << n)
#'      (pca1 <- PcaHubert(machines, k=3))
#'      plot(pca1, main="ROBPCA-MCD", off=0.03)
#'
#'  ## PCA with the projection algoritm of Hubert
#'      (pca2 <- PcaHubert(machines, k=3, mcd=FALSE))
#'      plot(pca2, main="ROBPCA-SD", off=0.03)
#'
#'  ## PCA with the adjusted for skewness algorithm of Hubert et al (2009)
#'      (pca3 <- PcaHubert(machines, k=3, mcd=FALSE, skew=TRUE))
#'      plot(pca3, main="ROBPCA-AO", off=0.03)
#'
#' @keywords datasets
NULL

#'
#' Skull dimensions of the wolf \emph{Canis lupus} L.
#'
#' A data set containing skull morphometric measurements on Rocky Mountain
#'  and Arctic wolves (\emph{Canis Lupus L.}). The tdata are published in Morrison (1990),
#'  originally from Jolicoeur (1959).
#'
#' @name wolves
#' @docType data
#' @usage data(wolves)
#' @format A data frame with 25 rows and 12 variables.
#' The variables are as follows (all measurements are in milimeters):
#'
#' \itemize{
#'      \item \code{class}: a factor presenting the combinations of \code{location}
#'          and \code{sex}. The levels are \code{arf} \code{arm} \code{rmf} and \code{rmm}
#'      \item \code{location}: a factor with levels \code{ar}=Arctic, \code{rm}=Rocky Mountain
#'      \item \code{sex}: a factor with levels \code{f}=female, \code{m}=male
#'      \item \code{x1}: palatal length
#'      \item \code{x2}: postpalatal length
#'      \item \code{x3}: zygomatic width
#'      \item \code{x4}: palatal width outside first upper molars
#'      \item \code{x5}: palatal width inside second upper molars
#'      \item \code{x6}: postglenoid foramina width
#'      \item \code{x7}: interorbital width
#'      \item \code{x8}: braincase width
#'      \item \code{x9}: crown length
#' }
#'
#' @source
#'  Jolicoeur, P. Multivariate geographical variation in the wolf \emph{Canis lupis L.},
#'  \emph{Evolution}, XIII, 283--299.
#'
#'  Morrison, D. F.  \emph{Multivariate Statistical Methods},  (3rd ed.), 1990.
#'  New York: McGraw-Hill, p. 288--289.
#'
#' @examples
#'
#'  data(wolves)
#'
#'  ## Remove the factors location and sex which we will not use for now
#'  x <- wolves[,-c(2:3)]
#'
#'  ## Plot a pairwise scaterplot matrix
#'  pairs(x[,2:10])
#'
#'  mcd <- CovMcd(x[, 2:10])
#'  plot(mcd, which="pairs")
#'
#'  lda <- LdaClassic(class~., data=x)
#'  lda@center
#'  lda@cov
#'
#'  predict(lda)
#'
#' @keywords datasets
NULL

#'
#' Fruit data set
#'
#' A data set that contains the spectra of six different cultivars of
#' the same fruit (cantaloupe - \emph{Cucumis melo} L. Cantaloupensis
#' group) obtained from Colin Greensill (Faculty of Engineering and Physical Systems, Central Queensland
#' University, Rockhampton, Australia). The total data set contained 2818 spectra measured in 256 wavelengths.
#' For illustrative purposes are considered only three cultivars out of it, named D, M and
#' HA with sizes 490, 106 and 500, respectively. Thus the data set thus contains 1096 observations.
#' For more details about this data set see the references below.
#' @name fruit
#' @docType data
#' @usage data(fruit)
#' @format A data frame with 1096 rows and 257 variables (one grouping variable -- \code{cultivar} -- and 256 measurement variables).
#' @source
#' Colin Greensill (Faculty of Engineering and Physical Systems, Central Queensland
#' University, Rockhampton, Australia).
#'
#' @references
#'  Hubert, M. and Van Driessen, K., (2004). Fast and robust discriminant analysis.
#'  \emph{Computational Statistics and Data Analysis}, \bold{45}(2):301--320.
#'  \doi{10.1016/S0167-9473(02)00299-2}.
#'
#'  Vanden Branden, K and Hubert, M, (2005). Robust classification in high dimensions based on the SIMCA Method.
#'  \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{79}(1-2), pp. 10--21.
#'  \doi{10.1016/j.chemolab.2005.03.002}.
#'
#'  Hubert, M, Rousseeuw, PJ and Verdonck, T, (2012). A Deterministic Algorithm for Robust Location and Scatter.
#'  \emph{Journal of Computational and Graphical Statistics}, \bold{21}(3), pp 618--637.
#'  \doi{10.1080/10618600.2012.672100}.
#'
#' @examples
#'
#'  data(fruit)
#'  table(fruit$cultivar)
#'
#' @keywords datasets
NULL
#' Johns Hopkins University Ionosphere database.
#'
#' ''This radar data was collected by a system in Goose Bay, Labrador.  This
#'   system consists of a phased array of 16 high-frequency antennas with a
#'   total transmitted power on the order of 6.4 kilowatts.  The targets
#'   were free electrons in the ionosphere.
#'   "good" radar returns are those showing evidence of some type of structure
#'   in the ionosphere.  "bad" returns are those that do not; their signals pass
#'   through the ionosphere.
#'   Received signals were processed using an autocorrelation function whose
#'   arguments are the time of a pulse and the pulse number.  There were 17
#"   pulse numbers for the Goose Bay system.  Instances in this databse are
#'   described by 2 attributes per pulse number, corresponding to the complex
#'   values returned by the function resulting from the complex electromagnetic
#'   signal.'' [UCI archive]
#'
#' @name ionosphere
#' @docType data
#' @usage data(ionosphere)
#' @format A data frame with 351 rows and 33 variables: 32 measurements and one
#'  (the last, \code{Class}) grouping variable: 225 \code{'good'} and 126 \code{'bad'}.
#'
#'  The original dataset at UCI contains 351 rows and 35 columns. The first 34
#'  columns are features, the last column contains the classification label of
#'  'g' and 'b'. The first feature is binary and the second one is only 0s,
#;  therefore these two features were removed. We remain with 32 featres and
#'  one grouping variable - factor with labels 'good' and 'bad'.
#'
#' @source
#'  Source: Space Physics Group; Applied Physics Laboratory; Johns Hopkins University; Johns Hopkins Road; Laurel; MD 20723
#'
#'  Donor: Vince Sigillito (vgs@aplcen.apl.jhu.edu)
#'
#'  The data have been taken from the UCI Repository Of Machine Learning Databases at
#'  \url{https://archive.ics.uci.edu/ml/datasets/ionosphere}
#'
#'  This data set, with the original 34 features is available in the package \pkg{mlbench}
#'  and a different data set (refering to the same UCI repository) is available in
#'  the package \code{dprep} (archived on CRAN).
#' @references
#'  Sigillito, V. G., Wing, S. P., Hutton, L. V., and Baker, K. B. (1989).
#'      Classification of radar returns from the ionosphere using neural
#'      networks. Johns Hopkins APL Technical Digest, 10, 262-266.
#' @examples
#'  data(ionosphere)
#'  pairs(ionosphere[, 1:6])
NULL
