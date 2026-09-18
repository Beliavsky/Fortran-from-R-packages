##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' Computes the Fowlkes and Mallows index
#' 
#' @name FowlkesMallowsIndex
#' @aliases print.FowlkesMallowsIndex
#' @description  Fowlkes-Mallows index is an external evaluation method
#'  that is used to determine the similarity between two clusterings
#'  (clusters obtained after a clustering algorithm). This measure of
#'  similarity could be either between two hierarchical clusterings or a
#'  clustering and a benchmark classification. A higher the value for the
#'  Fowlkes-Mallows index indicates a greater similarity between the clusters
#'  and the benchmark classifications.
#'  This index can be used to compare either two cluster label sets or a
#'  cluster label set with a true label set. The formula of the
#'  adjusted Fowlkes-Mallows index (ABk) is given in the details.
#'
#' @param c1 Labels of the first partition or contingency table. A numeric or 
#'  character vector containing the class labels of the first partition or 
#'  a 2-dimensional numeric matrix which contains the cross-tabulation 
#'  of cluster assignments.
#' @param c2 Labels of the second partition. A numeric or character vector 
#'  containing the class labels of the second partition. 
#'  The length of vector c2 must be equal to the length of vector c1. 
#'  This second input is required only if \code{c1} is not a 2-dimensional numeric matrix.
#'
#' @param noisecluster Label or number associated to the \emph{noise class} or \emph{noise level}.
#'  Number or character label which denotes the points which do not belong to any cluster.
#'  These points are not takern into account for the  computation of the Fowlkes and Mallows index
#'
#' @return A list containing the following components:
#'  \item{ABK}{Adjusted Fowlkes and Mallows index. A number between -1 and 1. The adjusted 
#'  Fowlkes and Mallows index is the corrected-for-chance version of the Fowlkes and Mallows index.}
#'  \item{BK}{Value of the Fowlkes and Mallows index. A number between 0 and 1. }
#'  \item{EBk}{Expected value of the index computed under the null hypothesis of no-relation.}
#'  \item{VarBk}{Variance of the Fowlkes and Mallows index. Variance of the index computed under the null
#'      hypothesis of no-relation.}
#'
#' @details The formula of the adjusted Fowlkes-Mallows index (ABk) is as follows:
#'
#'  \deqn{ABk=\frac{Bk-Expected~value~of~Bk}{Max~Index - Expected~value~of~Bk}}
#'
#' @references
#'  Fowlkes, E.B. and Mallows, C.L. (1983), A Method for Comparing Two Hierarchical Clusterings, 
#'  \emph{Journal of the American Statistical Association}, Vol. 78, pp. 553--569.
#'
#' @seealso \code{\link{randIndex}}
#' @examples
#'
#'  ##  1. FowlkesMallowsIndex (adjusted) with the two vectors as input
#'  c <- matrix(c(1, 1, 1, 2, 2, 1,2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3), 
#'          ncol=2, byrow=TRUE)
#'
#'  ##  c1 - numeric vector containing the labels of the first partition
#'  c1 <- c[, 1]
#'
#'  ##  c2 - numeric vector containing the labels of the second partition
#'  c2 <- c[, 2]
#'
#'  (FM <- FowlkesMallowsIndex(c1, c2))
#'
#'  ##  2. FM index (adjusted) with the contingency table as input.
#'  T <- matrix(c(1, 1, 0, 1, 2, 1, 0, 0, 4), ncol=3, byrow=TRUE)
#'  (FM <- FowlkesMallowsIndex(T))
#'
#'
#'  ##  3. Compare FM (unadjusted) for iris data (true classification against 
#'  ##      tclust classification).
#'
#'  ##  First partition c1 is the true partition
#'  c1 <- iris$Species
#'
#'  ##  Second partition c2 is the output of tclust clustering procedure
#'  out <- tclust(iris[, 1:4], k=3, alpha=0, restr.fact=100)
#'  c2<- out$cluster
#'
#'  (FM <- FowlkesMallowsIndex(c1, c2))
#'
#'  ##  4. Compare FM index (unadjusted) for iris data (exclude unassigned units from tclust).
#'
#'  ##  First partition c1 is the true partition
#'  c1 <- iris$Species
#'
#'  ##  Second partition c2 is the output of tclust clustering procedure
#'  out <- tclust(iris[, 1:4], k=3, alpha=0.1, restr.fact=100)
#'  c2<- out$cluster
#'
#'  ##  Units inside c2 which contain number 0 are referred to trimmed observations
#'  noisecluster <- 0
#'  (FM <- FowlkesMallowsIndex(c1, c2, noisecluster=noisecluster))
#'

FowlkesMallowsIndex <- function(c1, c2 = NULL, noisecluster = NULL) {

    if(!is.matrix(c1))
        c1 <- as.matrix(c1)
    if(!missing(c2) && !is.matrix(c2))
        c2 <- as.matrix(c2)
        
    if(missing(c2)) {
        if(ncol(c1) < 2) {
            stop("FowlkesMallowsIndex: Requires a contingency table with at least two columns")
        }
        
        M <- c1
    } else {
        if(min(dim(c1)) > 1 || min(dim(c2)) > 1) {
            stop("FowlkesMallowsIndex: Requires two vector arguments")
        }
    
        if(!is.null(noisecluster)) {
            boo1 <- if (is.character(c1) || is.factor(c1)) c1 == noisecluster else c1 == noisecluster
            boo2 <- if (is.character(c2) || is.factor(c2)) c2 == noisecluster else c2 == noisecluster
            keep <- !(boo1 | boo2)
            c1 <- c1[keep]
            c2 <- c2[keep]
        }
    
        if (is.character(c1) || is.factor(c1) || is.character(c2) || is.factor(c2)) {
            M <- table(c1, c2)
        } else {
            class_c1 <- sort(unique(c1))
            class_c2 <- sort(unique(c2))
            Rows <- length(class_c1)
            Cols <- length(class_c2)
            M <- matrix(0, nrow = Rows, ncol = Cols)
            for (i in seq_len(Rows)) {
                for (j in seq_len(Cols)) {
                    M[i, j] <- sum(c1 == class_c1[i] & c2 == class_c2[j])
                }
            }
            
        }    
    }

    midot <- rowSums(M)
    mdotj <- colSums(M)
    n <- sum(M)
    
    Tk <- sum(M^2) - n
    Pk <- sum(midot^2) - n
    Qk <- sum(mdotj^2) - n
    Bk <- Tk / sqrt(Pk * Qk)
    
    EBk <- sqrt(Pk * Qk) / (n * (n - 1))
    
    Pk2 <- sum(midot * (midot - 1) * (midot - 2))
    Qk2 <- sum(mdotj * (mdotj - 1) * (mdotj - 2))

    VarBk <- 2 / (n * (n - 1)) +
        4 * Pk2 * Qk2 / ((n * (n - 1) * (n - 2)) * Pk * Qk) +
        (Pk - 2 - 4 * Pk2 / Pk) * (Qk - 2 - 4 * Qk2 / Qk) /
        (n * (n - 1) * (n - 2) * (n - 3)) -
        Pk * Qk / (n^2 * (n - 1)^2)

    ABk <- (Bk - EBk) / (1 - EBk)

    return(list(ABk = ABk, Bk = Bk, EBk = EBk, VarBk = VarBk))
}
