#' Computing the statistics for each state
#' 
#' This utility computes the statistics (mean, sd, kurtosis, length) for each state.
#' It can be based on the local or global decoding result. The concept of asymptotic 
#' statistics can be applied by which the largest N observations (in absolute term)
#' can be dropped to avoid distortion from outliers. 
#' It is assumed the object already has come with filled data in 
#' \code{observations, states.prob, states.local, states.global} slots.
#'
#' @param object an ldhmm object that contains the observations.
#' @param drop numeric, an integer to drop the largest N observations, default is zero.
#' @param use.local logical, use local decoding result, default is \code{TURE}. 
#'                  Otherwise, use global decoding result.
#' @param x numeric, the observations.
#'
#' @return an ldhmm object containing results of decoding, based on data
#'
#' @keywords mllk
#'
#' @author Stephen H. Lihn
#'
#' @export ldhmm.calc_stats_from_obs
#' @export ldhmm.drop_outliers
#' 
#' @importFrom utils head
#' @importFrom stats sd
#' @importFrom moments skewness
#' @importFrom moments kurtosis
#'
### <======================================================================>
ldhmm.calc_stats_from_obs <- function(object, drop=0, use.local=TRUE)
{
    m <- object@m
    x <- object@observations

    calc_stats <- function(states) {
        s <- matrix(NA, nrow=m, ncol=6, byrow=TRUE)
        for (j in 1:m) {
            y <- x[states==j]
            if (drop > 0) y <- ldhmm.drop_outliers(y, drop)
            s[j,1] <- mean(y)
            s[j,2] <- stats::sd(y)
            s[j,3] <- moments::kurtosis(y)
            s[j,4] <- moments::skewness(y)
            s[j,5] <- length(y)
            s[j,6] <- s[j,1]/s[j,2]^2
        }
        colnames(s) <- c("mean", "sd", "kurtosis", "skewness", "length", "alloc")
        return(s)
    }
    
    if (use.local) {
        return(calc_stats(object@states.local))
    } else {
        return(calc_stats(object@states.global))
    }
}
### <---------------------------------------------------------------------->
#' @rdname ldhmm.calc_stats_from_obs
ldhmm.drop_outliers <- function(x, drop=1)
{
    if (drop==0) return(x)
    if (drop<0) stop("Drop must be a positive integer")
    
    `%notin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L
    z <- abs(x)
    z1 <- rev(z[order(z)])
    x[abs(x) %notin% utils::head(z1, drop)]
}
### <---------------------------------------------------------------------->
