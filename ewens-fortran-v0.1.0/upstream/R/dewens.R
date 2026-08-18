#' Probability mass function for the Ewens distribution
#'
#' @description
#' Gives the probability mass function for the Ewens distribution, as described in Ewens, Warren (1972). "The sampling theory of selectively neutral alleles". *Theoretical Population Biology*. 3: 87–112. doi:10.1016/0040-5809(72)90035-4.
#'
#' @details
#' The probability of a vector of counts \eqn{m_1, ..., m_n} is given by the expression \deqn{\frac{n!}{\theta (\theta + 1) ... (\theta + n - 1)}\prod_{j=1}^n \frac{\theta^{m_j}}{j^{m_j} m_j!}}
#' 
#' @param x A vector giving class memberships of each observation in the sample
#' @param theta A non-negative parameter governing the expected sample diversity.
#' @param log if TRUE, probabilities are given as log(p). Default is FALSE.
#' @returns A numeric vector giving a probability (or if `log = TRUE`, a log probability)
#' @examples
#' x <- sample(LETTERS, 120, replace = TRUE)
#' dewens(x, theta = 1)
#' dewens(x, theta = 0) ## returns NaN since vector incompatible with zero diversity
#'
#' @export
dewens <- function(x, theta = 1, log = FALSE) {
    if (theta < 0) {
        stop("Parameter theta must be non-negative. ")
    }
    ### Get counts
    tab <- tabulate(factor(x))
    mvec   <- tabulate(tab, nbins = max(tab))
    n      <- length(x)
    j_seq  <- seq_along(mvec)

    log_pmf <- lfactorial(n) -
        sum(log(theta + seq_len(n) - 1)) +   # log rising factorial
        sum(mvec * log(theta)                  # theta^{a_j}
            - mvec * log(j_seq)               # j^{a_j}
            - lfactorial(mvec))               # a_j!

    if (log) {
        return(log_pmf)
    } else {
        return(exp(log_pmf))
    }
}
