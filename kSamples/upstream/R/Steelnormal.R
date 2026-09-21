Steelnormal <- function(mu,sig0,sig,tau,Wvec,ni,
		alternative=c("greater","less","two-sided"),
		continuity.corr=TRUE,corr){
#
# This function evalutates the normal approximation significance probability
# of the Steel test statistic for any of the Mann-Whitney statistics in k comparisons 
# for each of k treatment samples Y with a common control sample X, based on the 
# mutiple comparison Steel test distribution,thus giving the adjusted p-values for 
# each of the Mann-Whitney statistics.
# The minimum of these adjusted p-values then serves as the overall p-value of the 
# Steel test.
# Wvec contains the k Mann-Whitney statistics, 
# mu,sig0,sig,tau contain parameters required for the normal approximation and 
# ni contains the sizes of the treatment samples.
#
    alternative <- match.arg(alternative)
    k <- length(ni)
    pval <- numeric(k)
    if(continuity.corr == FALSE) corr <- corr*0
    if (alternative == "greater") {
        for (j in 1:k) {
            S <- (Wvec[j] -corr[j] - mu[j])/tau[j]
            funz <- function(z, k, sig0, sig, tau, S, ni) {
                fac <- 1
                for (i in 1:k) {
                  fac <- fac * pnorm((S * tau[i] - ni[i] * sig0 *
                    z)/sig[i])
                }
                dnorm(z) * fac
            }
            pval[j] <- 1 - integrate(funz, -10, 10, k, sig0,
                sig, tau, S, ni)$value
        }
    }
    if (alternative == "less") {
        for (j in 1:k) {
            S <- (Wvec[j] + corr[j]- mu[j])/tau[j]
            funz <- function(z, k, sig0, sig, tau, S, ni) {
                fac <- 1
                for (i in 1:k) {
                  fac <- fac * (1 - pnorm((S * tau[i] - ni[i] *
                    sig0 * z)/sig[i]))
                }
                dnorm(z) * fac
            }
            pval[j] <- 1 - integrate(funz, -10, 10, k, sig0,
                sig, tau, S, ni)$value
        }
    }
    if (alternative == "two-sided") {
        for (j in 1:k) {
            S <- (abs(Wvec[j] -corr[j] - mu[j]))/tau[j]
            funz <- function(z, k, sig0, sig, tau, S, ni) {
                fac <- 1
                for (i in 1:k) {
                  fac <- fac * (pnorm((S * tau[i] - ni[i] * sig0 *
                    z)/sig[i]) - pnorm((-S * tau[i] - ni[i] *
                    sig0 * z)/sig[i]))
                }
                dnorm(z) * fac
            }
            pval[j] <- 1 - integrate(funz, -10, 10, k, sig0,
                sig, tau, S, ni)$value
        }
    }
    pval
}
# end of Steelnormal
