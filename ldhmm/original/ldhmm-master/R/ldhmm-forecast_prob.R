#' Computing the forecast probability distribution
#' 
#' This utility computes the forecast probability distribution (Zucchini, 5.3)
#'
#' @param object an ldhmm object
#' @param x numeric, the observations.
#' @param xf numeric, the future observations to be forecasted.
#' @param h integer, time steps to forecast.
#'
#' @return matrix of probabilities, size of h times size of xf.
#'
#' @keywords forecast
#'
#' @author Stephen H. Lihn
#'
#' @export
#' 
### <======================================================================>
ldhmm.forecast_prob <- function(object, x, xf, h=1)
{
    m <- object@m
    n <- length(x)
    nxf <- length(xf)
    
    pdf <- ldhmm.state_pdf(object, 1:m, x)
    if (is(pdf, "numeric")) pdf <- as.matrix(pdf, nrow=m, ncol=n)
    
    phi <- object@delta * pdf[,1]
    sum_phi <- sum(phi)
    log_phi_scale <- log(sum_phi)
    phi <- phi/sum_phi
    for (i in 2:n) {
        phi <- phi %*% object@gamma * pdf[,i]
        sum_phi <- sum(phi)
        log_phi_scale <- log_phi_scale + log(sum_phi)
        phi <- phi/sum_phi
    }
    # at this point, phi is phi_T = alpha_T/|alpha_T|
    phi_T <- phi
    # TODO separate this to another function below
    
    pdf_f <- ldhmm.state_pdf(object, 1:m, xf)
    if (is(pdf_f, "numeric")) pdf_f <- as.matrix(pdf_f, nrow=m, ncol=nxf)
    df <- matrix(0, nrow=h, ncol=nxf)
    for (i in 1:h) {
        phi <- phi %*% object@gamma
        for (j in 1:m) df[i,] <- df[i,] + phi[j] * pdf_f[j,]
    }
    return(df)
}
### <---------------------------------------------------------------------->
