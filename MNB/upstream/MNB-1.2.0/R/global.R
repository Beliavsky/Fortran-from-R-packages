#' Global influence
#'
#' @description  It performers influence analysis by a global influence to evaluate the impact on the parameter
#' estimates when we remove a particular observation.
#' @param star Initial values for the parameters to be optimized over.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @param plot TRUE or FALSE. Indicates if a graph should be plotted.
#' @details The function returns a list (L) with the generalized Cook distance, Likelihood displacement and
#' index plot.
#' @return L and graphics
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' }
#' @examples
#'
#' \donttest{
#'
#' data(seizures)
#' head(seizures)
#'
#' star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)
#' global.MNB(formula=Y ~ trt + period +
#' trt:period + offset(log(weeks)),star=star,dataSet=seizures,plot=FALSE)
#'
#' }
#'
#' @export
#' @import graphics

global.MNB <- function(formula,star,dataSet,plot=TRUE) {

    dataSet.orig <- dataSet

    op.glob <- fit.MNB(star = star, formula = formula,
        dataSet = dataSet,tab=FALSE)
    vp <- op.glob$par
    dfunc <- op.glob$value
    mhess <- op.glob$hessian

    dataSet.ind <- split(dataSet, f = dataSet$ind)
    n <- length(dataSet.ind)

    cook <- numeric()
    like <- numeric()
    for (i in 1:n) {
        ii <- i
        newdataSet <- subset(dataSet, dataSet$ind !=
            ii)
        op.del <- fit.MNB(star = star, formula = formula,
            dataSet = newdataSet,tab=FALSE)
        vp.del <- op.del$par
        dfunc.del <- op.del$value
        mhess.del <- op.del$hessian
        cook[i] <- t(vp.del - vp) %*% solve(-mhess) %*%
            (vp.del - vp)
        like[i] <- 2 * (dfunc - dfunc.del)
        dataSet <- dataSet.orig
    }

    result <- list()
    result$genCookDistance <- cook
    result$likeDisplacement <- like

    if (plot==TRUE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
      graphics::plot(cook, xlab = "Index", ylab = "Generalized Cook Distance",
            pch = 15, main = "", cex.axis = 1.2,
            cex.lab = 1.2, cex = 0.6, bg = 5)

      graphics::plot(abs(like), xlab = "Index", ylab = "|Likelihood Displacement|",
            pch = 15, main = "", cex.axis = 1.2,
            cex.lab = 1.2, cex = 0.6, bg = 5)

    }

    return(result)
}



