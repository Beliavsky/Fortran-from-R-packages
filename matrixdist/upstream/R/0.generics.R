#' New generic for simulating matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return A realization from the matrix distribution.
#'
setGeneric("sim", function(x, ...) standardGeneric("sim"))

#' New generic for moments of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Moment of the matrix distribution.
#'
setGeneric("moment", function(x, ...) standardGeneric("moment"))

#' New generic for Laplace transform of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Laplace transform of the matrix distribution.
#'
setGeneric("laplace", function(x, ...) standardGeneric("laplace"))

#' New generic for mgf of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Mgf of the matrix distribution.
#'
setGeneric("mgf", function(x, ...) standardGeneric("mgf"))

#' New generic for pgf of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{dph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Pgf of the matrix distribution.
#'
setGeneric("pgf", function(x, ...) standardGeneric("pgf"))

#' New generic for minimum of two matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("minimum", function(x1, x2, ...) standardGeneric("minimum"))

#' New generic for maximum of two matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("maximum", function(x1, x2, ...) standardGeneric("maximum"))

#' New generic for mixture of two matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{ph} and \linkS4class{dph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("mixture", function(x1, x2, ...) standardGeneric("mixture"))

#' New generic for N-fold convolution of two matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{ph} and \linkS4class{dph}.
#'
#' @param x1 An object of the class \linkS4class{dph}.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("Nfold", function(x1, x2, ...) standardGeneric("Nfold"))

#' New generic for the density of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Density from the matrix distribution.
#'
setGeneric("dens", function(x, ...) standardGeneric("dens"))

#' New generic for the hazard rate of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Hazard rate from the matrix distribution.
#'
setGeneric("haz", function(x, ...) standardGeneric("haz"))

#' New generic for the distribution of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return CDF from the matrix distribution.
#'
setGeneric("cdf", function(x, ...) standardGeneric("cdf"))

#' New generic for the quantile of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Quantile from the matrix distribution.
#'
setGeneric("quan", function(x, ...) standardGeneric("quan"))

#' New generic for estimating matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#'
setGeneric("fit", function(x, y, ...) standardGeneric("fit"))

#' New generic for regression with matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#' @export
#'
setGeneric("reg", function(x, y, ...) standardGeneric("reg"))

#' New generic for evaluating survival matrix distributions
#'
#' Methods are available for objects of class \linkS4class{sph}.
#'
#' @param x An object of the model class.
#' @param subject A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @export
#'
setGeneric("evaluate", function(x, subject, ...) standardGeneric("evaluate"))

#' New generic for likelihood ratio test between two matrix distribution models
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x,y Objects of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return A likelihood ratio test result.
#' @export
#'
setGeneric("LRT", function(x, y, ...) standardGeneric("LRT"))

#' New generic for obtaining the Fisher information of survival matrix distributions
#'
#' Methods are available for objects of class \linkS4class{sph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @export
#'
setGeneric("Fisher", function(x, ...) standardGeneric("Fisher"))

#' New generic for transformation via rewards of a matrix distribution
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("TVR", function(x, ...) standardGeneric("TVR"))

#' New generic for mixture-of-experts regression with matrix distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#' @export
#'
setGeneric("MoE", function(x, y, ...) standardGeneric("MoE"))

#' New generic for the marginals of multivariate matrix distributions
#'
#' Methods are available for objects of multivariate classes.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Marginal of the matrix distribution.
#'
setGeneric("marginal", function(x, ...) standardGeneric("marginal"))

#' New generic for linear combinations of multivariate matrix distributions
#'
#' Methods are available for objects of multivariate classes.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Marginal of the matrix distribution.
#'
setGeneric("linCom", function(x, ...) standardGeneric("linCom"))

setGeneric("mean")

setGeneric("var")

setGeneric("cor")
