### <======================================================================>
setClass("ghyp",
         representation(call = "call",
                        lambda = "numeric",
                        alpha.bar = "numeric",
                        chi = "numeric",
                        psi = "numeric",
                        mu = "numeric",
                        sigma = "matrix",
                        gamma = "numeric",
                        model = "character",
                        dimension = "numeric",
                        expected.value = "numeric",
                        variance = "matrix",
                        parametrization = "character",
                        data = "matrix"),
          prototype(call = call("ghyp"),
                    lambda = 1,
                    alpha.bar = 1,
                    chi = 1,
                    psi = 1,
                    mu = 0,
                    sigma = matrix(0),
                    gamma = 0,
                    model = "Symmetric Generalized Hyperbolic",
                    dimension = 1,
                    expected.value = numeric(0),
                    variance = matrix(0),
                    parametrization = "alpha.bar",
                    data = matrix(0))
)
### <---------------------------------------------------------------------->

### <======================================================================>
setClass("mle.ghyp",
         representation(n.iter = "numeric",
                        llh = "numeric",
                        converged = "logical",
                        error.code = "numeric",
                        error.message = "character",
                        fitted.params = "logical",
                        aic = "numeric",
                        parameter.variance = "matrix",
                        trace.pars = "list"),

         prototype(n.iter = numeric(0),
                   llh = numeric(0),
                   converged = FALSE,
                   error.code = 0,
                   error.message = character(0),
                   fitted.params = logical(0),
                   aic = numeric(0),
                   parameter.variance = matrix(0),
                   trace.pars = list()),

         contains = "ghyp"
)
### <---------------------------------------------------------------------->

### <======================================================================>
#' Class ghyp.attribution
#'
#'  The class \dQuote{ghyp.attribution} contains the Expected Shortfall of 
#'  the portfolio as well as the contribution of each asset to the total risk
#'  and the sensitivity of each Asset. The sensitivity gives an information 
#'  about the overall risk modification of the portfolio if the weight in a 
#'  given asset is marginally increased or decreased (1 percent).
#'  
#' @docType class
#' @section Objects from the Class:
#' Objects should only be created by calls to the constructors \code{\link{ESghyp.attribution}}.
#' 
#' @slot ES Portfolio's expected shortfall (ES) for a given confidence level. Class \code{matrix}.
#' @slot contribution Contribution of each asset to the overall ES. Class \code{matrix}.
#' @slot sensitivity Sensitivity of each asset. Class \code{matrix}.
#' @slot weights Weight of each asset.
#' 
#' @method \dQuote{plot} plot
#' @method \dQuote{weights} weights
#' @method \dQuote{contribution} contribution
#' @method \dQuote{sensitivity} sensitivity
#' 
#' @author Marc Weibel
#' 
#' @note When showing special cases of the generalized hyperbolic distribution
#'the corresponding fixed parameters are not printed.
#'
#' @examples 
#' \dontrun{
#' data(smi.stocks)
#' multivariate.fit <- fit.ghypmv(data = smi.stocks,
#' opt.pars = c(lambda = FALSE, alpha.bar = FALSE),
#' lambda = 2)
#' 
#' portfolio <- ESghyp.attribution(0.01, multivariate.fit)
#' summary(portfolio)
#' }
#' @keywords classes
#' @name ghyp.attribution-class
#' @rdname ghyp.attribution-class
#' @exportClass 
setClass("ghyp.attribution",
         representation(weights = 'matrix',
                        ES = 'matrix',                                               
                        contribution = 'matrix',                                                                      
                        sensitivity = 'matrix'),
         prototype(weights = matrix(0),
                   ES = matrix(0),                    
                   contribution = matrix(0),                    
                   sensitivity = matrix(0))		 	
)
### <---------------------------------------------------------------------->
