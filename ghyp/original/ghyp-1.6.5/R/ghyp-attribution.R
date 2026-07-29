### <======================================================================>
#' Risk attribution.
#' 
#'  Functions to get the \emph{contribution} of each asset to the portfolio's 
#' \emph{Expected Shortfall} based on multivariate generalized hyperbolic 
#' distributions as well as the expected shortfall \emph{sensitivity} to marginal
#' changes in portfolio allocation.
#' 
#' 
#' @param alpha a vector of confidence levels for ES.
#' @param object a multivariate fitted ghyp object inheriting from class \code{\link[=ghyp-class]{ghyp}}.
#' @param distr whether the ghyp-object specifies a return or a loss-distribution (see \bold{Details}).
#' @param weights vector of portfolio weights. Default is an equally-weighted portfolio.
#' @param ... optional arguments passed from \linkS4class{ghyp.attribution} to \code{\link{qghyp}} and \code{\link{integrate}}.
#' 
#' @details The parameter \code{distr} specifies whether the ghyp-object 
#' describes a return or a loss-distribution.  In case of a return
#' distribution the expected-shortfall on a confidence level
#' \eqn{\alpha}{alpha} is defined as \eqn{\hbox{ES}_\alpha := 
#' \hbox{E}(X| X \leq F^{-1}_X(\alpha))}{ES_alpha := E(X | X <= F^-1(alpha))}
#' while in case of a loss distribution it is defined on a confidence
#' level \eqn{\alpha}{alpha} as \eqn{\hbox{ES}_\alpha := \hbox{E}(X | X
#' > F^{-1}_X(\alpha))}{ES_alpha := E(X | X > F^-1(alpha))}.\cr  
#'
#'  
#' @return \code{ESghyp.attribution} is an object of class \linkS4class{ghyp.attribution}.
#'
#'
#' @author Marc Weibel
#'
#' @seealso \code{\link{contribution,ghyp.attribution-method}}, \code{\link{sensitivity,ghyp.attribution-method}} and \code{\link[=weights]{weights}} for Expected Shortfall.
#'   
#'
#' @examples
#' \dontrun{
#' data(smi.stocks)
#' ## Fit a NIG model to Novartis, CS and Nestle log-returns
#' assets.fit <- fit.NIGmv(smi.stocks[, c("Novartis", "CS", "Nestle")], silent = TRUE)
#' ## Define Weights of the Portfolio
#' weights <- c(0.2, 0.5, 0.3)
#' ## Confidence level for Expected Shortfall
#' es.levels <- c(0.01)
#' 
#' portfolio.attrib <- ESghyp.attribution(alpha=es.levels, object=assets.fit, weights=weights)
#' }
#' @keywords risk attribution
#' @export
"ESghyp.attribution" <- function(alpha, object = ghyp(), distr = c("return", "loss"), weights = NULL, ...)
{
  distr <- match.arg(distr)
  
  .test.ghyp(object, case = "multivariate")
  
  ## Check if vector weights' length == dim 
  dim.object <- ghyp.dim(object)
  
  ## If weights=NULL, use an equally-weighted Portfolio
  if(is.null(weights)) 
  {
    weights <- rep(1/dim.object, dim.object)
    print('no weights have been assigned : use equally-weighted Portfolio')
  }	
  
  ## Check if vector weights corresponds to object dimension
  if(length(as.vector(weights)) != dim.object)
    stop(paste("Weights must be a vector of Dimension '", dim.object, "' !"), sep="")
  
  alpha <- .check.data(alpha, na.rm = FALSE, fit = FALSE, dim = 1)	
  ES.contrib <- ES.sensi <- matrix(0, dim.object, length(alpha))
  
  ## Check if Data are near 0, as contribution may
  ## not be computed. Rescale Parameters if necessary
  scale.factor <- 1
  if(sum(coef(object)$sigma<0.01)!=0)
  {
    scale.factor <- 100
    attributes(object)$mu <- coef(object)$mu * scale.factor
    attributes(object)$sigma <- coef(object)$sigma * scale.factor^2
    attributes(object)$gamma <- coef(object)$gamma * scale.factor
  }
  
  ##  ATTRIBUTION 
  
  ## Portfolio : object %*% weights
  portfolio.object <- transform(object, multiplier=weights)
  ES <- ESghyp(alpha, portfolio.object)
  
  if(.is.gaussian(object)){
    for (j in 1:length(alpha)){
      if(distr == "return"){					
        ES.contrib[ ,j] <-  (mean(object) - 1/sqrt(vcov(portfolio.object)) * 
                               as.vector(vcov(object)%*%weights) * 
                               dnorm(qnorm(1 - alpha[j])) / (alpha[j])) * weights		 						
      }else{                          # For losses
        ES.contrib[ ,j] <- (mean(object) + 1/sqrt(vcov(portfolio.object)) * 
                              as.vector(vcov(object)%*%weights) * 
                              dnorm(qnorm(alpha[j])) / (1-alpha[j])) * weights	
      }
    }		
  }else if(.is.student.t(object, symmetric = TRUE) & object@parametrization == "alpha.bar"){
    nu <- coef(object)$nu
    sigma.t <- sqrt((nu - 2) / nu) * object@sigma
    sigma.t.port <- portfolio.object@sigma
    for (j in 1:length(alpha)){
      if(distr == "return"){								
        ES.contrib[ ,j] <-  (object@mu - (1/sigma.t.port) * as.vector(sigma.t%*%weights) * dt(qt(1-alpha[j], df = nu), df = nu) /
                               (alpha[j]) * (nu + qt(1-alpha[j], df = nu)^2) / (nu - 1)) * weights		
      }else{				
        ES.contrib[ ,j] <-  (object@mu + 1/sigma.t.port * as.vector(sigma.t%*%weights) * dt(qt(alpha[j], df = nu), df = nu) /
                               (1 - alpha[j]) * (nu + qt(alpha[j], df = nu)^2) / (nu - 1)) * weights		
      }
    }		
  }else{
    value.raw <- qghyp(alpha, portfolio.object, ...)
    
    if(all(is.na(value.raw))){
      return(value.raw)
    }
    
    value.var <- matrix(value.raw[!is.na(value.raw)], ncol = 1)
    
    ## arguments for function 'integrate'
    ctrl <- list(subdivisions=100,				  		# the maximum number of subintervals.
                 rel.tol=1e-03, 						# relative accuracy requested.
                 abs.tol=1e-03,						# absolute accuracy requested.
                 stop.on.error=TRUE,
                 keep.xy=FALSE,
                 aux=NULL)
    namesctrl <- names(ctrl)				
    ctrl[(namesArgs <- names(list(...)))] <- list(...)				
    if (length(noNames <- namesArgs[!namesArgs %in% namesctrl])) 
      warning("unknown names in optional parameters: ", paste(noNames, collapse = ", "))		
    
    ## ES Contribution for Each Asset
    for (j in 1:length(alpha))
    {
      for (i in 1:dim.object)
      {	
        B <- rbind(rep(0,dim.object), weights)
        B[1,i] <- weights[i]
        object.attrib <- transform(object, multiplier=B)				
        
        if(distr == "return"){
          
          tmp.contrib <- try(1/alpha[j] * integrate(function(y) 
          { 
            sapply(y, function(y) 
            {
              integrate(function(x) {
                sapply(x, function(x) dghyp(c(y,x),object.attrib)*y)
              }, lower=-Inf, upper=value.var[j], 
              subdivisions=ctrl$subdivisions, 
              rel.tol=ctrl$rel.tol, 
              abs.tol=ctrl$abs.tol,
              stop.on.error=ctrl$stop.on.error,
              keep.xy=ctrl$keep.xy,
              aux=NULL)$value
            })
          }, lower=-Inf, upper=Inf, subdivisions=ctrl$subdivisions, 
          rel.tol=ctrl$rel.tol, 
          abs.tol=ctrl$abs.tol,
          stop.on.error=ctrl$stop.on.error,
          keep.xy=ctrl$keep.xy,
          aux=NULL)$value, silent=TRUE)												
        }else{
          tmp.contrib <- try(1/(1-alpha[j]) * integrate(function(y) 
          { 
            sapply(y, function(y) 
            {
              integrate(function(x) {
                sapply(x, function(x) dghyp(c(y,x),object.attrib)*y)
              }, lower=value.var[j], upper=Inf, 
              subdivisions=ctrl$subdivisions, 
              rel.tol=ctrl$rel.tol, 
              abs.tol=ctrl$abs.tol,
              stop.on.error=ctrl$stop.on.error,
              keep.xy=ctrl$keep.xy,
              aux=NULL)$value
            })
          }, lower=-Inf, upper=Inf, subdivisions=ctrl$subdivisions, 
          rel.tol=ctrl$rel.tol, 
          abs.tol=ctrl$abs.tol,
          stop.on.error=ctrl$stop.on.error,
          keep.xy=ctrl$keep.xy,
          aux=NULL)$value, silent=TRUE)	
        }			
        
        if(is(tmp.contrib, "try-error")){
          stop("Failed to determine contribution \n")
          ES.contrib[i, j] <- NA
        } else {
          ES.contrib[i, j] <- tmp.contrib
        }	
      }
    }        
  }
  
  ## Scale data back to original 
  ## and scale contribution to overall ES	
  ES.contrib <- t(t(ES.contrib) * ES/colSums(ES.contrib))
  ES.contrib <- ES.contrib / scale.factor
  
  ## Compute Sensitivity of each Asset
  ES.sensi <- ES.contrib / weights
  
  ## Name Columns and Rows 
  colnames(ES.contrib) <- colnames(ES.sensi) <- paste(alpha*100, "%", sep="")
  object.names <- rownames(coef(object)$sigma)
  if(!is.null(object.names)) 
    rownames(ES.contrib) <- rownames(ES.sensi) <- names(weights) <- object.names	
  
  # Name Weights and ES Vectors
  weights <- as.matrix(weights)		
  ES <- t(as.matrix(colSums(as.matrix(ES.contrib))))
  colnames(weights) <- rownames(ES) <- 'Portfolio'
  
  # Returns ES, Contribution and Sensitivity
  return(new("ghyp.attribution", weights=weights, ES = ES, 
             contribution=ES.contrib, sensitivity=ES.sensi))
}
### <---------------------------------------------------------------------->

