
#' @title Elastic Net Penalized Exponentially Distributed Response Variables
#'
#' @description \code{git.glmGammaNet} Fit glmnet model for exponentiall distributed response data.
#'
#' @param A The matrix of independent variables.
#' @param b The vector of response variables.
#' @param alpha.EN The coefficient of elastic net regularizer (1 means lasso).
#' @param num_lambda Size of the lambda grid.
#' @param glm_type Type of glm model, 1 is exponential, 2 is gamma (not implemented yet).
#' @param max_iter Max number of iteration for the prox grad descent optimizer.
#' @param abs_tol Absolute error threshold for the pgd optimizer.
#' @param rel_tol Relative error threshold for the pgd optimizer (not used for vanilla PGD).
#' @param normalize_grad Swtich for whether to normalize the gradient or not.
#' @param k_fold The number of folds for cross validation.
#' @param k_fold_iter The number of iterations for the cross-validation.
#' @param has_intercept Parameter to determine if there is an intercept (TRUE) or not (FALSE).
#' @param ... Additional Parameters.
#'
#' @return Vector of optimal coefficient for the glm model.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @examples
#' # Function to return the periodogram of data series
#' myperiodogram <- function (data, max.freq = 0.5, 
#'                            twosided = FALSE, keep = 1){
#'   data.fft <- fft(data)
#'   N <- length(data)
#'   tmp <- Mod(data.fft[2:floor(N/2)])^2/N
#'   tmp <- sapply(tmp, function(x) max(1e-05, x))
#'   freq <- ((1:(floor(N/2) - 1))/N)
#'   tmp <- tmp[1:floor(length(tmp) * keep)]
#'   freq <- freq[1:floor(length(freq) * keep)]
#'   if (twosided) {
#'     tmp <- c(rev(tmp), tmp)
#'     freq <- c(-rev(freq), freq)
#'   }
#'   return(list(spec = tmp, freq = freq))
#' }
#' 
#' # Function to compute the standard error based the periodogram of 
#' # the influence functions time series
#' SE.Exponential <- function(data, d = 7, alpha = 0.5, keep = 1){
#'   N <- length(data)
#'   # Compute the periodograms
#'   my.periodogram <- myperiodogram(data)
#'   my.freq <- my.periodogram$freq
#'   my.periodogram <- my.periodogram$spec
#'   # Remove values of frequency 0 as it does not contain information 
#'   # about the variance
#'   my.freq <- my.freq[-1]
#'   my.periodogram <- my.periodogram[-1]
#'   # Implement cut-off
#'   nfreq <- length(my.freq)
#'   my.freq <- my.freq[1:floor(nfreq*keep)]
#'   my.periodogram <- my.periodogram[1:floor(nfreq*keep)]
#'   # GLM with BFGS optimization
#'   # Create 1, x, x^2, ..., x^d
#'   x.mat <- rep(1,length(my.freq))
#'   for(col.iter in 1:d){
#'     x.mat <- cbind(x.mat,my.freq^col.iter)
#'   }
#'   # Fit the Exponential model
#'   res <- glmnet_exp(x.mat, my.periodogram, alpha.EN = alpha)
#'   # Return the estimated variance
#'   return(sqrt(exp(res[1])/N))
#' }
#' 
#' # Examples using PerformanceAnalytics data
#' if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
#'   # Loading hedge fund data from PA
#'   data(edhec, package = "PerformanceAnalytics")
#'   colnames(edhec)
#'   
#'   # Computing the expected shortfall for the time series of returns
#'   # library(RPEIF)
#'   # test.mat <- apply(edhec, 2, IF.ES)
#'   # test.mat <- apply(test.mat, 2, as.numeric)
#'   
#'   # Returning the standard errors from the Exponential distribution fit
#'   # apply(test.mat, 2, SE.Exponential)
#' }
#' 
glmnet_exp <- function(A,
                      b,
                      alpha.EN = 0.5,
                      num_lambda = 100L,
                      glm_type = 1L,
                      max_iter = 100L,
                      abs_tol = 1.0e-4,
                      rel_tol = 1.0e-2,
                      normalize_grad = FALSE,
                      k_fold = 5L,
                      has_intercept = TRUE,
                      k_fold_iter = 5L,
                      ...){
  
  # Testing input for independent variables
  if(!is.matrix(A))
    stop("A must be a matrix.")
  # Testing input for response variable
  if(!is.vector(b))
    stop("b must be a vector of responses.")
  
  return(fitGlmCv(A,
         b,
         alpha = alpha.EN,
         num_lambda = num_lambda,
         glm_type = glm_type,
         max_iter = max_iter,
         abs_tol = abs_tol,
         rel_tol = rel_tol,
         normalize_grad = normalize_grad,
         k_fold = k_fold,
         has_intercept = has_intercept,
         k_fold_iter = k_fold_iter))
}
