#' Computing the MLEs
#' 
#' Computing the MLEs using \code{nlm} package
#'
#' @param object an ldhmm object that can supply m, param.nbr and stationary.
#' @param x numeric, the observations.
#' @param min.gamma numeric, a minimum transition probability added to gamma to avoid singularity, default is \code{1e-6}.
#' @param decode logical, run decoding after optimization, default is \code{FALSE}.
#' @param plot.fn name of the function that takes ldhmm object. It will be called occasionally to track the progress
#'                of the fit, mainly by plotting the time series and states. E.g. When one fits the SPX index,
#'                the function \code{ldhmm.oxford_man_plot_obs} can be used to show the expected volatility vs 
#'                Oxford-Man realized volatility. Default is \code{NULL}.
#' @param plot.interval a positive integer, specifying how often to invoke plot function, default is 200 iterations.
#' @param ssm.fn name of the function that takes ldhmm object. This function is called after the MLLK call.
#'               The purpose is to generate an additional score for optimization. E.g. It can be used to separate
#'               the states into predefined intervals, modeling a state space model. Default is \code{NULL}.
#' @param print.level numeric, this argument determines the level of printing 
#'                    which is done during the minimization process. 
#'                    The default value of 0 means that no printing occurs, 
#'                    a value of 1 means that initial and final details are printed 
#'                    and a value of 2 means that full tracing information is printed.
#' @param iterlim numeric, a positive integer specifying the maximum number of iterations 
#'                to be performed before the program is terminated.
#' @param ... additional parameters passed to the MLE optimizer 
#'
#' @return an ldhmm object containg results of MLE optimization
#'
#' @keywords mle
#'
#' @author Stephen H. Lihn
#' 
#' @importFrom stats nlm
#' @importFrom optimx optimx
#'
#' @export
#' 
#' @examples
#' \dontrun{
#'     param0 <- matrix(c(0.003, 0.02, 1, -0.006, 0.03, 1.3), 2, 3, byrow=TRUE)
#'     gamma0 <- ldhmm.gamma_init(m=2, prob=c(0.9, 0.1, 0.1, 0.9))
#'     h <- ldhmm(m=2, param=param0, gamma=gamma0)
#'     spx <- ldhmm.ts_log_rtn()
#'     ldhmm.mle(h, spx$x)
#' }
### <======================================================================>
ldhmm.mle <- function(object, x, min.gamma=1e-6, decode=FALSE, 
                      plot.fn=NULL, plot.interval=200, ssm.fn=NULL, 
                      print.level=0, iterlim=1000, ...)
{
    m  <- object@m
    object@gamma <- ldhmm.gamma_init(m, prob=object@gamma, min.gamma=min.gamma) # normalize gamma
    x <- as.numeric(x)
    
    v <- ldhmm.n2w(object)
    np  <- length(v)
    n   <- sum(!is.na(x))
    
    mu.scale <- mean(abs(x), na.rm=TRUE)
    iter_cnt <- 0
    # prepare function for optimizer 
    fn <- function(v) {
        iter_cnt <<- iter_cnt + 1
        object <- ldhmm.w2n(object, v, mu.scale=mu.scale)
        object@iterations <- iter_cnt 
        N <- length(v)
        if (!is.null(plot.fn)) {
            if (iter_cnt %% plot.interval == 0 & iter_cnt > N*2) plot.fn(object) 
        }
        mllk <- ldhmm.mllk(object, x=x, 
                           mllk.print.level=print.level)
        if (is.null(ssm.fn)) return (mllk)
        return (mllk + ssm.fn(object))
    }
    # call optimizer 
    if (object@mle.optimizer == "nlm") 
    {
        v <- ldhmm.n2w(object)
        mod <- stats::nlm(fn, p=v,
                          print.level=print.level, iterlim=iterlim)

        # construct result object
        p <- ldhmm.w2n(object, mod$estimate, mu.scale=mu.scale)
        p@mllk <- mod$minimum
        p@return.code <- mod$code
        p@iterations <- mod$iterations 
    } 
    else if (object@mle.optimizer == "BFGS") 
    {
        v <- ldhmm.n2w(object)
        opt.out <- optimx::optimx(v, fn, 
                                  method = c("BFGS"),
                                  itnmax = iterlim, 
                                  control = list(trace=TRUE))
        
        # construct result object
        v.out <- as.numeric(opt.out[1,])[1:np]
        p <- ldhmm.w2n(object, v.out, mu.scale=mu.scale)
        p@mllk <- opt.out$value
        p@return.code <- opt.out$convcode
        p@iterations <- iter_cnt # opt.out$niter 
        
    } 
    else stop(paste("Error: Unknown MLE optimizer", object@mle.optimizer))
    
    p@AIC <- 2*(p@mllk + np)
    p@BIC <- 2*p@mllk + np*log(n)
    
    if (decode) return(ldhmm.decoding(p,x))
    else return (p)
}
### <---------------------------------------------------------------------->
