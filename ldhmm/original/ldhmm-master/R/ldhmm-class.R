#' The ldhmm class
#' 
#' This S4 class is the major object class for ldhmm package
#'
#' @name ldhmm-class
#'
#' @slot call The match.call slot
#' @slot m numeric, length 1, number of states
#' @slot param.nbr numeric, number of parameters (2 or 3) for each ecld object
#' @slot param matrix, natural parameters for ecld objects, size of states times param.nbr.
#'             Each row can be 2-parameter sequences, or 3-parameter sequences.
#'             Three-parameter unit \code{(mu, sigma, lambda)} forms an ecld object
#'             representing a leptokurtic symmetric lambda distribution.
#'             On the other hand, to provide compatibility to a normal distribution HMM,
#'             two-parameter unit \code{(mu, sigma)} forms an ecld object with lambda=1.
#' @slot gamma matrix, the transition probability matrix, must be m by m.
#' @slot delta numeric, the initial distribution for each state, default is \code{NULL}.
#' @slot stationary logical, specify whether the initial distribution is stationary or not,
#'                  default is \code{TRUE}.
#' @slot mle.optimizer character, the MLE optimizer. Currently it is just set to "nlm".
#' @slot return.code numeric, the return code from the MLE optimizer.
#' @slot iterations numeric, number of iterations MLE optimizer takes.
#' @slot mllk numeric, the final mllk value.
#' @slot AIC numeric, the final AIC.
#' @slot BIC numeric, the final BIC.
#' @slot observations numeric, stores the observations post optimization
#' @slot states.prob matrix, stores the state probabilities post optimization
#' @slot states.local numeric, stores the local decoding states post optimization
#' @slot states.global numeric, stores the global decoding states post optimization (Viterbi)
#' @slot states.local.stats matrix, stores the statistics of local states post optimization
#' @slot states.global.stats matrix, stores the statistics of global states post optimization
#'
#' @keywords class constructor
#'
#' @include ldhmm-package.R
#' @include ldhmm-numericOrNull-class.R
#'
#' @exportClass ldhmm
setClass("ldhmm",
         representation(call = "call",
                        m  = "numeric",
                        param.nbr = "numeric",
                        param = "matrix",
                        gamma = "matrix",
                        delta = "numericOrNull",
                        stationary = "logical",
                        mle.optimizer = "character",
                        return.code = "numericOrNull",
                        iterations = "numericOrNull",
                        mllk = "numericOrNull",
                        AIC = "numericOrNull",
                        BIC = "numericOrNull",
                        observations = "numeric",
                        states.prob = "matrix",
                        states.local = "numeric",
                        states.global = "numeric",
                        states.local.stats = "matrix",
                        states.global.stats = "matrix"),
          prototype(call = call("ldhmm"),
                    m = 2,
                    delta = NULL,
                    stationary = TRUE,
                    mle.optimizer = "nlm",
                    return.code = NULL,
                    iterations = NULL,
                    mllk = NULL,
                    AIC = NULL,
                    BIC = NULL)
)
### <---------------------------------------------------------------------->
