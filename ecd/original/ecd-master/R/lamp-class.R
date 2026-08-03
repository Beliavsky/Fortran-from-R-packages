#' An S4 class to represent the lambda process
#' 
#' The \code{lamp} class serves as an object-oriented interface for the lambda process.
#' The main purpose of the class is to store all the parameters required for simulation.
#'
#' @slot call   the match.call slot.
#' @slot lambda numeric, lambda index of lambda process, which is 2/alpha.
#' @slot alpha  numeric, stable alpha. This is derived from lambda for convenience reason.
#' @slot beta   numeric, stable beta.
#' @slot pm     numeric, parameterization, default to 1.
#' @slot rnd.walk numeric, Random walk method. Default is 1.
#' @slot sd     numeric, standard deviation adjustment. No adjustment if NaN.
#' @slot sd.method numeric, methodology of sd adjustment. 0 means in scale parameter, 1 means in Levy sums.
#' @slot T.inf  numeric, the infinite bound to cut off the Levy sums.
#' @slot rnd.n  numeric, the length of one rnd call.
#' @slot N.lower  numeric, the lower bound of N to truncate the boundary effect. Default is 0.
#' @slot N.upper  numeric, the upper bound of N to limit the outliers. Default is 1000.
#' @slot use.mpfr logical, use Mpfr for high precision sums.
#' @slot file  character, file path to save the object and simulation result.
#' @slot tau   numeric, storage for the stable random variables.
#' @slot tau_i numeric, for internal use, length or index of tau.
#' @slot Z_i  numeric, length of Z.
#' @slot Z    numeric, simulation result of the lambda process, Z.
#' @slot B    numeric, simulation result of the binomial process, B.
#' @slot N    numeric, simulation result of the count process, N.
#' @slot tm   POSIXct, timestamp of simulation.
#' @include ecd-class.R
#' @keywords lamp
#' 
#' @author Stephen H. Lihn
#'
#' @exportClass lamp
setClass("lamp",
         representation(call = "call",
                        lambda = "numeric",
                        alpha = "numeric",
                        beta = "numeric",
                        pm = "numeric",
                        rnd.walk = "numeric",
                        sd = "numeric",
                        sd.method = "numeric",
                        T.inf = "numericMpfr",
                        rnd.n = "numeric",
                        N.lower = "numeric",
                        N.upper = "numeric",
                        use.mpfr = "logical",
                        file = "character",
                        tau = "numeric",
                        tau_i = "numeric",
                        Z_i = "numeric",
                        Z = "numeric",
                        B = "numeric",
                        N = "numeric",
                        tm = "POSIXct"),
          prototype(call = call("lamp"),
                    lambda = 4,
                    alpha = 1/2,
                    beta = 1,
                    pm = 1,
                    rnd.walk = NaN,
                    sd = NaN,
                    sd.method = 0,
                    T.inf = 86400*1000,
                    rnd.n = 1000000,
                    N.lower = 0,
                    N.upper = 1000,
                    use.mpfr = FALSE,
                    file = character(0),
                    tau = numeric(0),
                    tau_i = 0,
                    Z_i = 0,
                    Z = numeric(0),
                    B = numeric(0),
                    N = numeric(0),
                    tm = Sys.time())
)
### <---------------------------------------------------------------------->
