#' Creates Stochastic mcnode Objects
#'
#' Creates a \code{\link{mcnode}} object using a random generating function.
#'
#' Any function that accepts vectors/matrices as arguments may be used (notably all current
#' random generators of the \pkg{stats} package). Arguments may be sent classically, but it
#' is strongly recommended to use consistent \samp{mcnode}s if arguments should be recycled,
#' since complex recycling is handled for \samp{mcnode} but not for vectors.
#'
#' Compatibility rules for \samp{mcnode} arguments:
#' \describe{
#'   \item{\samp{type="V"}}{accepts \samp{"0" mcnode} of dim \samp{(1 x 1 x nvariates)} or
#'     \samp{(1 x 1 x 1)} (recycled) and \samp{"V" mcnode} of dim \samp{(nsv x 1 x nvariates)}
#'     or \samp{(nsv x 1 x 1)} (recycled).}
#'   \item{\samp{type="U"}}{accepts \samp{"0" mcnode} of dim \samp{(1 x 1 x nvariates)} or
#'     \samp{(1 x 1 x 1)} (recycled) and \samp{"U" mcnode} of dim \samp{(1 x nsu x nvariates)}
#'     or \samp{(1 x nsu x 1)} (recycled).}
#'   \item{\samp{type="VU"}}{accepts \samp{"0"}, \samp{"V"}, \samp{"U"} and \samp{"VU"} nodes
#'     with appropriate recycling rules.}
#'   \item{\samp{type="0"}}{accepts \samp{"0" mcnode} of dim \samp{(1 x 1 x nvariates)} or
#'     \samp{(1 x 1 x 1)} (recycled).}
#' }
#'
#' If \samp{rtrunc=TRUE}, the distribution is truncated on \samp{(linf, lsup]}. The function
#' \samp{func} should have a \samp{q} form and a \samp{p} form.
#'
#' If \samp{lhs=TRUE}, a Latin Hypercube Sampling is used on \samp{nsv} and \samp{nsu}.
#' The function \samp{func} should have a \samp{q} form. Not allowed with multivariate
#' distributions.
#'
#' @param func a function providing random data or its name as character.
#' @param type the type of \samp{mcnode} to be built. By default, a \samp{"V"} node.
#'   See \code{\link{mcnode}} for details.
#' @param ... all other arguments but the size of the sample to be passed to \samp{func}.
#'   These arguments should be vectors or \samp{mcnode}s (arrays prohibited).
#' @param nsv the number of simulations in the variability dimension.
#' @param nsu the number of simulations in the uncertainty dimension.
#' @param nvariates the number of variates of the output.
#' @param outm the output of the \samp{mcnode} for multivariate nodes. May be \samp{"each"}
#'   (default), \samp{"none"}, or a vector of function names (as character strings) applied
#'   on the variates dimension before output (e.g. \samp{"mean"}, \samp{c("min","max")}).
#' @param nsample the name of the parameter of the function giving the size of the vector.
#'   By default \samp{"n"}, as in most distributions of the \pkg{stats} library.
#' @param seed the random seed used for the evaluation. If \samp{NULL} the seed is unchanged.
#' @param rtrunc should the distribution be truncated? See \code{\link{rtrunc}}.
#' @param linf if truncated: lower limit. May be a scalar, an array or a \samp{mcnode}.
#' @param lsup if truncated: upper limit. May be a scalar, an array or a \samp{mcnode}.
#'   \samp{lsup} must be pairwise strictly greater than \samp{linf}.
#' @param lhs should a Latin Hypercube Sampling be used? See \code{\link{lhs}}.
#' @return An \samp{mcnode} object.
#' @seealso \code{\link{mcnode}} for a description of \samp{mcnode} objects.
#'   \code{\link{Ops.mcnode}} for operations on \samp{mcnode} objects.
#'   \code{\link{rtrunc}} for important warnings on the \samp{rtrunc} option.
#' @keywords methods
#' @examples
#' Oldnvar <- ndvar()
#' Oldnunc <- ndunc()
#' ndvar(5)
#' ndunc(4)
#'
#' ## Compatibility with mcdata as arguments
#' x0 <- mcstoc(runif, type = "0")
#' xV <- mcstoc(runif, type = "V")
#' xU <- mcstoc(runif, type = "U")
#' xVU <- mcstoc(runif, type = "VU")
#'
#' ## "V" accepts "0" and "V" mcdata
#' mcstoc(rnorm, type = "V", mean = x0, sd = xV)
#'
#' ## "VU" accepts "0", "U", "V" and "VU" with correct recycling
#' mcstoc(rnorm, type = "VU", mean = xV, sd = xU)
#'
#' ndvar(Oldnvar)
#' ndunc(Oldnunc)
#' @export
mcstoc <- function(func=runif, type=c("V","U","VU","0"), ..., nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each", nsample="n", seed=NULL, rtrunc=FALSE, linf=-Inf, lsup=Inf, lhs=FALSE)
{

    func <- match.fun(func)
    if(!is.null(seed)) set.seed(seed)

    if(!is.character(outm)  ||
       !(all(outm %in% c("none","each"))) && !all(sapply(outm, exists, mode="function")))
      stop("outm should be 'none','each' or a vector of name(s) of valid function(s)")

    type <- match.arg(type)
    argsd <- list(...)
    dimf <- switch(type, "V"=c(nsv,1,nvariates),"U"=c(1,nsu,nvariates),"VU"=c(nsv,nsu,nvariates),"0"=c(1,1,nvariates))
    nsv <- dimf[1]
    nsu <- dimf[2]
    nva <- dimf[3]

    if(rtrunc) argsd <- c(argsd,list(linf=linf),list(lsup=lsup))          # launch linf and lsup in the process

    largsd <- length(argsd)

#### A function to deal mcnodes (including linf and lsup) as arguments

    LAFUNC <- function(argsd,typethismc){
          if(!is.null(typethismc)){                                    #mcnode as arguments

            if(!(type=="VU" || typethismc=="0" || typethismc==type)) stop("Incompatible type of nodes") # incompatible node

            dimm <- dim(argsd)
            if ((typethismc == "V" && dimm[1] != nsv) ||
                (typethismc == "U" && dimm[2] != nsu) ||
                (typethismc == "VU" && (dimm[1] != nsv || dimm[2] != nsu)))
                stop("Nodes of incompatible dimensions")
 # incompatible dimension

            if(maxdim3 > 1){  #at least one multivariate node as parameter, need recycling on the third dimension
              if(typethismc=="U") argsd <- apply(argsd, 3, matrix, nrow=maxdim1, ncol=maxdim2, byrow=TRUE)  # recycling U as matrix (maxdim1*maxdim2) x nvariates
                else {
                  if(maxdim1 ==1 && maxdim2 ==1) argsd <- matrix(argsd, nrow=1) # Very special case to be added
                  else argsd <- apply(argsd, 3, matrix, nrow = maxdim1, ncol = maxdim2)
                  }        # recycling 0, V, VU as matrix (maxdim1*maxdim2) x nvariates
            }

            else { dim(argsd) <- NULL    # as vector
                   if(typethismc == "U" && maxdim1!=1) argsd <- rep(argsd, each = maxdim1)     #recycling U as vector nsv*nsu
                    }
          }
          else if(is.array(argsd)) stop("Array prohibited in mcstoc as parameter. Use an mcnode instead")
      return(unclass(argsd))
    }

####

      typemc <- lapply(argsd, attr, which = "type")
      yamc <-   !is.null(unlist(typemc))                                        # At least one mcnode
      if(yamc){
            # evaluate the minimal common to build the minimal recycling level...
            maxdim1 <- unlist(lapply(argsd, function(x) dim(x)[1]))
            maxdim1 <- ifelse(is.null(maxdim1), 1, max(maxdim1))
            maxdim2 <- unlist(lapply(argsd, function(x) dim(x)[2]))
            maxdim2 <- ifelse(is.null(maxdim2), 1, max(maxdim2))
            maxdim3 <- unlist(lapply(argsd, function(x) dim(x)[3]))
            maxdim3 <- ifelse(is.null(maxdim3), 1, max(maxdim3))
      }

    if(largsd != 0){
      argsd <- mapply(LAFUNC, argsd, typemc, SIMPLIFY=FALSE)
      }
#####################################
# If lhs or rtrunc, redefine the function to draw random variables
#########
  #keep a copy of original function
  funcorigin <- func

  if(lhs || rtrunc){                                             #define good function for the random sampling
    distr <- as.character(match.call()$func)                     #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                              #remove the r
    qfun <- paste("q",distr,sep="")                              #define "qfunc"

  if(rtrunc){
       pfun <- paste("p",distr,sep="")                            #define pfunc

       func <- function(...){
          argsd <- list(...)
          nnfin <- argsd[[nsample]]
          linf <- if(length(argsd$linf) <= nnfin) as.vector(argsd$linf) else rep(argsd$linf, length.out=nnfin)      # solve a problem when linf was multivariate
          lsup <- if(length(argsd$lsup) <= nnfin) as.vector(argsd$lsup) else rep(argsd$lsup, length.out=nnfin)      # solve a problem when linf was multivariate
          lmax <- max(length(linf),length(lsup))
          if(any(rep(linf, length.out = lmax) >= rep(lsup, length.out = lmax))) stop("linf should be < lsup")  #recycle vectors
          argsd$linf <- argsd$lsup <- argsd[[nsample]] <- NULL
          #find the p of the limit
          pinf <- as.vector(do.call(pfun,c(list(q=linf),argsd),quote=TRUE))
          psup <- as.vector(do.call(pfun,c(list(q=lsup),argsd),quote=TRUE))
          #sample uniformely between the limits
          if(!lhs) lesp <- runif(nnfin,min=pinf,max=psup)
          else     lesp <- lhs(distr="runif", nsv=dimf[1], nsu=dimf[2], nvariates=dimf[3], min=pinf, max=psup)
          #get the q
          data <- (do.call(qfun,c(list(p=lesp),argsd)))[1:nnfin]
          data[pinf==0 & data > lsup] <- NaN          #ex: rtrunc("lnorm",10,linf=-2,lsup=-1)
          data[psup==1 & data < linf] <- NaN          #ex: rtrunc("unif",10,linf=2,lsup=4,max=1)
          data[is.na(linf) | is.na(lsup)] <- NaN      #ex: rtrunc("norm",10,sd=-2)

          #Two tests for extreme situations. None Catch all possibilities. THe error is first to avoid the warning
          if(any(data <= linf | data > lsup, na.rm=TRUE)) stop("Error in rtrunc: some values are not in the expected range (maybe due to rounding errors)")
          if(isTRUE(all.equal(pinf,1)) | isTRUE(all.equal(psup,0)) ) warning("Warning: check the results from rtrunc. It may have reached rounding errors")
          return(data)
       }#end redefinition func
    }
    else func <- function(...) {                      # LHS only
          argsd <- list(...)
          argsd[[nsample]] <- NULL
          lesp <- lhs(distr="runif", nsv=dimf[1], nsu=dimf[2], nvariates=dimf[3], min=0, max=1)
          return(do.call(qfun,c(list(p=lesp),argsd)))}
  }

  # do a try to test the length if nvariates != 1
  if(nvariates != 1){
    if(largsd != 0) argsdtest <- mapply(function(x,typemc){
                                          if(is.null(typemc)) return(unclass(x))
                                          if(is.matrix(x))    return(x[1,,drop=FALSE])     # mc (they have been unclassed)
                                          return(x[1])},
                                        argsd, typemc, SIMPLIFY=FALSE)
      else argsdtest <- vector(mode="list",length=0)

    argsdtest[[nsample]] <- 1
    if(rtrunc) argsdtest$linf <- argsdtest$lsup <- NULL
    dimf <- c(1,1,1)
    data <- do.call(funcorigin,argsdtest,quote=TRUE)
    l <- length(data)
    if(l == nvariates) {
        if(rtrunc | lhs) stop("mcstoc does not handle rtrunc and lhs for multivariate distributions")
          dimf <- c(nsv,nsu,1)}      # If it returns a vector
     else if(l == 1) dimf <- c(nsv,nsu,nvariates)      # if it returns a number
          else stop("the function should return a vector of size 1 or nvariates if",nsample,"=1")
        argsd[[nsample]] <- prod(dimf)
        data <- do.call(func, argsd, quote = TRUE)

       #Post Production, multivariate
       if (yamc){
          if(l==1){ # univariate distribution
                if(maxdim1 == 1 && maxdim2 == 1) data <- aperm(array(data, dim = c(nvariates, nsv, nsu)), c(2, 3, 1))
           else if(maxdim1 == 1 && nsv!=1)       data <- aperm(array(data, dim = c(nsu, nvariates, nsv)), c(3, 1, 2))
           else if(maxdim2 == 1 && nsu!=1)       data <- aperm(array(data, dim = c(nsv, nvariates, nsu)), c(1, 3, 2))
           else data <- array(data, dim = c(nsv, nsu, nvariates))
          }
          else {     # l != 1 : multivariate
                if(maxdim1 == 1 && nsv != 1)     data <- aperm(array(data, dim = c(nsu, nsv, nvariates)), c(2, 1, 3))
           else data <- array(data, dim = c(nsv, nsu, nvariates))
           }
        }
        else data <- array(data, dim = c(nsv, nsu, nvariates))
      } #end multivariates
   else{  # univariate
        argsd[[nsample]] <- prod(dimf)
        data <- do.call(func, argsd, quote = TRUE)

        if (yamc && maxdim1 == 1 && nsv != 1) data <- aperm(array(data, dim = c(nsu, nsv, nvariates)), c(2, 1, 3))
        else data <- array(data, dim = c(nsv, nsu, nvariates))
   }

      class(data) <- "mcnode"
      attr(data,"type") <- type
      attr(data,"outm") <- outm
      return(data)
    }
