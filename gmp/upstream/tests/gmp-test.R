library(gmp)

##
##' @title Test a unary (if unary=TRUE) or *binary* function
##' @param FUN a function, such as add.bigq() ...
##' @param x a list of "numbers"
##' @param out string determining output class; if "str", use characters, otherwise double
##' @return
##' @author Antoine Lucas (& Martin Maechler)
##' @examples test(as.bigq, 0)
test <- function(FUN, x, xlabs, out = "str", unary = FALSE)
{
  if(missing(xlabs))
    xlabs <- if(is.character(names(x))) names(x) else sapply(x, formatN)
  stopifnot(is.function(FUN), is.list(x),
	    (n <- length(x)) >= 1, length(xlabs) == n)
  if(out == "str") {
      sortie <- as.character
      res <- ""
      error <- "error"
  } else {
      sortie <- as.double
      res <- 0
      error <- NA
  }
  nr <- if(unary) 1 else n
  xlabs <- gsub(" ", "", xlabs)
  res <- matrix(res, nr, n,
                dimnames = list(if(!unary) abbreviate(xlabs, 11, named=FALSE), xlabs))
  for(i in 1:nr){
    classNameI = class(x[[i]])
    for(j in 1:n) {
      classNameJ = class(x[[j]])

      e <- if(unary) tryCatch(FUN(x[[j]]),        error=identity) else
                     tryCatch(FUN(x[[i]],x[[j]]), error=identity)
      if(inherits(e, "error"))
        e <- error
      else if(length(e) == 0)
        e <- numeric()
      ## we don't test standard R floating operations.
      if( (classNameI[1] == "numeric" || classNameI[1] == "integer") && ( classNameJ[1] == "numeric" || classNameJ[1] == "integer") && class(e)[1] == "numeric") e <- "-"

      ## ## now, for some functions also compute the corresponding numeric values
      if(length(e) > 0 && is.double(e[1]) && is.finite(e[1]))
        e <- format(signif(e[1], digits=14), digits=7) # signif(), not round()

      res[i,j] <- sortie(e)[1]
    }
  }
  res ## for printing, the user may prefer as.data.frame(.)
}## end{test}


allfunctionid <- c("as.bigz","+","-","*",
		   "divq.bigz","/","%%","^",
		   "inv.bigz", "gcd.bigz", "gcdex", "lcm.bigz",
		   "as.bigq",
		   "chooseZ",
		   "max","min","|","&","xor","c","cbind","rbind")
unaryfunctionid <- c("log","log2","log10","c",
		     "isprime","nextprime", "factorialZ",
		     "sizeinbase","fibnum","fibnum2","lucnum","lucnum2",
		     "factorize","abs","!")
numericFunName <- function(gmpName) {
  if(gmpName != (r <- sub("[ZQ]$","", gmpName)) &&
     r!="as" && existsFunction(r)) # e.g. chooseZ
    return(r)
  if(gmpName != (r <- sub("\\.big[zq]$","", gmpName)) &&
     r!="as" && r!="sub" && existsFunction(r))
    return(r)
  ttt <- c("add" = "+",
           "sub" = "-",
           "mul" = "*",
           "pow" = "^",
           "div" = "/",
           "divq" = "%/%",
           "mod" = "%%")
  if(!is.na(t.r <- ttt[r]))
    t.r[[1L]]
  else ## return argument
    gmpName
}


options(width = 140, nwarnings = 10000)

sapply(allfunctionid,   numericFunName)
sapply(unaryfunctionid, numericFunName)


ex <- expression(23,as.bigz(23),as.bigq(23),c(3,23),as.bigz(c(3,23)),as.bigq(c(3,23)), "25", 2.3, -4, 4L, 0, as.bigz(34),
                 as.bigq(32,7), as.bigz(31,45), NULL,NA, -3L)## TODO:  as.bigz(3)^700
x <- lapply(ex, eval)

## Those "numbers" in x for which arithmetic should also work in double precision:
## not modulo-arithmetic, not larger than double.prec
useN <- sapply(x, function(u) is.null(u[1]) || is.na(u[1]) ||
               (is.finite(as.numeric(u[1])) && (!inherits(u[1], "bigz") || is.null(modulus(u[1])))))
names(x) <- vapply(ex, format, "")
if(FALSE)## shorter & easier {but *not* the original calls from 'ex'}
    names(x) <- sapply(x, formatN)
str(x)
x. <- x[useN]
nx <- lapply(x., as.numeric)
gmp.NS <- asNamespace("gmp")# also get namespace *hidden* functions, i.e. methods:
for(fid in allfunctionid)
  {
    cat ("------------------------------------------\n", fid," ", sep="")
    FUN <- get(fid, envir = gmp.NS, mode="function")
    rc   <- test(FUN, x )
    res  <- test(FUN, x. , out = "numeric")
    if((nfid <- numericFunName(fid)) != fid || existsFunction(nfid, where=baseenv())) {
      FUN <- get(nfid, envir = gmp.NS, mode="function")
      if(nfid != fid) cat("-> num.fn.:", nfid)
      nres <- test(FUN, nx, out = "numeric")
      cat("\n-> all.equal(target = res, current = F(<numeric x>)): ",
          all.equal(res, nres), "\n")
    } else cat("\n\n")
    print(as.data.frame(rc)); cat("\n")
    ##    ^^^^^^^^^^^^^ (for now, to diminuish difference to last version )
  }

summary(warnings()) # ideally *not* platform dependent

##==============================================================================

for(fid in unaryfunctionid)
  {
    cat ("------------------------------------------\n", fid, "\n\n", sep="")
    FUN <- get(fid, envir = gmp.NS, mode="function")
    print(as.data.frame(test(FUN, x, unary=TRUE)))
  }

##==============================================================================

###----------- matrix -----------------------------
x  <- matrix(1:6,3)
stopifnot(identical(as.bigz(x), matrix(as.bigz(as.vector(x)), 3)),
          dim(x) == 3:2,
          dim(x) == dim(ym   <- as.bigz(x, 6:1)),
          dim(x) == dim(ymr  <- as.bigz(x, 4:6)),
          dim(x) == dim(ymc  <- as.bigz(x, 4)),
          dim(x) == dim(ymq  <- as.bigq(x)),
          dim(x) == dim(y    <- as.bigq(x, 6:1))
          ,
          apply(ym,1,max) == 1:3,
          apply(ym,2,min) == c(1,0))

x %*% t(x)

ym %*% t(ym)
ym %*% t(ymr)
ymc %*% t(ymc)
ymq %*% t(ymq)
y %*% t(y)

dd <- dim(D  <- diag(1:4))
stopifnot(dd == dim(Dmq  <- as.bigq(D)),
          dd == dim(Dz <- as.bigz(D)),
          dd == dim(Dm   <- as.bigz(D,6:1)),
          dd == dim(Dmr  <- as.bigz(D,7)),
          dd == dim(Dmc  <- as.bigz(D,4)),
          TRUE)
solve(D)
solve(Dmq)
solve(Dmr)
tools::assertError(solve(Dmc))# Error: argument has no inverse
tools::assertError(solve(Dm)) # Error: System is singular

(D.D <- D %*% t(Dm))# now [>= Jan.2012] works too
vq <- as.bigq(1:4, 4)
r41 <- cbind(as.bigq((1:4)^2, 4))
stopifnot(identical(D.D, tcrossprod(D,Dm)),
	  dim(r41) == c(4,1),
	  identical(r41, Dz %*% vq), ## bigz %*% bigq - used to fail
	  identical(r41, crossprod(Dz, vq))## ditto
	  )

##
## some specific tests

factorize("33162879029270137")

factorize(15959989)

## assignation
x = as.bigz(1:8)
x[3:2] = 9:10
x

x = as.bigz(matrix(1:12,3))
x[3:2,] = 1:8
x
x[,2] = 0
x

tools::assertError(x[,5])

