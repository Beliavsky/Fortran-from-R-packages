is.bigz <- function(x) is.raw(x) && inherits(x, "bigz")
is.bigq <- function(x) is.raw(x) && inherits(x, "bigq")

setGeneric("asNumeric", useAsDefault = function(x) {
    if(is.numeric(x)) x else if(is.atomic(x)) {
        storage.mode(x) <- "numeric"; x }
    else as(x, "numeric")
})

#----------------------------------------------------------
#
#  Author        : Immanuel Scholz (immanuel.scholz@gmx.de)
#		   Technische Universitaet Dresden
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL (>= 2)
#
#----------------------------------------------------------

add.bigz <- function(e1, e2) {
    if(inherits(e2, "bigq"))
        .Call(bigrational_add, e1, e2)
    else .Call(biginteger_add, e1, e2)
}

sub.bigz <- function(e1, e2=NULL)
{
    if(is.null(e2))
        .Call(biginteger_sub, 0, e1)
    ## else if(inherits(e2, "bigq"))
    ##     .Call(bigrational_sub, e1, e2)
    else
        .Call(biginteger_sub, e1, e2)
}

mul.bigz <- function(e1, e2) {
    if(inherits(e2, "bigq"))
        .Call(bigrational_mul, e1, e2)
    else .Call(biginteger_mul, e1, e2)
}

## divq : integer division
"%/%.bigz" <- divq.bigz <- function(e1, e2) {
   if(inherits(e2, "bigq")) {
       if(!all(is.whole(e2[is.finite(e2)])))
           e2 <- as.bigz(e2)
       else
           stop("In 'n %/% d', d must be integer")
   }
   .Call(biginteger_divq, e1, e2)
}

## div : division of integers -> either rational or (mod) integer division
div.bigz <- function(e1, e2) {
    if(inherits(e2, "bigq"))
        .Call(bigrational_div, e1, e2)
    else .Call(biginteger_div, e1, e2)
}

"%%.bigz" <- mod.bigz <- function(e1, e2) {
   if(inherits(e2, "bigq")) {
       if(!all(is.whole.bigq(e2[is.finite(e2)])))
           e2 <- as.bigz(e2)
       else
           stop("In 'n %% d', d must be integer")
   }
   .Call(biginteger_mod, e1, e2)
}
.mod.bigz <- function(e1, e2) .Call(biginteger_mod, e1, e2)

pow.bigz <- function(e1, e2,...) {
    if(inherits(e2, "bigq"))
        pow.bigq(e1, e2)
    else .Call(biginteger_pow, e1, e2)
}

##' Inverse:  inv(a,b) := (1 / a) (modulo b)
inv.bigz <- function(a,b,...) .Call(biginteger_inv,a,b)

"!.bigz" <- function(a) a == 0


## as.boolean(x): x != 0

"|.bigz" <- function(a,b) {
	 a1 = a != 0
	 b1 = b != 0
	 a1 | b1
}

"&.bigz" <- function(a,b) {
	 a1 = a != 0
	 b1 = b != 0
	 a1 & b1
}

"xor.bigz" <- function(x,y) {
	 a1 = x != 0
	 b1 = y != 0
	 xor(a1 , b1)
}

gcd <- function(a,b)
      UseMethod("gcd")
gcd.default <- function(a,b) as.integer(gcd.bigz(a,b))
gcd.bigz <- function(a,b) .Call(biginteger_gcd,a,b)

## just because lcm() is a trivial function in 'graphics' .. hmm
##lcm <- function(a,b)
##      UseMethod("lcm")

lcm.default <- function(a,b)
  as.integer(lcm.bigz(a,b))
lcm.bigz <- function(a,b) .Call(biginteger_lcm,a,b)

print.bigz <- function(x, quote = FALSE, initLine = is.null(modulus(x)), ...)
{
  if((n <- length(x)) > 0) {
    if(initLine) {
      cat("Big Integer ('bigz') ")
      kind <- if(!is.null(nr <- attr(x, "nrow")))
        sprintf("%d x %d matrix", nr, n/nr)
      else if(n > 1) sprintf("object of length %d", n) else ""
      cat(kind,":\n", sep="")
    }
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigz(0)\n")
  invisible(x)
}


as.bigz <- function(a, mod = NA)
{
  if(isZ <- missing(mod) && inherits(a, "bigz"))
    mod <- modulus(a) # possibly NULL
  if(is.null(mod)) mod <- NA
  if(!isZ && inherits(a, "bigq"))
    as.bigz.bigq(a, mod)
  else
    .Call(biginteger_as, a, mod)
}

## the .as*() functions are exported for Rmpfr
.as.bigz <- function(a, mod = NA) {
  if(inherits(a, "bigq")) as.bigz.bigq(a, mod) else .Call(biginteger_as, a, mod)
}
..as.bigz <- function(a, mod = NA) .Call(biginteger_as, a, mod)

.as.char.bigz <-
as.character.bigz <-
    function(x, b = 10L, ...) .Call(biginteger_as_character, x, b)


##' format() Numbers such as to distinguish  bigz, integer, double, mpfr, etc
formatN <- function(x, ...) UseMethod("formatN")
formatN.integer <- function(x, ...) paste0(as.character(x, ...), "L")
formatN.bigz    <- function(x, ...) {
    r <- as.character(x, ...)
    if(any(noMod <- is.null(modulus(x))))
	r[noMod] <- paste0(r[noMod],"_Z")
    r
}
formatN.double	<- function(x, ...) {
    r <- vapply(x, format, "", ...)
    if(any(intLike <- !grepl("[^-0-9]",r)))
	r[intLike] <- paste0(r[intLike],".")
    r
}
##' Default Method: Use the standard format() --- e.g. for complex
formatN.default <- function(x, ...) format(x, ...)


as.double.bigz  <- function(x,...) .Call(biginteger_as_numeric, x)
as.integer.bigz <- function(x,...) .Call(biginteger_as_integer, x)

.bigz2num <- function(x) {
    r <- .Call(biginteger_as_numeric, x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}
setMethod("asNumeric", "bigz", .bigz2num)


length.bigz <- function(x) .Call(biginteger_length, x)

"length<-.bigz"<- function(x, value) .Call(biginteger_setlength, x, value)

modulus      <- function(a) UseMethod("modulus")
modulus.bigz <- function(a) attr(a, "mod")

`modulus<-`      <- function(a, value) UseMethod("modulus<-")
`modulus<-.bigz` <- function(a, value) as.bigz(a, value)


## inv <- function(a,...) UseMethod("inv")

## pow <- function(a,...) UseMethod("pow")

powm <- function(x,y, n) .Call(biginteger_powm, x,y,n)

## <op>.bigz(): *not* used
lt.bigz  <- function(e1, e2) .Call(biginteger_lt, e1, e2)
gt.bigz  <- function(e1, e2) .Call(biginteger_gt, e1, e2)
lte.bigz <- function(e1, e2) .Call(biginteger_lte, e1, e2)
gte.bigz <- function(e1, e2) .Call(biginteger_gte, e1, e2)
eq.bigz  <- function(e1, e2) .Call(biginteger_eq, e1, e2)
neq.bigz <- function(e1, e2) .Call(biginteger_neq, e1, e2)

lt.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_lt, e1, e2)
    else
        .Call(bigrational_lt, e1, e2)
}
gt.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_gt, e1, e2)
    else
        .Call(bigrational_gt, e1, e2)
}

lte.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_lte, e1, e2)
    else
        .Call(bigrational_lte, e1, e2)
}
gte.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_gte, e1, e2)
    else
        .Call(bigrational_gte, e1, e2)
}
eq.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_eq, e1, e2)
    else
        .Call(bigrational_eq, e1, e2)
}
neq.big <- function(e1, e2) {
    if(!is.bigq(e1) && !is.bigq(e2)) # try integer
        .Call(biginteger_neq, e1, e2)
    else
        .Call(bigrational_neq, e1, e2)
}


is.whole <- function(x) UseMethod("is.whole")
is.whole.default <- function(x) {
    n <- length(x)
    if(is.atomic(x)) {
        if(is.integer(x) || is.logical(x))
            return(rep.int(TRUE, n))
        if(is.numeric(x))
            return(x == floor(x))
        if(is.complex(x))
            return(x == round(x))
    }
    ## else:
    logical(n) ## == rep.int(FALSE, length(x))
}


is.na.bigz <- function(x) .Call(biginteger_is_na, x)
is.finite.bigz <- function(x) !is.na.bigz(x) # otherwise all are finite
is.whole.bigz <- function(x) rep.int(TRUE, length(x))
is.infinite.bigz <- function(x) rep.int(FALSE, length(x))

frexpZ <- function(x) .Call(bigI_frexp, x)

##' @title log2(Inverse of frexpZ(a))
##' @param L list(d = ., exp = .)
##' @return numeric vector
##' @author Martin Maechler
lg2.invFrexp <- function(L) {
    stopifnot(is.list(L), is.numeric(d <- L$d), is.numeric(ex <- L$exp),
	      (length(d)) == length(ex))
    ex + log2(d)
}

###------------------------- 'Math' S3 group ------------------------------

## o 'abs', 'sign', 'sqrt',
##   'floor', 'ceiling', 'trunc',
##   'round', 'signif'
## o 'exp', 'log', 'expm1', 'log1p',
##   'cos', 'sin', 'tan',
##   'acos', 'asin', 'atan'
##   'cosh', 'sinh', 'tanh',
##   'acosh', 'asinh', 'atanh'
## o 'lgamma', 'gamma', 'digamma', 'trigamma'
## o 'cumsum', 'cumprod', 'cummax', 'cummin'

## Most 'Math' group functions should go via CRAN package 'Rmpfr' :
Math.bigz <- function(x, ...) {
    if(requireNamespace("Rmpfr", quietly=TRUE)) {
        NextMethod(Rmpfr::.bigz2mpfr(x)) # FIXME use ..bigz2mpfr (two '.') in future
    }
    else
        stop("Math group method ", dQuote(.Generic),
             "is available via CRAN R package 'Rmpfr'.\n",
             "Install it and try again")

}

abs.bigz <- function(x) .Call(biginteger_abs,x)
sign.bigz <- function(x) .Call(biginteger_sgn,x)

floor.bigz <- ceiling.bigz <- function(x) x
trunc.bigz <- function(x, ...) x
round.bigz <- function(x, digits=0) {
    ## round(x * 10^d) / 10^d
    stopifnot(length(digits) == 1L)
    if(digits >= 0)
        x
    else { # digits < 0
        p10 <- as.bigz(10) ^ -digits # still bigz
        round(x / p10) * p10
    }
}

gamma.bigz <- function(x) factorialZ(x-1)

cumsum.bigz <- function(x) .Call(biginteger_cumsum, x)
## TODO: add cummax(), cummin(), cumprod()


log2.bigz <- function(x) .Call(biginteger_log2, x)
## not exported:
ln.bigz	 <- function(x) .Call(biginteger_log, x)
log.bigz <- function(x, base=exp(1))
{
    if(missing(base))
	ln.bigz(x)
    else
	ln.bigz(x)/log(base)
}

log10.bigz <- function(x) ln.bigz(x) / log(10)

##------------end{'Math'} group -------------------------------------


###------------------------- 'Summary' S3 group ------------------------------
##---- "max"   "min"   "range" "prod"  "sum"   "any"   "all"  -----

max.bigz <- function(...,na.rm=FALSE)
{
    X <- c.bigz(...)
    if(inherits(X, "bigq"))
        .Call(bigrational_max, X, na.rm)
    else .Call(biginteger_max, X, na.rm)

}

min.bigz <- function(...,na.rm=FALSE)
{
    X <- c.bigz(...)
    if(inherits(X, "bigq"))
        .Call(bigrational_min, X, na.rm)
    else .Call(biginteger_min, X, na.rm)
}

## range(): works automatically via  range.default() and the above min(), max()

prod.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
    if(inherits(X, "bigq"))
        .Call(bigrational_prod, if(na.rm) X[!is.na(X)] else X)
    else .Call(biginteger_prod, if(na.rm) X[!is.na(X)] else X)
}

sum.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
    if(inherits(X, "bigq"))
        .Call(bigrational_sum, if(na.rm) X[!is.na(X)] else X)
    else .Call(biginteger_sum, if(na.rm) X[!is.na(X)] else X)
}
##------------end{Summary group}------------------------------------

## FIXME: implement faster in C
setMethod("which.max", "bigz", function(x) which.max(x == max(x)))
setMethod("which.min", "bigz", function(x) which.max(x == min(x)))

##' to be applied e.g. to the result of  lapply(<bigz>, Fn)
c_bigz <- function(L) .Call(biginteger_c, L)

c.bigz <- function(..., recursive = FALSE)
{
    argL <- list(...)
    if(any(vapply(argL, inherits, NA, what="bigq")))
        c_bigq(argL)
    else
        c_bigz(argL)
}

## This is practically identical to  grid :: rep.unit :
rep.bigz <- function(x, times=1, length.out=NA, each=1, ...) {
    ## if (length(x) == 0)
    ##	   stop("invalid 'unit' object")
    if(!missing(times) && missing(length.out) && missing(each))
        .Call(biginteger_rep, x, times)
    else {
	## Determine an appropriate index, then call subsetting code
	x[ rep(seq_along(x), times=times, length.out=length.out, each=each) ]
    }
}

duplicated.bigz <- function(x, incomparables = FALSE, ...) {
    x <- as.character(x) # lazy and inefficient --> TODO in C++
    NextMethod("duplicated", x)
}

unique.bigz <- function(x, incomparables = FALSE, ...)
    x[!duplicated(x, incomparables, ...)]


all.equal.bigz <- function(target, current, ...) {
    if(is.bigq(target))
        return(all.equal.bigq(target, current, ...))
    ## Using  tolerance = 0  --- implicitly below
    if(length(target) != length(current))
        return("lengths differ")
    if(!identical(dim(target), dim(current)))
        return("dimensions differ")
    target <- as.vector(target)
    current <- as.vector(current)
    ina <- is.na(target)
    if(any(ina != is.na(current)))
	paste("'is.NA' value mismatch:", sum(is.na(current)),
                     "in current", sum(ina), "in target")
    else if(all(ina | target == current)) # equal NAs _or_ numbers
        TRUE
    else
        "'target'(bigz) and 'current' differ"
}


# Isprime, return:
#   0 if not prime
#   1 if probably prime
#   2 if prime
isprime <- function(n,reps=40)
  {
    .Call(biginteger_is_prime, n, as.integer(reps))
  }

nextprime <- function(n) .Call(biginteger_nextprime, n)

gcdex <- function(a, b) .Call(biginteger_gcdex, a, b)

urand.bigz <- function(nb=1, size=200, seed=0)
  {
    ok <- (seed != 0)
    .Call(biginteger_rand_u,
          as.integer(nb),
          as.integer(size),
          seed,
          as.integer(ok))
  }

sizeinbase <- function(a, b=10)
{
    if(as.integer(b) < 2)
        stop("base must be >= 2")
    .Call(biginteger_sizeinbase, a, as.integer(b))
}

factorialZ  <- function(n) .Call(bigI_factorial, as.integer(n))
chooseZ  <- function(n, k) .Call(bigI_choose,  n, as.integer(k))

fibnum  <- function(n) .Call(bigI_fibnum,  as.integer(n))
fibnum2 <- function(n) .Call(bigI_fibnum2, as.integer(n))

lucnum  <- function(n) .Call(bigI_lucnum,  as.integer(n))
lucnum2 <- function(n) .Call(bigI_lucnum2, as.integer(n))

factorize <- function(n) .Call(factorR, as.bigz(n))

solve.bigz <- function(a, b,...)
  {
    if(missing(b))
      .Call(inverse_z, a)
    else
      .Call(solve_z, a, b)
  }


`[[.bigz` <- function(x, i=NA) .Call(biginteger_get_at, x, i)

`[[<-.bigz` <- function(x, i=NA, value)
    .Call(biginteger_set_at, x, i, value)

`[.bigz` <- function(x, i=NULL, j=NULL, drop=TRUE)
{

  mdrop <- missing(drop)
  Narg <- nargs() - (!mdrop)
  # matrix access [i,j] [,j] [i,]
  # vector access [i]
  matrixAccess = Narg > 2
  has.j <- !missing(j)
  if(!is.null(attr(x, "nrow")) & matrixAccess) { ## matrix
    .Call(matrix_get_at_z, x, i,j)
  } else { ## non-matrix
    if(has.j) stop("invalid vector subsetting")

    r <- .Call(biginteger_get_at, x, i)
    attr(r,"nrow") <- NULL
    r
  }
}

`[<-.bigz` <- function(x, i=NULL, j=NULL, value)
{

  # matrix access [i,j] [,j] [i,]
  # vector access [i]
  matrixAccess = nargs() > 3

  has.j <- !missing(j)
  if(!is.null(attr(x, "nrow")) & matrixAccess) { ## matrix
      .Call(matrix_set_at_z, x, value, i,j)
  } else { ## non-matrix 
    if(has.j) stop("invalid vector subsetting")
    r <- .Call(biginteger_set_at, x, i, value) 
    attr(r,"nrow") <- attr(x, "nrow")
    r
  }
}


