##' Least Squares and Quadratic Programming under Nonnegativity
##' Constraints
##'
##'
##' These functions are particularly useful for solving least squares
##' or quadratic programming problems when some or all of the solution
##' values are subject to nonnegativity constraint. One may further
##' restrict the NN-restricted coefficients to have a fixed positive
##' sum.
##'
##' Function \code{nnls} solves the least squares problem under
##' nonnegativity (NN) constraints. It is an R interface to the NNLS
##' function that is described in Lawson and Hanson (1974, 1995). Its
##' Fortran implementation is public domain and available at
##' \url{https://www.netlib.org/lawson-hanson/} (with slight
##' modifications by Yong Wang for compatibility with the lastest
##' Fortran compiler.)
##'
##' Given matrix \code{a} and vector \code{b}, \code{nnls} solves the
##' nonnegativity least squares problem:
##'
##' \deqn{\mathrm{minimize \ \ } || a x - b ||,}{minimize || a x - b ||,}
##' \deqn{\mathrm{\ \ \ subject\ to\ \ } x \ge 0.}{ subject to x >= 0.}
##'
##' Function \code{pnnls} also solves the above nonnegativity least
##' squares problem when \code{k=0}, but it may also leave the first
##' \code{k} coefficients unrestricted. The output value of \code{k}
##' can be smaller than the input one, if \code{a} has linearly
##' dependent columns. If \code{sum} is a positive value, \code{pnnls}
##' solves the problem by further restricting that the NN-restricted
##' coefficients must sum to the given value.
##'
##' Function \code{pnnqp} solves the quadratic programming problem
##'
##' \deqn{\mathrm{minimize\ \ } \frac12 x^T q x + p^T x,}{minimize 0.5 x^T q x +
##' p^T x,}
##'
##' when only some or all coefficients are restricted by
##' nonnegativity. The quadratic programming problem is solved by
##' transforming the problem into a least squares one under the same
##' constraints, which is then solved by function
##' \code{pnnls}. Arguments \code{k} and \code{sum} have the same
##' meanings as for \code{pnnls}.
##'
##' Functions \code{nnls}, \code{pnnls} and \code{pnnqp} are able to
##' return any zero-valued solution as 0 exactly. This differs from
##' functions \code{lsei} and \code{qp}, which may produce very small
##' values for exactly 0s, thanks to numerical errors.
##'
##'@aliases nnls pnnls pnnqp
##'@param a Design matrix.
##'@param b Response vector.
##'@param k Integer, meaning that the first \code{k} coefficients are not
##'NN-restricted.
##'@param sum = NULL, if NN-restricted coefficients are not further restricted
##'to have a fixed sum;
##'
##'= a positive value, if NN-restricted coefficients are further restricted to
##'have a fixed positive sum.
##'@param q Positive semidefinite matrix of numeric values for the quadratic
##'term of a quadratic programming problem.
##'@param p Vector of numeric values for the linear term of a quadratic
##'programming problem.
##'@param tol Tolerance used for calculating pseudo-rank of \code{q}.
##'@return
##'
##'\item{x}{Solution}
##'
##'\item{r}{The upper-triangular matrix \code{Q*a}, pivoted by variables in the
##'order of \code{index}, when \code{sum=NULL}. If \code{sum > 0}, \code{r} is
##'for the transformed \code{a}.}
##'
##'\item{b}{The vector \code{Q*b}, pivoted by variables in the order of
##'\code{index}, when \code{sum=NULL}. If \code{sum > 0}, \code{b} is for the
##'transformed \code{b}.}
##'
##'\item{index}{Indices of the columns of \code{r}; those unrestricted and in
##'the positive set are first given, and then those in the zero set.}
##'
##'\item{rnorm}{Euclidean norm of the residual vector.}
##'
##'\item{mode}{= 1, successful computation;
##'
##'= 2, bad dimensions of the problem;
##'
##'= 3, iteration count exceeded (more than 3 times the number of variables
##'iterations).}
##'
##'\item{k}{Number of the first few coefficients that are truly not
##'NN-restricted.}
##' 
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{lsei}}, \code{\link{hfti}}.
##'@references
##'
##'Lawson and Hanson (1974, 1995). Solving Least Squares Problems. Englewood
##'Cliffs, N.J., Prentice-Hall.
##'
##'Dax (1990). The smallest point of a polytope. Journal of Optimization Theory
##'and Applications, 64, pp. 429-432.
##'
##'Wang (2010). Fisher scoring: An interpolation family and its Monte Carlo
##'implementations. Computational Statistics and Data Analysis, 54, pp.
##'1744-1755.
##' 
##'@keywords array algebra
##'@examples
##'
##'a = matrix(rnorm(40), nrow=10)
##'b = drop(a %*% c(0,1,-1,1)) + rnorm(10)
##'nnls(a, b)$x                     # constraint x >= 0
##'pnnls(a, b, k=0)$x               # same as nnls(a, b)
##'pnnls(a, b, k=2)$x               # first two coeffs are not NN-constrained
##'pnnls(a, b, k=2, sum=1)$x        # NN-constrained coeffs must sum to 1
##'pnnls(a, b, k=2, sum=2)$x        # NN-constrained coeffs must sum to 2
##'q = crossprod(a)
##'p = -drop(crossprod(b, a))
##'pnnqp(q, p, k=2, sum=2)$x        # same solution
##'
##'pnnls(a, b, sum=1)$x             # zeros found exactly
##'pnnqp(q, p, sum=1)$x             # zeros found exactly
##'lsei(a, b, rep(1,4), 1, lower=0) # zeros not so exact
##'
##'@usage 
##'nnls(a, b) 
##'pnnls(a, b, k=0, sum=NULL) 
##'pnnqp(q, p, k=0, sum=NULL, tol=1e-20)
##'
##'@export nnls
##'@export pnnls
##'@export pnnqp
##'
##'@useDynLib lsei, .registration = TRUE

nnls = function(a, b) {
  if(!is.vector(b)) b = drop(b)
  if(!is.matrix(a)) stop("a not matrix")
  m = nrow(a)
  n = ncol(a)
  if(length(b) != m) stop("length(b) != ncol(a)")
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  .Fortran("nnls",r=a,m,m,n,b=b,x=x,rnorm=rnorm,w,zz,index=index,
           mode=mode,PACKAGE="lsei")[c("x","r","b","index","rnorm","mode")]
}

pnnls = function(a, b, k=0, sum=NULL) {
  if(!is.vector(b)) b = drop(b)
  if(!is.matrix(a)) stop("a not matrix")
  m = nrow(a)
  n = ncol(a)
  if(!is.null(sum)) {
    if(sum <= 0) stop("Argument 'sum' must be positive or NULL")
    if(k<n) a[,(k+1):n] = a[,(k+1):n] * sum - b
    else stop("k >= ncol(a) (null simplex)")
    a = rbind(a, c(double(k), rep(1, n-k)))
    b = c(double(m), 1)
    m = as.integer(m+1)
  }
  if(length(b) != m) stop("length(b) != ncol(a)")
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  k = as.integer(k)
  r = .Fortran("pnnls",r=a,m,m,n,b=b,x=x,rnorm=rnorm,w,zz,index=index,
      mode=mode,k=k,PACKAGE="lsei")
  r$r = r$r[1:min(m,n),]
  if(!is.null(sum)) {
    t = sum(r$x[(r$k+1):n])
    r$x = r$x / t
    r$x[(r$k+1):n] = r$x[(r$k+1):n] * sum
    r$rnorm = sqrt( pmax((r$rnorm/t)^2 - (1 - 1/t)^2, 0) )
  }
  r[c("x","r","b","index","rnorm","mode","k")]
}

# --------------------------------- #
# Least distance programming (LDP): #
#                                   #
#        Minimize   ||x||           #
#        Suject to  e x >= f        #
# --------------------------------- #

## Special treatment: Relax boundaries slightly to ignore negligible
## incompatibility that may be produced numerically.

## Remove constraints with effectively zero coefficients???

ldp = function(e, f) {
  if(!is.vector(f)) f = drop(f)
  if(is.vector(e)) dim(e) = c(1, length(e))

  f = f - max(abs(f)) * 5e-15   # relax boundaries slightly
  m = nrow(e)
  n = ncol(e)                         # number of variables
  storage.mode(e) = "double"
  storage.mode(f) = "double"
  x = double(n)                       # only for output
  xnorm = double(1)                   # only for output
  w = double((n+1)*(m+2) + 2*m)       # working space
  index = integer(m)                  # only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  r = .Fortran("ldp",e,m,m,n,f,x=x,xnorm,w,index,mode=mode,
               PACKAGE="lsei")[c("x","mode")]
  if(r$mode != 1) stop("Incompatible constraints in ldp()")
  r$x
}

ldp2 = function(e, f, tol=1e-15) {
  if(!is.vector(f)) f = drop(f)
  if(is.vector(e)) dim(e) = c(1, length(e))
  f = f - max(abs(f)) * 5e-15   # relax boundaries slightly
  k = ncol(e)                            # number of variables
  E = rbind(t(e), t(f))
  h = c(double(k), 1)
  r = E %*% nnls(E, h)$x - h             # residuals
  if(sqrt(sum(r^2)) <= tol) stop("Incompatible inequalities in ldp()")
  as.vector(-r[1:k] / r[k+1])
}

# # Example from Lawson and Hanson. (1974), p.171:
# e = cbind(c(-.207,-.392,.599), c(2.558, -1.351, -1.206))
# f = c(-1.3,-.084,.384)
# ldp(e, f)           # Solution: 0.1268538 -0.2554018

# G = matrix(rnorm(12), nrow=4); h = rnorm(4); x = ldp(G, h); print(x); G %*% x - h

# --------------------------------------------------------------- #
# Least squares problem with linear inequality constraints (LSI): #
#                                                                 #
#         Minimize       || a x - b ||                            #
#         Subject to     e x >= f                                 #
# --------------------------------------------------------------- #

lsi = function(a, b, e=NULL, f=NULL, lower=-Inf, upper=Inf) {
  if(is.vector(e)) dim(e) = c(1, length(e))
  if(any(lower != -Inf, upper != Inf)) {
    k0 = ncol(a)
    lower = rep(lower, len=k0)
    upper = rep(upper, len=k0)
    jl = lower != -Inf
    ju = upper != Inf
    e = rbind(e, diag(1, k0)[jl,], diag(-1, k0)[ju,])
    f = c(f, lower[jl], -upper[ju])
  }
  if(is.null(e)) return(pnnls(a, b, k=ncol(a))$x)
  a.svd = svdrs(a, b)
  k = sum(a.svd$d > max(a.svd$d) * 1e-14)     # pseudo-rank
  P1b = a.svd$uTb[1:k,,drop=FALSE]
  Q1 = a.svd$v[,1:k,drop=FALSE]
  et = e %*% sweep(Q1, 2, a.svd$d[1:k], "/")
  ft = f - et %*% P1b
  z = ldp2(et, ft)
  drop(Q1 %*% ((z + P1b) / a.svd$d[1:k]))
}

# # Example from Lawson and Hanson. (1974), p.170:
# a = cbind(c(.25,.5,.5,.8),rep(1,4)); b = c(.5,.6,.7,1.2); e = cbind(c(1,0,-1),c(0,1,-1)); f = c(0,0,-1); lsi(a, b, e, f)    # Solution: 0.6213152 0.3786848

# E = matrix(rnorm(24), nrow=8); f = rnorm(8); G = matrix(rnorm(12), nrow=4); h = rnorm(4); x = lsi(E, f, G, h); G %*% x - h

# -------------------------------------------------------------------------- #
# Least squares problem with both linear equalities and inequalities (LSEI): #
#                                                                            #
#          Minimize      || a x - b ||                                       #
#          Subject to    c x = d,  e x >= f                                  #
# -------------------------------------------------------------------------- #

##'Least Squares and Quadratic Programming under Equality and Inequality Constraints
##'
##' These functions can be used for solving least squares or quadratic
##' programming problems under general equality and/or inequality
##' constraints.
##' 
##'The \code{lsei} function solves a least squares problem under both equality
##'and inequality constraints. It is an implementation of the LSEI algorithm
##'described in Lawson and Hanson (1974, 1995).
##'
##'The \code{lsi} function solves a least squares problem under inequality
##'constraints. It is an implementation of the LSI algorithm described in
##'Lawson and Hanson (1974, 1995).
##'
##'The \code{ldp} function solves a least distance programming problem under
##'inequality constraints. It is an R wrapper of the LDP function which is in
##'Fortran, as described in Lawson and Hanson (1974, 1995).
##'
##'The \code{qp} function solves a quadratic programming problem, by
##'transforming the problem into a least squares one under the same equality
##'and inequality constraints, which is then solved by function \code{lsei}.
##'
##'The NNLS and LDP Fortran implementations used internally is downloaded from
##'\url{https://www.netlib.org/lawson-hanson/}.
##'
##'
##'Given matrices \code{a}, \code{c} and \code{e}, and vectors \code{b},
##'\code{d} and \code{f}, function \code{lsei} solves the least squares problem
##'under both equality and inequality constraints:
##'
##'\deqn{\mathrm{minimize\ \ } || a x - b ||,}{minimize || a x - b ||,}
##'\deqn{\mathrm{subject\ to\ \ } c x = d, e x \ge f.}{subject to c x = d, e x
##'>= f.}
##'
##'Function \code{lsi} solves the least squares problem under inequality
##'constraints:
##'
##'\deqn{\mathrm{minimize\ \ } || a x - b ||,}{minimize || a x - b ||,}
##'\deqn{\mathrm{\ \ \ subject\ to\ \ } e x \ge f.}{subject to e x >= f.}
##'
##'Function \code{ldp} solves the least distance programming problem under
##'inequality constraints:
##'
##'\deqn{\mathrm{minimize\ \ } || x ||,}{minimize || x ||,} \deqn{\mathrm{\ \ \
##'subject\ to\ \ } e x \ge f.}{subject to e x >= f.}
##'
##'Function \code{qp} solves the quadratic programming problem:
##'
##'\deqn{\mathrm{minimize\ \ } \frac12 x^T q x + p^T x,}{minimize 0.5 x^T q x +
##'p^T x,} \deqn{\mathrm{subject\ to\ \ } c x = d, e x \ge f.}{subject to c x =
##'d, e x >= f.}
##'
##'@aliases lsei lsi ldp qp
##'@param a Design matrix.
##'@param b Response vector.
##'@param c Matrix of numeric coefficients on the left-hand sides of equality
##'constraints. If it is NULL, \code{c} and \code{d} are ignored.
##'@param d Vector of numeric values on the right-hand sides of equality
##'constraints.
##'@param e Matrix of numeric coefficients on the left-hand sides of inequality
##'constraints. If it is NULL, \code{e} and \code{f} are ignored.
##'@param f Vector of numeric values on the right-hand sides of inequality
##'constraints.
##'@param q Matrix of numeric values for the quadratic term of a quadratic
##'programming problem.
##'@param p Vector of numeric values for the linear term of a quadratic
##'programming problem.
##'@param lower,upper Bounds on the solutions, as a way to specify such simple
##'inequality constraints.
##'@param tol Tolerance, for calculating pseudo-rank in \code{qp}.
##'@return A vector of the solution values
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{nnls}},\code{\link{hfti}}.
##'@references Lawson and Hanson (1974, 1995). Solving least squares problems.
##'Englewood Cliffs, N.J., Prentice-Hall.
##'@keywords array algebra
##'@examples
##'
##'beta = c(rnorm(2), 1)
##'beta[beta<0] = 0
##'beta = beta / sum(beta)
##'a = matrix(rnorm(18), ncol=3)
##'b = a %*% beta + rnorm(3,sd=.1)
##'c = t(rep(1, 3))
##'d = 1
##'e = diag(1,3)
##'f = rep(0,3)
##'lsei(a, b)                        # under no constraint
##'lsei(a, b, c, d)                  # under eq. constraints
##'lsei(a, b, e=e, f=f)              # under ineq. constraints
##'lsei(a, b, c, d, e, f)            # under eq. and ineq. constraints
##'lsei(a, b, rep(1,3), 1, lower=0)  # same solution
##'q = crossprod(a)
##'p = -drop(crossprod(b, a))
##'qp(q, p, rep(1,3), 1, lower=0)    # same solution
##'
##'## Example from Lawson and Hanson (1974), p.140
##'a = cbind(c(.4302,.6246), c(.3516,.3384))
##'b = c(.6593, .9666)
##'c = c(.4087, .1593)
##'d = .1376
##'lsei(a, b, c, d)   # Solution: -1.177499  3.884770
##'
##'## Example from Lawson and Hanson (1974), p.170
##'a = cbind(c(.25,.5,.5,.8),rep(1,4))
##'b = c(.5,.6,.7,1.2)
##'e = cbind(c(1,0,-1),c(0,1,-1))
##'f = c(0,0,-1)
##'lsi(a, b, e, f)      # Solution: 0.6213152 0.3786848
##'
##'## Example from Lawson and Hanson (1974), p.171:
##'e = cbind(c(-.207,-.392,.599), c(2.558, -1.351, -1.206))
##'f = c(-1.3,-.084,.384)
##'ldp(e, f)            # Solution: 0.1268538 -0.2554018
##'
##'@usage
##'lsei(a, b, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf)
##'lsi(a, b, e=NULL, f=NULL, lower=-Inf, upper=Inf)
##'ldp(e, f)
##'qp(q, p, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf, tol=1e-15)
##' 
##'@export lsei
##'@export lsi
##'@export ldp
##'@export qp 

lsei = function(a, b, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf) {
  if(is.null(c) | length(c) == 0) return(lsi(a, b, e, f))
  if(is.vector(c)) dim(c) = c(1, length(c))
  c.qr = qr(t(c))
  L = t(qr.R(c.qr))
  k = c.qr$rank              # number of effective equality constraints
  pvt = c.qr$pivot           # pivoting for equality constraints
  k1 = 1:k
  if(nrow(c) <= ncol(c) && k == nrow(c)) {
    y1 = forwardsolve(L, d)
    at = t(qr.qty(c.qr, t(a)))
  }
  else {
    L1 = L[k1,k1,drop=FALSE]
    L2 = L[-k1,k1,drop=FALSE]
    pk1 = pvt[k1]
    d1 = d[pk1]         # effective equality constraints
    d2 = d[-pk1]        # redundant equality constraints (consistent?)
    y1 = forwardsolve(L1, d1)
    if(max(abs(L2 %*% y1 - d2)) > max(abs(diag(L1))) * 1e-14)
      stop("Inconsistent equality constraints in lsei()")
    at = t(qr.qty(c.qr, t(a)))
  }
  if(any(lower != -Inf, upper != Inf)) {
    k0 = ncol(a)
    lower = rep(lower, len=k0)
    upper = rep(upper, len=k0)
    jl = lower != -Inf
    ju = upper != Inf
    e = rbind(e, diag(1, k0)[jl,], diag(-1, k0)[ju,])
    f = c(f, lower[jl], -upper[ju])
  }
  if(is.null(e)) 
    y2 = pnnls(at[,-k1,drop=FALSE], b - at[,k1,drop=FALSE] %*% y1, ncol(at)-k)$x
  else {
    if(is.vector(e)) dim(e) = c(1, length(e))
    et = t(qr.qty(c.qr, t(e)))
    y2 = lsi(at[,-k1,drop=FALSE], b - at[,k1,drop=FALSE] %*% y1,
             et[,-k1,drop=FALSE], f - et[,k1,drop=FALSE] %*% y1)
  }
  qr.qy(c.qr, c(y1, y2))
}

# beta = c(rnorm(2), 1); beta[beta<0] = 0; beta = beta/sum(beta)
# a = matrix(rnorm(18), ncol=3); b = a %*% beta + rnorm(3,sd=.1); c = matrix(rep(1, 3), nrow=1); d = 1; e = diag(rep(1,3)); f = rep(0,3); lsei(a, b, c, d, e, f)

# # c = matrix(rnorm(6), ncol=3); d = rnorm(2); a = matrix(rnorm(24), nrow=8); b = rnorm(8); e = matrix(rnorm(12), nrow=4); f = rnorm(4); x = lsei(a, b, c, d, e, f); print(x); print(c %*% x - d); e %*% x - f

# ------------------------------------------------------- #
# Least squares solution using Householder transformation #
# ------------------------------------------------------- #



##'Least Squares Solution using Householder Transformation
##'
##'Solves the least squares problem using Householder forward triangulation
##'with column interchanges. It is an R interface to the HFTI function that is
##'described in Lawson and Hanson (1974, 1995). Its Fortran implementation is
##'public domain and is available at \url{https://www.netlib.org/lawson-hanson/}.
##'
##'Given matrix \code{a} and vector \code{b}, \code{hfti} solves the least
##'squares problem:
##'
##'\deqn{\mathrm{minimize\ \ } || a x - b ||.}{minimize || a x - b ||.}
##'
##'@param a Design matrix.
##'@param b Response vector or matrix.
##'@param tol Tolerance for determining the pseudorank.
##'@return \item{b}{first \code{krank} elements contains the solution}
##'\item{krank}{psuedo-rank} \item{rnorm}{Euclidean norm of the residual
##'vector.}
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@seealso \code{\link{lsei}}, \code{\link{nnls}}.
##'@references Lawson and Hanson (1974, 1995). Solving least squares problems.
##'Englewood Cliffs, N.J., Prentice-Hall.
##'@keywords array algebra
##'@examples
##'
##'a = matrix(rnorm(10*4), nrow=10)
##'b = a %*% c(0,1,-1,1) + rnorm(10)
##'hfti(a, b)
##'
##'@export hfti

hfti = function(a, b, tol = 1e-7) {
  if(is.vector(b)) b = as.matrix(b)
  if(!(is.matrix(a) & is.matrix(b))) stop("a or b not a matrix")
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  if(m !=  dim(b)[1]) stop("dim(a)[1] != dim(b)[1]")
  nb = as.integer(dim(b)[2])
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  krank = as.integer(0)
  rnorm = double(nb)                  # only for output
  h = g = double(n)                   # n-vector of working space
  ip = rep.int(0,n)                   # m-vector of working space
  .Fortran("hfti",a=a,m,m,n,b=b,m,nb,tol,krank=krank,rnorm=rnorm,h,g,
           ip=ip,PACKAGE="lsei")[c("b","krank","rnorm")]
}

# ---------------------------------------------------------- #
# svdrs: singular value decomposition with right side vector #
#                                                            #
# For the least squares problem                              #
#                                                            #
#                ||a x - b||                                 #
# ---------------------------------------------------------- #

svdrs = function(a, b) {
  m1 = nrow(a)
  n1 = ncol(a)
  if(m1 < n1) a = rbind(a, matrix(0, nrow=n1-m1, ncol=n1))
  mda = nrow(a)
  k = min(m1, n1)
  missing.b = FALSE
  if( missing(b) ) {b = diag( rep(1.0, mda), nrow=mda ); missing.b = TRUE}
  if( is.vector(b) ) b = as.matrix(b)
  nb = ncol(b)
  if( nrow(b) < mda ) b = rbind(b, matrix(0, nrow=mda-nrow(b), ncol=ncol(b)))
  s = double(n1)
  work = double(2*n1)
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  r = .Fortran("svdrs",a=a,mda,m1,n1,b=b,nrow(b),nb,s=s,
      work,PACKAGE="lsei")[c("s","a","b")]
  if( missing.b )
    list(d=r$s[1:k], u=t(r$b[1:min(k,nrow(b)), 1:min(m1,nb),drop=FALSE]),
         v=r$a[1:n1,1:k,drop=FALSE])
  else list(d=r$s[1:k], uTb=r$b[1:min(k,nrow(b)), 1:min(m1,nb),drop=FALSE],
            v=r$a[1:n1,1:k,drop=FALSE])
}

# x = matrix(rnorm(6), nrow=3)
# r = svdrs(x)
# r$u %*% diag(r$d) %*% t(r$v) - x
# svdrs(x, 1:3)

# Quadratic programming: x^T q x / 2 + p^T x

qp = function(q, p, c=NULL, d=NULL, e=NULL, f=NULL,
    lower=-Inf, upper=Inf, tol=1e-15) {
  eq = eigen(q)
  v2 = sqrt(eq$values[eq$values >= eq$values[1] * tol])
  kr = length(v2)
  a = t(eq$vectors[,1:kr,drop=FALSE]) * v2
  b = - colSums(eq$vectors[,1:kr,drop=FALSE] * p / rep(v2, each=length(p)))
  lsei(a, b, c, d, e, f, lower, upper)
}

# partial nonnegativity quadratic programming

pnnqp = function(q, p, k=0, sum=NULL, tol=1e-20) {
  eq = eigen(q)
  v2 = sqrt(eq$values[eq$values >= eq$values[1] * tol])
  kr = length(v2)
  a = t(eq$vectors[,1:kr,drop=FALSE]) * v2
  b = - crossprod(p, eq$vectors[,1:kr,drop=FALSE])[1,] / v2
  pnnls(a, b, k, sum)
}

