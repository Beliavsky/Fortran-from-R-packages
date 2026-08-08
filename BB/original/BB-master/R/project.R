
projectLinear <- function(par, A, b, meq) {
# A projection function to incorporate linear equalities and inequalities in nonlinear optimization using `spg'
#
# The inequalities are defined such that:  A %*% x - b > 0 

n <- length(par)
if (meq > 0 | any(b - c(A %*% par) > 0)){
   par <- par +  quadprog::solve.QP(Dmat=diag(1,n), dvec=rep(0, n), Amat=t(A),
		      bvec = b - c(A %*% par), meq=meq, factorized=TRUE)$solution
}
par
}

