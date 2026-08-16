#' Nonlinear optimization using augmented Lagrange method (original version)
#'
#' @param pars an numeric vector of decision variables (length n).
#' @param fun the objective function (must return a scalar).
#' @param eqfun an optional function for calculating equality constraints.
#' @param eqB a vector of the equality bounds (if eq_fn provided).
#' @param ineqfun an optional function for calculating inequality constraints.
#' @param ineqLB the lower bounds for the inequality (must be finite)
#' @param ineqUB the upper bounds for the inequalitiy (must be finite)
#' @param LB lower bounds on decision variables
#' @param UB upper bounds on decision variables
#' @param control a list of solver control parameters (see details).
#' @param ... additional arguments passed to the supplied functions (common to all functions supplied).
#' @returns An list with the following slot:
#' \describe{
#'  \item{pars}{Optimal Parameters.}
#'  \item{convergence }{Indicates whether the solver has converged (0) or not (1 or 2).}
#'  \item{values}{Vector of function values during optimization with last one the
#'    value at the optimal.}
#'  \item{lagrange}{The vector of Lagrange multipliers.}
#'  \item{hessian}{The Hessian of the augmented problem at the optimal solution.}
#'  \item{ineqx0}{The estimated optimal inequality vector of slack variables used
#'    for transforming the inequality into an equality constraint.}
#'  \item{nfuneval}{The number of function evaluations.}
#'  \item{elapsed}{Time taken to compute solution.}
#'}
#' @details
#' The optimization problem solved by \code{csolnp} is formulated as:
#' \deqn{
#' \begin{aligned}
#' \min_{x \in \mathbb{R}^n} \quad & f(x) \\
#' \text{s.t.} \quad & g(x) = b \\
#' & h_l \le h(x) \le h_u\\
#' & x_l \le x \le x_u\\
#' \end{aligned}
#' }
#'
#' where \eqn{f(x)} is the objective function, \eqn{g(x)} is the vector of equality constraints
#' with target value \eqn{b}, \eqn{h(x)} is the vector of inequality constraints bounded
#' by \eqn{h_l} and \eqn{h_u}, with parameter bounds \eqn{x_l} and \eqn{x_u}. Internally,
#' inequality constraints are converted into equality constraints using slack variables
#' and solved using an augmented Lagrangian approach.
#' The control is a list with the following options:
#'\describe{
#'  \item{rho}{This is used as a penalty weighting scaler for infeasibility in the
#'  augmented objective function. The higher its value the more the weighting to
#'  bring the solution into the feasible region (default 1). However, very high
#'  values might lead to numerical ill conditioning or significantly slow down
#'  convergence.}
#'  \item{outer.iter}{Maximum number of major (outer) iterations (default 400).}
#'  \item{inner.iter}{Maximum number of minor (inner) iterations (default 800).}
#'  \item{delta}{Relative step size in forward difference evaluation (default 1.0e-7).}
#'  \item{tol}{ Relative tolerance on feasibility and optimality (default 1e-8).}
#'  \item{trace}{The value of the objective function and the parameters is printed
#'  at every major iteration (default 1).}
#'}
#' @examples
#' {
#' # From the original paper by Y.Ye
#' # see the unit tests for more....
#' # POWELL Problem
#' fn1 = function(x)
#' {
#'     exp(x[1] * x[2] * x[3] * x[4] * x[5])
#' }
#' eqn1 = function(x){
#'     z1 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] + x[5] * x[5]
#'     z2 = x[2] * x[3] - 5 * x[4] * x[5]
#'     z3 = x[1] * x[1] * x[1] + x[2] * x[2] * x[2]
#'     return(c(z1, z2, z3))
#' }
#' x0 = c(-2, 2, 2, -1, -1)
#' }
#' powell = solnp(x0, fun = fn1, eqfun = eqn1, eqB = c(10, 0, -1))
#' @keywords optimize
#' @rdname solnp
#' @author Alexios Galanos
#' @export
#'
solnp = function(pars, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL, control = list(), ...)
{
	# start timer
	tic = Sys.time()
	xnames = names(pars)
	# get environment
	.solnpenv <- environment()
	assign("xnames", xnames, envir = .solnpenv)
	# initiate function count
	assign(".solnp_nfn", 0, envir = .solnpenv)
	assign(".solnp_errors", 0, envir = .solnpenv)

	# index of function indicators
	# [1] length of pars
	# [2] has function gradient?
	# [3] has hessian?
	# [4] has ineq?
	# [5] ineq length
	# [6] has jacobian (inequality)
	# [7] has eq?
	# [8] eq length
	# [9] has jacobian (equality)
	# [10] has upper / lower bounds
	# [11] has either lower/upper bounds or ineq


	ind = rep(0, 11)
	np = ind[1]  = length(pars)
	# lower parameter bounds - indicator
	# lpb[1]=1 means lower/upper bounds present
	# lpb[2]=1 means lower/upper bounds OR inequality bounds present

	# do parameter and LB/UB checks
	check1 = .checkpars(pars, LB, UB, .solnpenv)

	# .LB and .UB assigned

	.LB = get(".LB", envir = .solnpenv)
	.UB = get(".UB", envir = .solnpenv)


	if( !is.null(.LB) || !is.null(.UB) ) ind[10] = 1

	# do function checks and return starting value
	funv = .checkfun(pars, fun, .solnpenv, ...)
	#.solnp_fun assigned
	.solnp_fun = get(".solnp_fun", envir = .solnpenv)

	# Analytical Gradient Functionality not yet implemented in subnp function

	# gradient and hessian checks
	#if(!is.null(grad)){
	#	gradv = .checkgrad(pars, grad, .solnpenv, ...)
	#	ind[2] = 1
	#} else{
	#	.solnp_gradfun = function(pars, ...) .fdgrad(pars, fun = .solnp_fun, ...)
		ind[2] = 0
	#	gradv = .solnp_gradfun(pars, ...)
	#}
	# .solnp_gradfun(pars, ...) assigned

	.solnp_hessfun = NULL
	ind[3] = 0
	#hessv = NULL
	# .solnp_hessfun(pars, ...) assigned

	# do inequality checks and return starting values

	if(!is.null(ineqfun)){
		ineqv 	= .checkineq(pars, ineqfun, ineqLB, ineqUB, .solnpenv, ...)
		ind[4] 	= 1
		nineq 	= length(ineqLB)
		ind[5] 	= nineq

		# check for infinities/nans
		.ineqLBx = .ineqLB
		.ineqUBx = .ineqUB
		.ineqLBx[!is.finite(.ineqLB)] = -1e10
		.ineqUBx[!is.finite(.ineqUB)] =  1e10
		ineqx0 	= (.ineqLBx + .ineqUBx)/2
		#if(!is.null(ineqgrad)){
		#	ineqjacv = .cheqjacineq(pars, gradineq, .ineqUB, .ineqLB, .solnpenv, ...)
		#	ind[6] = 1
		#} else{
		# .solnp_ineqjac = function(pars, ...) .fdjac(pars, fun = .solnp_ineqfun, ...)
		ind[6] = 0
		#ineqjacv = .solnp_ineqjac(pars, ...)
		#}
	} else{
		.solnp_ineqfun = function(pars, ...) .emptyfun(pars, ...)
		# .solnp_ineqjac = function(pars, ...) .emptyjac(pars, ...)
		ineqv 	= NULL
		ind[4] 	= 0
		nineq 	= 0
		ind[5] 	= 0
		ind[6] 	= 0
		ineqx0 	= NULL
		.ineqLB = NULL
		.ineqUB = NULL
	}
	# .solnp_ineqfun and .solnp_ineqjac assigned
	# .ineqLB and .ineqUB assigned
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .solnpenv)
	.ineqLB = get(".ineqLB", envir = .solnpenv)
	.ineqUB = get(".ineqUB", envir = .solnpenv)


	# equality checks
	if(!is.null(eqfun)){
		eqv 	= .checkeq(pars, eqfun, eqB, .solnpenv, ...)
		ind[7] 	= 1
		.eqB = get(".eqB", envir = .solnpenv)
		neq 	= length(.eqB)
		ind[8] 	= neq
		#if(!is.null(eqgrad)){
		#	eqjacv = .cheqjaceq(pars, gradeq, .solnpenv, ...)
		#	ind[9] = 1
		#} else{
		#	.solnp_eqjac = function(pars, ...) .fdjac(pars, fun = .solnp_eqfun, ...)
		#	eqjacv = .solnp_eqjac(pars, ...)
			ind[9] = 0
		#}
	} else {
		eqv = NULL
		#eqjacv = NULL
		.solnp_eqfun = function(pars, ...) .emptyfun(pars, ...)
		#.solnp_eqjac = function(pars, ...) .emptyjac(pars, ...)
		ind[7] 	= 0
		neq 	= 0
		ind[8] 	= 0
		ind[9] 	= 0
	}
	# .solnp_eqfun(pars, ...) and .solnp_eqjac(pars, ...) assigned
	# .solnp_eqB assigned
	.solnp_eqfun = get(".solnp_eqfun", envir = .solnpenv)

	if( ind[ 10 ] || ind [ 4 ]) ind[ 11 ] = 1

	# parameter bounds (pb)
	pb  = rbind( cbind(.ineqLB, .ineqUB), cbind(.LB, .UB) )

	# check control list
	ctrl  = .solnpctrl( control )
	rho   = ctrl[[ 1 ]]
	# maxit = outer iterations
	maxit = ctrl[[ 2 ]]
	# minit = inner iterations
	minit = ctrl[[ 3 ]]
	delta = ctrl[[ 4 ]]
	tol   = ctrl[[ 5 ]]
	trace = ctrl[[ 6 ]]

	# total constraints (tc) = no.inequality constraints + no.equality constraints
	tc = nineq + neq

	# initialize fn value and inequalities and set to NULL those not needed
	j  = jh = funv
	tt = 0 * .ones(3, 1)

	if( tc > 0 ) {
		# lagrange multipliers (lambda)
		lambda = 0 * .ones(tc, 1)
		# constraint vector = [1:neq 1:nineq]
		constraint = c(eqv, ineqv)
		if( ind[4] ) {
			tmpv = cbind(constraint[ (neq + 1):tc ] - .ineqLB, .ineqUB - constraint[ (neq + 1):tc ] )
			testmin = apply( tmpv, 1, FUN = function( x ) min(x[ 1 ], x[ 2 ]) )
			if( all(testmin > 0) ) ineqx0 = constraint[ (neq + 1):tc ]
			constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - ineqx0
		}
		tt[ 2 ] = .vnorm(constraint)
		if( max(tt[ 2 ] - 10 * tol, nineq, na.rm = TRUE) <= 0 ) rho = 0
	} else{
		lambda = 0
	}
	# starting augmented parameter vector
	p  = c(ineqx0, pars)
	hessv  = diag(np + nineq)
	mu = np
	.solnp_iter = 0
	ob = c(funv, eqv, ineqv)

	while( .solnp_iter < maxit ){
		.solnp_iter = .solnp_iter + 1
		.subnp_ctrl = c(rho, minit, delta, tol, trace)

		# make the scale for the cost, the equality constraints, the inequality
		# constraints, and the parameters
		if( ind[7] ) {
			# [1 neq]
			vscale = c( ob[ 1 ], rep(1, neq) * max( abs(ob[ 2:(neq + 1) ]) ) )
		} else {
			vscale = 1
		}

		if( !ind[ 11 ] ) {
			vscale = c(vscale, p)
		} else {
			# [ 1 neq np]
			vscale = c(vscale, rep( 1, length.out = length(p) ) )
		}

		vscale = apply( matrix(vscale, ncol = 1), 1, FUN = function( x ) min( max( abs(x), tol ), 1/tol ) )

		res   = .subnp(pars = p, yy = lambda, ob = ob, hessv = hessv, lambda = mu, vscale = vscale,
				ctrl = .subnp_ctrl, .env = .solnpenv, ...)
		if(get(".solnp_errors", envir =  .solnpenv) == 1){
			maxit = .solnp_iter
		}
		p  = res$p
		lambda  = res$y
		hessv  = res$hessv
		mu = res$lambda
		temp = p[ (nineq + 1):(nineq + np) ]
		funv = .safefunx(temp, .solnp_fun, .env = .solnpenv, ...)
		ctmp = get(".solnp_nfn", envir =  .solnpenv)
		assign(".solnp_nfn", ctmp + 1, envir = .solnpenv)

		tempdf = cbind(temp, funv)

		if( trace ){
			.report(.solnp_iter, funv, temp)
		}

		eqv = .solnp_eqfun(temp, ...)
		ineqv = .solnp_ineqfun(temp, ...)

		ob = c(funv, eqv, ineqv)

		tt[ 1 ] = (j - ob[ 1 ]) / max(abs(ob[ 1 ]), 1)
		j = ob[ 1 ]

		if( tc > 0 ){
			constraint = ob[ 2:(tc + 1) ]

			if( ind[ 4 ] ){
				tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ],
				              pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )

				if( min(tempv) > 0 ) {
					p[ 1:nineq ] = constraint[ (neq + 1):tc ]
				}

				constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - p[ 1:nineq ]
			}

			tt[ 3 ] = .vnorm(constraint)

			if( tt[ 3 ] < 10 * tol ) {
				rho = 0
				mu  = min(mu, tol)
			}

			if( tt[ 3 ] < 5 * tt[ 2 ]) {
				rho = rho/5
			}

			if( tt[ 3 ] > 10 * tt[ 2 ]) {
				rho = 5 * max( rho, sqrt(tol) )
			}

			if( max( c( tol + tt[ 1 ], tt[ 2 ] - tt[ 3 ] ) ) <= 0 ) {
				lambda = 0
				hessv = diag( diag ( hessv ) )
			}

			tt[ 2 ] = tt[ 3 ]
		}

		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			maxit = .solnp_iter
		}

		jh = c(jh, j)
	}

	if( ind[ 4 ] ) {
		ineqx0 = p[ 1:nineq ]
	}

	p = p[ (nineq + 1):(nineq + np) ]

	if(get(".solnp_errors", envir =  .solnpenv) == 1){
		convergence = 2
		if( trace ) cat( paste( "\nsolnp--> Solution not reliable....Problem Inverting Hessian.\n", sep="" ) )
	} else{
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			convergence = 0
			if( trace ) cat( paste( "\nsolnp--> Completed in ", .solnp_iter, " iterations\n", sep="" ) )
		} else{
			convergence = 1
			if( trace ) cat( paste( "\nsolnp--> Exiting after maximum number of iterations\n",
							"Tolerance not achieved\n", sep="" ) )
		}
	}
	# end timer
	ctmp = get(".solnp_nfn", envir =  .solnpenv)
	toc = Sys.time() - tic
	names(p) = xnames
	ans = list(pars = p, convergence = convergence, values = as.numeric(jh), lagrange = lambda,
			hessian = hessv, ineqx0 = ineqx0, nfuneval = ctmp, outer.iter = .solnp_iter,
			elapsed = toc, vscale = vscale)
	return( ans )
}
