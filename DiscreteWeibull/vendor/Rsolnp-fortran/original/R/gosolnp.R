#' Random Initialization and Multiple Restarts of the solnp solver.
#'
#' When the objective function is non-smooth or has many local minima, it is hard
#' to judge the optimality of the solution, and this usually depends critically on
#' the starting parameters. This function enables the generation of a set of
#' randomly chosen parameters from which to initialize multiple restarts of the
#' solver (see note for details).
#'
#' @param pars The starting parameter vector. This is not required unless the fixed option is
#' also used.
#' @param fixed The numeric index which indicates those parameters which should stay fixed
#' instead of being randomly generated.
#' @param fun The main function which takes as first argument the parameter vector and returns
#' a single value.
#' @param eqfun (Optional) The equality constraint function returning the vector of evaluated
#' equality constraints.
#' @param eqB (Optional) The equality constraints.
#' @param ineqfun (Optional) The inequality constraint function returning the vector of evaluated
#' inequality constraints.
#' @param ineqLB (Optional) The lower bound of the inequality constraints.
#' @param ineqUB (Optional) The upper bound of the inequality constraints.
#' @param LB The lower bound on the parameters. This is not optional in this function.
#' @param UB The upper bound on the parameters. This is not optional in this function.
#' @param control (Optional) The control list of optimization parameters. The \code{eval.type}
#' option in this control list denotes whether to evaluate the function as is and
#' exclude inequality violations in the final ranking (default, value = 1),
#' else whether to evaluate a penalty barrier function comprised of the objective
#' and all constraints (value = 2). See \code{solnp} function documentation for
#' details of the remaining control options.
#' @param distr A numeric vector of length equal to the number of parameters, indicating the
#' choice of distribution to use for the random parameter  generation. Choices are
#' uniform (1), truncated normal (2), and normal (3).
#' @param distr.opt If any choice in \code{distr} was anything other than uniform (1), this is a
#' list equal to the length of the parameters with sub-components for the mean and
#' sd, which are required in the truncated normal and normal distributions.
#' @param n.restarts The number of solver restarts required.
#' @param n.sim The number of random parameters to generate for every restart of the solver.
#' Note that there will always be significant rejections if inequality bounds are
#' present. Also, this choice should also be motivated by the width of the upper
#' and lower bounds.
#' @param cluster If you want to make use of parallel functionality, initialize and pass a cluster
#' object from the parallel package (see details), and remember to terminate it!
#' @param rseed (Optional) A seed to initiate the random number generator, else system time will
#' be used.
#' @param ... (Optional) Additional parameters passed to the main, equality or inequality
#' functions
#' @details
#' Given a set of lower and upper bounds, the function generates, for those
#' parameters not set as fixed, random values from one of the 3 chosen
#' distributions. Depending on the \code{eval.type} option of the \code{control}
#' argument, the function is either directly evaluated for those points not
#' violating any inequality constraints, or indirectly via a penalty barrier
#' function jointly comprising the objective and constraints. The resulting values
#' are then sorted, and the best N (N = random.restart) parameter vectors
#' (corresponding to the best N objective function values) chosen in order to
#' initialize the solver. Since version 1.14, it is up to the user to prepare and
#' pass a cluster object from the parallel package for use with gosolnp, after
#' which the parLapply function is used. If your function makes use of additional
#' packages, or functions, then make sure to export them via the \code{clusterExport}
#' function of the parallel package. Additional arguments passed to the solver via the
#' \dots option are evaluated and exported by gosolnp to the cluster.
#' @return A list containing the following values:
#' \item{pars}{Optimal Parameters.}
#' \item{convergence }{Indicates whether the solver has converged (0) or not (1).}
#' \item{values}{Vector of function values during optimization with last one the
#' value at the optimal.}
#' \item{lagrange}{The vector of Lagrange multipliers.}
#' \item{hessian}{The Hessian at the optimal solution.}
#' \item{ineqx0}{The estimated optimal inequality vector of slack variables used
#' for transforming the inequality into an equality constraint.}
#' \item{nfuneval}{The number of function evaluations.}
#' \item{elapsed}{Time taken to compute solution.}
#' \item{start.pars}{The parameter vector used to start the solver}
#' @references
#' Y.Ye, \emph{Interior algorithms for linear, quadratic, and linearly constrained
#' non linear programming}, PhD Thesis, Department of EES Stanford University,
#' Stanford CA.\cr
#' Hu, X. and Shonkwiler, R. and Spruill, M.C. \emph{Random Restarts in Global
#' Optimization}, 1994, Georgia Institute of technology, Atlanta.
#' @author Alexios Galanos and Stefan Theussl\cr
#' Y.Ye (original matlab version of solnp)
#' @note
#' The choice of which distribution to use for randomly sampling the parameter
#' space should be driven by the user's knowledge of the problem and confidence or
#' lack thereof of the parameter distribution. The uniform distribution indicates
#' a lack of confidence in the location or dispersion of the parameter, while the
#' truncated normal indicates a more confident choice in both the location and
#' dispersion. On the other hand, the normal indicates perhaps a lack of knowledge
#' in the upper or lower bounds, but some confidence in the location and dispersion
#' of the parameter. In using choices (2) and (3) for \code{distr},
#' the \code{distr.opt} list must be supplied with \code{mean} and \code{sd} as
#' subcomponents for those parameters not using the uniform (the examples section
#' hopefully clarifies the usage).
#' @examples
#' \dontrun{
#' # [Example 1]
#' # Distributions of Electrons on a Sphere Problem:
#' # Given n electrons, find the equilibrium state distribution (of minimal Coulomb
#' # potential) of the electrons positioned on a conducting sphere. This model is
#' # from the COPS benchmarking suite. See http://www-unix.mcs.anl.gov/~more/cops/.
#' gofn = function(dat, n)
#' {
#'
#' 	x = dat[1:n]
#' 	y = dat[(n+1):(2*n)]
#' 	z = dat[(2*n+1):(3*n)]
#' 	ii = matrix(1:n, ncol = n, nrow = n, byrow = TRUE)
#' 	jj = matrix(1:n, ncol = n, nrow = n)
#' 	ij = which(ii<jj, arr.ind = TRUE)
#' 	i = ij[,1]
#' 	j = ij[,2]
#' 	#  Coulomb potential
#' 	potential = sum(1.0/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2))
#' 	potential
#' }
#'
#' goeqfn = function(dat, n)
#' {
#' 	x = dat[1:n]
#' 	y = dat[(n+1):(2*n)]
#' 	z = dat[(2*n+1):(3*n)]
#' 	apply(cbind(x^2, y^2, z^2), 1, "sum")
#' }
#'
#' n = 25
#' LB = rep(-1, 3*n)
#' UB = rep(1,  3*n)
#' eqB = rep(1, n)
#' ans = gosolnp(pars  = NULL, fixed = NULL, fun = gofn, eqfun = goeqfn, eqB = eqB,
#' LB = LB, UB = UB, control = list(outer.iter = 100, trace = 1),
#' distr = rep(1, length(LB)), distr.opt = list(), n.restarts = 2, n.sim = 20000,
#' rseed = 443, n = 25)
#' # should get a function value around 243.813
#'
#' # [Example 2]
#' # Parallel functionality for solving the Upper to Lower CVaR problem (not properly
#' # formulated...for illustration purposes only).
#'
#' mu =c(1.607464e-04, 1.686867e-04, 3.057877e-04, 1.149289e-04, 7.956294e-05)
#' sigma = c(0.02307198,0.02307127,0.01953382,0.02414608,0.02736053)
#' R = matrix(c(1, 0.408, 0.356, 0.347, 0.378,  0.408, 1, 0.385, 0.565, 0.578, 0.356,
#' 0.385, 1, 0.315, 0.332, 0.347, 0.565, 0.315, 1, 0.662, 0.378, 0.578,
#' 0.332, 0.662, 1), 5,5, byrow=TRUE)
#' # Generate Random deviates from the multivariate Student distribution
#' set.seed(1101)
#' v = sqrt(rchisq(10000, 5)/5)
#' S = chol(R)
#' S = matrix(rnorm(10000 * 5), 10000) %*% S
#' ret = S/v
#' RT = as.matrix(t(apply(ret, 1, FUN = function(x) x*sigma+mu)))
#' # setup the functions
#' .VaR = function(x, alpha = 0.05)
#' {
#' 	VaR = quantile(x, probs = alpha, type = 1)
#' 	VaR
#' }
#'
#' .CVaR = function(x, alpha = 0.05)
#' {
#' 	VaR = .VaR(x, alpha)
#' 	X = as.vector(x[, 1])
#' 	CVaR = VaR - 0.5 * mean(((VaR-X) + abs(VaR-X))) / alpha
#' 	CVaR
#' }
#' .fn1 = function(x,ret)
#' {
#' 	port=ret%*%x
#' 	obj=-.CVaR(-port)/.CVaR(port)
#' 	return(obj)
#' }
#'
#' # abs(sum) of weights ==1
#' .eqn1  = function(x,ret)
#' {
#' 	sum(abs(x))
#' }
#'
#' LB=rep(0,5)
#' UB=rep(1,5)
#' pars=rep(1/5,5)
#' ctrl = list(delta = 1e-10, tol = 1e-8, trace = 0)
#' cl = makePSOCKcluster(2)
#' # export the auxilliary functions which are used and cannot be seen by gosolnp
#' clusterExport(cl, c(".CVaR", ".VaR"))
#' ans = gosolnp(pars, fun = .fn1, eqfun = .eqn1, eqB = 1, LB = LB, UB = UB,
#' n.restarts = 2, n.sim=500, cluster = cl, ret = RT)
#' ans
#' # don't forget to stop the cluster!
#' stopCluster(cl)
#' }
#' @export
gosolnp = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
		ineqUB = NULL, LB = NULL, UB = NULL, control = list(), distr = rep(1, length(LB)), distr.opt = list(),
		n.restarts = 1, n.sim = 20000, cluster = NULL, rseed = NULL, ...)
{
    # allowed distributions:
    # 1: uniform (no confidence in the location of the parameter...somewhere in LB-UB space)
    # 2: truncnorm (high confidence in the location of the parameter)
    # 3: normal (Uncertainty in Lower and Upper bounds, but some idea about the dispersion about the location)
    # ...
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(control$trace)) trace = FALSE else trace = as.logical(control$trace)
	if(is.null(control$eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(control$eval.type),2))
	if(parmethod == 0) parmethod = 1
	control$eval.type = NULL
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\ngosolnp-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\ngosolnp-->error: the function requires upper parameter bounds\n", call. = FALSE)
	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\ngosolnp-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = 1:n

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\ngosolnp-->error: fixed indices out of bounds\n", call. = FALSE)
	}
	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
				eval(parse(text=paste(xl[i],"=list(...)[[i]]",sep="")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		clusterEvalQ(cluster, require(Rsolnp))
	}
	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, trace = trace, rseed = rseed,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, rseed = rseed, trace = trace,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...))

	gosolnp_rndpars = gosolnp_rndpars[,1:n, drop = FALSE]
	# initiate solver restarts
	if( trace ) cat("\ngosolnp-->Starting Solver\n")
	solution = vector(mode = "list", length = n.restarts)
	if( !is.null(cluster) )
	{
		clusterExport(cluster, c("gosolnp_rndpars"), envir = environment())
		solution = parLapply(cluster, as.list(1:n.restarts), fun = function(i) {
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	} else {
		solution = lapply(as.list(1:n.restarts), FUN = function(i){
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	}
	if(n.restarts>1){
		best = sapply(solution, FUN = function(x) if(x$convergence!=0) NA else x$values[length(x$values)])
		if(all(is.na(best)))
			stop("\ngosolnp-->Could not find a feasible starting point...exiting\n", call. = FALSE)
		nb = which(best == min(best, na.rm = TRUE))[1]
		solution = solution[[nb]]
		if( trace ) cat("\ngosolnp-->Done!\n")
		solution$start.pars = gosolnp_rndpars[nb,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	} else{
		solution = solution[[1]]
		solution$start.pars = gosolnp_rndpars[1,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	}
	return(solution)
}


#' Generates and returns a set of starting parameters by sampling the parameter
#' space based on the evaluation of the function and constraints.
#'
#' A simple penalty barrier function is formed which is then evaluated at randomly
#' sampled points based on the upper and lower parameter bounds
#' (when \code{eval.type} = 2), else the objective function directly for values not
#' violating any inequality constraints (when \code{eval.type} = 1). The sampled
#' points can be generated from the uniform, normal or truncated normal
#' distributions.
#'
#' @param pars The starting parameter vector. This is not required unless the fixed option is
#' also used.
#' @param fixed The numeric index which indicates those parameters which should stay fixed
#' instead of being randomly generated.
#' @param fun The main function which takes as first argument the parameter vector and returns
#' a single value.
#' @param eqfun (Optional) The equality constraint function returning the vector of evaluated
#' equality constraints.
#' @param eqB (Optional) The equality constraints.
#' @param ineqfun (Optional) The inequality constraint function returning the vector of evaluated
#' inequality constraints.
#' @param ineqLB (Optional) The lower bound of the inequality constraints.
#' @param ineqUB (Optional) The upper bound of the inequality constraints.
#' @param LB The lower bound on the parameters. This is not optional in this function.
#' @param UB The upper bound on the parameters. This is not optional in this function.
#' @param distr A numeric vector of length equal to the number of parameters, indicating the
#' choice of distribution to use for the random parameter generation. Choices are
#' uniform (1), truncated normal (2), and normal (3).
#' @param distr.opt If any choice in \code{distr} was anything other than uniform (1), this is a
#' list equal to the length of the parameters with sub-components for the mean and
#' sd, which are required in the truncated normal and normal distributions.
#' @param bestN The best N (less than or equal to n.sim) set of parameters to return.
#' @param n.sim The number of random parameter sets to generate.
#' @param cluster If you want to make use of parallel functionality, initialize and pass a cluster
#' object from the parallel package (see details), and remember to terminate it!
#' @param rseed (Optional) A seed to initiate the random number generator, else system time will
#' be used.
#' @param eval.type Either 1 (default) for the direction evaluation of the function (excluding
#' inequality constraint violations) or 2 for the penalty barrier method.
#' @param trace (logical) Whether to display the progress of the function evaluation.
#' @param ... (Optional) Additional parameters passed to the main, equality or inequality
#' functions
#' @details
#' Given a set of lower and upper bounds, the function generates, for those
#' parameters not set as fixed, random values from one of the 3 chosen
#' distributions. For simple functions with only inequality constraints, the direct
#' method (\code{eval.type} = 1) might work better. For more complex setups with
#' both equality and inequality constraints the penalty barrier method
#' (\code{eval.type} = 2)might be a better choice.
#' @return A matrix of dimension bestN x (no.parameters + 1). The last column is the
#' evaluated function value.
#' @author Alexios Galanos and Stefan Theussl\cr
#' @note
#' The choice of which distribution to use for randomly sampling the parameter
#' space should be driven by the user's knowledge of the problem and confidence or
#' lack thereof of the parameter distribution. The uniform distribution indicates a
#' lack of confidence in the location or dispersion of the parameter, while the
#' truncated normal indicates a more confident choice in both the location and
#' dispersion. On the other hand, the normal indicates perhaps a lack
#' of knowledge in the upper or lower bounds, but some confidence in the location
#' and dispersion of the parameter. In using choices (2) and (3) for \code{distr},
#' the \code{distr.opt} list must be supplied with \code{mean} and \code{sd} as
#' subcomponents for those parameters not using the uniform.
#' @examples
#' \dontrun{
#' library(Rsolnp)
#' library(parallel)
#' # Windows
#' cl = makePSOCKcluster(2)
#' # Linux:
#' # makeForkCluster(nnodes = getOption("mc.cores", 2L), ...)
#'
#' gofn = function(dat, n)
#' {
#'
#' 	x = dat[1:n]
#' 	y = dat[(n+1):(2*n)]
#' 	z = dat[(2*n+1):(3*n)]
#' 	ii = matrix(1:n, ncol = n, nrow = n, byrow = TRUE)
#' 	jj = matrix(1:n, ncol = n, nrow = n)
#' 	ij = which(ii<jj, arr.ind = TRUE)
#' 	i = ij[,1]
#' 	j = ij[,2]
#' 	#  Coulomb potential
#' 	potential = sum(1.0/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2))
#' 	potential
#' }
#'
#' goeqfn = function(dat, n)
#' {
#' 	x = dat[1:n]
#' 	y = dat[(n+1):(2*n)]
#' 	z = dat[(2*n+1):(3*n)]
#' 	apply(cbind(x^2, y^2, z^2), 1, "sum")
#' }
#' n = 25
#' LB  = rep(-1, 3*n)
#' UB  = rep( 1, 3*n)
#' eqB = rep( 1,   n)
#'
#' sp = startpars(pars = NULL, fixed = NULL, fun = gofn , eqfun = goeqfn,
#' eqB = eqB, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = LB, UB = UB,
#' distr = rep(1, length(LB)), distr.opt = list(), n.sim = 2000,
#' cluster = cl, rseed = 100, bestN = 15, eval.type = 2, n = 25)
#' #stop cluster
#' stopCluster(cl)
#' # the last column is the value of the evaluated function (here it is the barrier
#' # function since eval.type = 2)
#' print(round(apply(sp, 2, "mean"), 3))
#' # remember to remove the last column
#' ans = solnp(pars=sp[1,-76],fun = gofn , eqfun = goeqfn , eqB = eqB, ineqfun = NULL,
#' ineqLB = NULL, ineqUB = NULL, LB = LB, UB = UB, n = 25)
#' # should get a value of around 243.8162
#' }
#' @export
startpars = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL,
		ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL,
		distr = rep(1, length(LB)), distr.opt = list(), n.sim = 20000, cluster = NULL,
		rseed = NULL, bestN = 15, eval.type = 1, trace = FALSE, ...)
{
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(eval.type),2))
	if(parmethod == 0) parmethod = 1
	eval.type = NULL
	#trace = FALSE
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\nstartpars-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\nstartpars-->error: the function requires upper parameter bounds\n", call. = FALSE)

	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\nstartpars-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = seq_len(n)

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\nstartpars-->error: fixed indices out of bounds\n", call. = FALSE)
	}

	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}

	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
			  eval(parse(text = paste(xl[i], "=list(...)", "[[" , i, "]]", sep = "")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		if( !is.null(names(list(...))) ) parallel::clusterExport(cluster, names(list(...)), envir = environment())
		clusterEvalQ(cluster, require(Rsolnp))
	}

	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...))
	return(gosolnp_rndpars)
}


.randpars = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB,
		LB, UB, distr, distr.opt, n.restarts, n.sim, trace = TRUE, rseed,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	m = length(nf)
	set.seed(rseed)
	for(i in 1:m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}

	if( trace ) cat("ok!\n")

	if(!is.null(ineqfun)){
		if( trace ) cat("\nExcluding Inequality Violations...\n")
		ineqv = matrix(NA, ncol = length(ineqLB), nrow = n.restarts*n.sim)
		# ineqv = t(apply(rndpars, 1, FUN = function(x) ineqfun(x)))
		if(length(ineqLB) == 1){
			ineqv = apply(gosolnp_rndpars, 1, FUN = function(x){
						names(x) = gosolnp_parnames
						ineqfun(x, ...)} )
			lbviol = sum(ineqv<ineqLB)
			ubviol = sum(ineqv>ineqUB)
			if( lbviol > 0 | ubviol > 0 ){
				vidx = c(which(ineqv<ineqLB), which(ineqv>ineqUB))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)
			} else{
				vidx = 0
				lvx = 0
			}
		} else{
			ineqv = t(apply(gosolnp_rndpars, 1, FUN = function(x){
								names(x) = gosolnp_parnames
								ineqfun(x, ...)} ))

			# check lower and upper violations
			lbviol = apply(ineqv, 1, FUN = function(x) sum(any(x<ineqLB)))
			ubviol = apply(ineqv, 1, FUN = function(x) sum(any(x>ineqUB)))
			if( any(lbviol > 0) | any(ubviol > 0) ){
				vidx = c(which(lbviol>0), which(ubviol>0))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)

			} else{
				vidx = 0
				lvx = 0
			}
		}
		if( trace ) cat(paste("\n...Excluded ", lvx, "/",n.restarts*n.sim, " Random Sequences\n", sep = ""))
	}
	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", ".safefun"), envir = environment())
		evfun = parLapply(cluster, as.list(1:nx), fun = function(i){
					.safefun(gosolnp_rndpars[i, ], fun, gosolnp_parnames, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x) .safefun(x, fun, gosolnp_parnames, ...))
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	if( trace ) cat("ok!\n")
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}

# form a barrier function before passing the parameters
.randpars2 = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB, LB,
		UB, distr, distr.opt, n.restarts, n.sim, rseed, trace = TRUE,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_idx = "a"
	gosolnp_R = NULL
	if(!is.null(ineqfun) && is.null(eqfun) ){
		gosolnp_idx = "b"
		gosolnp_R = 100
	}
	if( is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "c"
		gosolnp_R = 100
	}
	if(!is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "d"
		gosolnp_R = c(100,100)
	}
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	gosolnp_m = length(nf)
	set.seed(rseed)
	for(i in 1:gosolnp_m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}
	if( trace ) cat("ok!\n")
	# Barrier Function
	pclfn = function(x){
		z=x
		z[x<=0] = 0
		z[x>0] = (0.9+z[x>0])^2
		z
	}
	.lagrfun = function(pars, m, idx, fun, eqfun = NULL, eqB = 0, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, ...)
	{
		fn = switch(idx,
				"a" = fun(pars[1:m], ...),
				"b" = fun(pars[1:m], ...) + pars[m+1]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ),
				"c" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]),
				"d" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]) + pars[m+2]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ) )
		return(fn)
	}

	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", "gosolnp_m", "gosolnp_idx",
						"gosolnp_R"), envir = environment())
		clusterExport(cluster, c("pclfn", ".lagrfun"), envir = environment())
		evfun = parallel::parLapply(cluster, as.list(1:nx), fun = function(i){
					.lagrfun(c(gosolnp_rndpars[i,], gosolnp_R), gosolnp_m,
							gosolnp_idx, fun, eqfun, eqB, ineqfun, ineqLB,
							ineqUB, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x){
					.lagrfun(c(x,gosolnp_R), gosolnp_m, gosolnp_idx, fun, eqfun,
							eqB, ineqfun, ineqLB, ineqUB, ...)})
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	#distmat = dist(evfun, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}


.distr1 = function(LB, UB, n)
{
	runif(n, min = LB, max = UB)
}

.distr2 = function(LB, UB, n, mean, sd)
{
	rtruncnorm(n, a = as.double(LB), b = as.double(UB), mean = as.double(mean), sd = as.double(sd))
}

.distr3 = function(n, mean, sd)
{
	rnorm(n, mean = mean, sd = sd)
}

.safefun = function(pars, fun, gosolnp_parnames, ...){
	# gosolnp_parnames = get("gosolnp_parnames", envir = .env)
	names(pars) = gosolnp_parnames
	v  = fun(pars, ...)
	if(is.na(v) | !is.finite(v) | is.nan(v)) {
		warning(paste("\ngosolnp-->warning: ", v , " detected in function call...check your function\n", sep = ""), immediate. = FALSE)
		v = 1e24
	}
	v
}
