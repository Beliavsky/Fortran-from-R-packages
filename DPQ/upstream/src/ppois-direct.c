/*
 *  DIRECT computation of ppois() {cumulative Poisson distribution function}
 *  ------
 *  instead of R's code :
 *
 *    ppois(x, lambda, low, log_p) := pgamma(lambda, x + 1, 1., !low, log_p)
 */

#include <float.h> /* DBL_MIN etc */

// if(verbose) :
#include <R_ext/Print.h>

#include "DPQpkg.h"

// Called from R's  okLongDouble() :
SEXP chk_LDouble(SEXP lambda_, SEXP verbose_, SEXP tol_) {
    int verbose = asLogical(verbose_);
    if(verbose == NA_LOGICAL)
	error("'verbose' must be TRUE or FALSE but is NA");
    double lambda = asReal(lambda_);
    if(lambda < 0) error("'lambda' must be >=0");
    double tol = asReal(tol_),
	rt_lam = sqrt(lambda), eps = exp(-rt_lam);
    long double
	ldlam = (long double) lambda,
	llam  = logl(ldlam),
	f0    = expl(-ldlam),// e^{-lambda}  [in long double]
	logf0 = logl(f0),    // log(e^{-lambda}) ~= -lambda  if "things work"
	explg = expl(llam),  // e^log(ldlam) ~= lambda       if "things work"
	lg1p  = log1pl(eps), // log(1 + eps) ~= eps - eps^2/2 + .. (for small eps := e^-sqrt(lam))
	RErr_log = 1 - logf0/-ldlam,
	RErr_exp = 1 - explg/ldlam,
	RErr_lg1p= 1 - lg1p/(eps*(1-eps/2))
	;
    Rboolean
	eq_log = fabsl(RErr_log) <= tol,
	eq_exp = fabsl(RErr_exp) <= tol,
	eq_lg1p= fabsl(RErr_lg1p)<= tol
	;
    if(verbose) {
	Rprintf("lambda=%g; eps = e^-sqrt(l.) = %g  ==>  logl(ldlam)=%" PR_g_
		"; expl(-ldlam)=%" PR_g_
		";\n logl(expl(-ldlam))= %" PR_g_ "~= -ldlam? rel.err=%g: %s"
		";\n expl(logl( ldlam))= %" PR_g_ " ~= ldlam? rel.err=%g: %s"
		";\n log1pl(eps)= %" PR_g_ "~= eps(1-eps/2)?  rel.err=%g: %s"
		"\n"
 		, lambda, eps, llam, f0
		, logf0, (double)RErr_log, (eq_log ? "TRUE" : "FALSE")
		, explg, (double)RErr_exp, (eq_exp ? "TRUE" : "FALSE")
		, lg1p,  (double)RErr_lg1p,(eq_lg1p? "TRUE" : "FALSE")
	    );
    }
    return ScalarLogical(eq_log && eq_exp);
}


/** Direct computation of the cumulative Poisson distribution function  ppois()
 *
 * Really "two functions", depending on  all_from_0:
 *
 * 1)  all_from_0 = FALSE : 'x' a *vector* (of length n)
 *
 *     ppoisD(x, lambda, 0) := ppois(x, lambda, low=TRUE, log_p=FALSE)
 *                         { = pgamma(lambda, x + 1, 1., !low, log_p) }

 * 2)  all_from_0 =  TRUE :  'x' must be "scalar" (length 1)
 *
 *     ppoisD(x, lambda, 1) := ppois(0:x, lambda, low=TRUE, log_p=FALSE)
 *                                   ~~~
 */
SEXP ppoisD(SEXP X, SEXP lambda_, SEXP all_from_0, SEXP verbose_)
{
    if(!isReal(X))
	error("'x' must be a \"double\" numeric vector");
    double
	*x = REAL(X),
	lam = asReal(lambda_), jI;
    if(ISNAN(lam)) error("lambda is NA -- invalid here");
    if(lam <= 0.)  error("lambda <= 0 is invalid here"); // == 0 : just here
    int from_0 = asLogical(all_from_0);
    if(from_0 == NA_LOGICAL)
	error("'all.from.0' must be TRUE or FALSE but is NA");
    int verbose = asInteger(verbose_);
    if(verbose == NA_INTEGER || verbose < 0)
	error("'verbose' must be in {0,1,2,..} but is %d", verbose);
    R_xlen_t i, nx = XLENGTH(X), n;
    if(from_0) { // return all probabilities at 0:x, i.e., 0:(n-1)
	jI = 1 + floor(x[0] + 1e-7);
	if(jI > R_XLEN_T_MAX)
	    error(_("x (= %g) is too large here"), jI);
	n = (R_xlen_t) jI;
    } else {
	n = nx; // return probabilities  at x[i], i = 0,..,(n-1)
	jI = ceil(lam) - 1; // := arg max_j  { dpois(j, lam) } = arg max_j{lam^j / j!}
	if(jI >= INT_MAX)
	    error("ceiling(lambda) > INT_MAX is invalid here");
    }

    SEXP Prob = PROTECT(allocVector(REALSXP, n)); // the result
    double *prob = REAL(Prob);
    long double f, P,
	exp_arg = 0, // "= 0": -Wall
	ldlam = (long double) lam,
	llam  = logl(ldlam),
	f0    = expl(-ldlam); // e^{-lambda}  [in long double]
	// f0 = f_0; where  f_j := e^{-lam} lam^j / j!  for j = 0,1,...
    if (f0 == 0.L) {
	exp_arg = -ldlam;
	if(verbose)
	    REprintf("ppoisD(*, lambda=%g): expl(-ldlam)=%" PR_g_ "= 0 ==> llam=%" PR_g_
		     ", exp_arg=%" PR_g_ "\n", lam, f0, llam, exp_arg);
    } else if(verbose)
	    REprintf("ppoisD(*, lambda=%g): ldlam = %" PR_g_ ", expl(-ldlam)=%" PR_g_
		     "; llam=%" PR_g_ ", exp_arg=%" PR_g_ "\n",
		     lam, ldlam, f0, llam, exp_arg);

    for(i = 0; i < n; i++) { /* prob[i] := ppois(xi, lambda) */

	/* Compute Prob = sum_{j = 0}^xi f_j,  where f_j := exp(-L) L^j / j! {L:= lambda}

	 * Numerically smartly, would mean to always
	 * sum starting with small terms to largest, which would mean two parts :
	 *
	 *	S1 = sum_{j =  0:j_    }  f_j	[forwards];
	 *	S2 = sum_{j = xi:(j_+1)}  f_j   [backwards].
	 *
	 * However,  if(from_0) for speed reasons we will just go _forwards_
	 */
	if(from_0) {
	    if(i == 0) {
		P = f = f0; // NB: here too, f = f_i := e^{-lam} lam^i / i!
	    } else if(f > 4*LDBL_MIN) { // i >= 1
		f *= ldlam/i;
		// ==>     f == f_i := e^{-lam} lam^i / i!
		P += f; // P == sum_{m=0, i} f_m
	    } else { // i >= 1, f = f0 = 0:
		//  f := e^-lam lam^i / i! = exp(-lam + i*log(lam) - log(i!))
		//       and  log(i!) = log(i* (i-1)!) = log(i) + log((i-1)!)
		exp_arg += llam - logl((long double)i); // ... not accurate ?
		// exp_arg = -ldlam + i*llam - lgammal((long double)(i+1));
		if((f = expl(exp_arg)) > 0) {
		    P += f; // P == sum_{m=0, i} f_m
		    if(verbose)
			REprintf(" ..>> i=%lld, finally new f = expl(exp_arg = %"PR_g_
				 ") = %"PR_g_ " > 0\n", (long long)i, exp_arg, f);
		} else if(verbose >= 2)
			REprintf(" .. i=%lld, f = expl(exp_arg = %"PR_g_") = %"PR_g_"\n",
				 (long long)i, exp_arg, f);
	    }
	    prob[i] = (double) P;
	}
	else {
	    // xi = x[i] .. x[] may be out of order, so we "start from scratch"
	    double xi = floor(x[i] + 1e-7);
	    if (ISNAN(xi))
		error("x[%lld] is NA -- invalid here", (long long)i+1);
	    if (xi < 0)		{ prob[i]= 0.; break; } // incl -Inf
	    if (!R_FINITE(xi))	{ prob[i]= 1.; break; } //      +Inf
	    Rboolean sml_x = (xi <= jI);
	    int j_ = (int)(sml_x ? xi : jI); // the maximal f_j term is at j = j_ <= jI < INT_MAX
	    f = f0; // == e^{-lambda}  [in long double]
	    // = f_0;  will always be f = f_j := e^{-lam} lam^j / j!
	    if (f0 == 0.L)
		exp_arg = -ldlam;
	    long double S1 = f;
	    for(int j = 1; j <= j_; j++) { // now, f = f_{j-1}
		if(f > LDBL_MIN)
		    f *= ldlam/j; // ==> f = f_j = f_{j-1} * (lam / j)
		else { // S1 = f = f0 = 0
		    exp_arg += llam - logl((long double)j);
		    f = expl(exp_arg);
		}
		S1 += f;
	    }

	    if(sml_x) {
		prob[i] = (double) S1;
	    }
	    else { // !sml_x : xi > jI := ceil(lam) - 1  ==> start summation at xi :
		/* f := f_xi = e^-lam * lam^xi / xi! = e^-lam * lam^xi / gamma(xi+1)
		 *           = e^{-lam + xi*log(lam) - log(gamma(xi+1)) } */
		if(verbose >= 2) {
		    f = expl(-ldlam + xi*llam - lgammal((long double)(xi+1)));
		    if(f == 0L) {
			REprintf("ppoisD(x=%g, lambda=%g, expl(-ldlam)=%"PR_g_ "=0 ==> log(lam)=%"PR_g_
				 ", exp_arg=%" PR_g_ "\n",
				 xi, lam, f0, llam, exp_arg);
			xi--;
			while((f = expl(-ldlam + xi*logl(ldlam) - lgammal((long double)(xi+1)))) == 0L
			      && xi > j_+1)
			    xi--;
		    }
		} else { // not verbose
		    while((f = expl(-ldlam + xi*logl(ldlam) - lgammal((long double)(xi+1)))) == 0L
			  && xi > j_+1)
			xi--;
		}
		long double S2 = f;
		for(int j = (int)xi; j > j_+1; j--) {// backwards; now, f = f_j
		    // f := f_{j-1} = f_j * (j / lam)
		    if(f > LDBL_MIN)
			f *= j/ldlam;
		    else { // f == 0 (underflow) or  f "subnormal" -- cannot accurately update
			f = expl(-ldlam + j*llam - lgammal((long double)(j+1)));
			if(verbose >=2 && f)
			    REprintf(" .. j=%d, finally new f = expl(.) = %" PR_g_ " > 0\n", j, f);
		    }
		    S2 += f;
		}
		prob[i] = (double) (S1 + S2);
	    }
	}

    } /* i = 0..(n-1) */
    UNPROTECT(1);
    return Prob;
}
