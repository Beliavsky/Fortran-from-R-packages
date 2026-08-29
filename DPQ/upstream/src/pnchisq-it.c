/*
 *  Computes the noncentral chi-squared distribution function with
 *  positive real degrees of freedom f and nonnegative noncentrality
 *  parameter theta
 */
#include "DPQpkg.h" // before stdio (for MINGW_...)

#include <float.h> // for DBL_*

/* --> below:  Pnchisq_it()  called from R
 *             ------------
 * FIXME:  1) additionally vectorize !?

 * DONE:  o added 'verbose' argument, instead of  '#ifdef ...'
 */

/**  This is a version of R (1.9.0)'s src/nmath/pnchisq.c  [R 1.9.0 @ 2004-04-12 => svn r28809]
 *                                    ~~~~~~~~~~~~~~~~~~~ */
double pnchisq_it(double x, double f, double theta,
		  double errmax, double reltol, int itrmax, Rboolean verbose,
		  int *i_, int *n_terms, double *terms)
{
    double ans, lam, u, v, x2, f2, t, term, f_x_2n, f_2n, lt;
    double lu = -1., l_lam = -1., l_x = -1., bound = -1.; /* initialized for -Wall */
    int n;
    Rboolean lamSml, tSml, is_r, is_b, is_it = FALSE;

    static const double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;
    /*= -708.3964 for IEEE double precision */

    i_[0] = 0;
    *n_terms = 0;

    if (x <= 0.)	return 0.;
    if(!R_FINITE(x))	return 1.;

    if(verbose)
	REprintf("pnchisq(x=%g, f=%g, theta=%g): ",x,f,theta);
    lam = .5 * theta;
    lamSml = (-lam < _dbl_min_exp); // <==> lam > 1022*log(2) ~= 708.3964
    if(lamSml) {
	/* MATHLIB_ERROR(
	   "non centrality parameter (= %g) too large for current algorithm",
	   theta) */
        u = 0;
        lu = -lam;/* == ln(u) */
        l_lam = log(lam);
    } else {
	u = exp(-lam);
    }

    /* evaluate the first term */
    v = u;
    x2 = .5 * x;
    f2 = .5 * f;
    f_x_2n = f - x;

    if(verbose)
	REprintf("-- v=exp(-th/2)=%g, x/2= %g, f/2= %g\n",v,x2,f2);

    if(f2 * DBL_EPSILON > 0.125 && /* very large f and x ~= f: probably needs */
       fabs(t = x2 - f2) <         /* other algorithm anyway */
       sqrt(DBL_EPSILON) * f2) {
	/* evade cancellation error */
	/* t = exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2*M_PI*(f2 + 1));*/
        lt = (1 - t)*(2 - t/(f2 + 1)) - 0.5 * log(2*M_PI*(f2 + 1));
	if(verbose)
	    REprintf(" (case I) ==> ");
    }
    else {
	/* Usual case 2: careful not to overflow .. : */
	lt = f2*log(x2) -x2 - lgammafn(f2 + 1);
    }
    if(verbose)
	REprintf(" lt= %g", lt);

    tSml = (lt < _dbl_min_exp);
    if(tSml) {           //  5 : R only warns for x > E[X] + 5*sigma
	if (x > f + theta +  3* sqrt( 2*(f + 2*theta))) {
 	    /* x > E[X] + 3* sigma(X) */
	    warning("x > E[X] + 3*sigma(X) -- result may not be good");
	    // return 1.;
	    /* better than 0 --- but definitely "FIXME" */
	}
	// else {
	l_x = log(x);
	ans = term = t = 0.;
	i_[0]++; /* initial  0-term that we leave away */
	// }
    }
    else {
	t = exp(lt);
	if(verbose)
	    REprintf(", t=exp(lt)= %g\n", t);
	ans = term = v * t;
    }
    terms[0] = term;

    for (n = 1, f_2n = f + 2., f_x_2n += 2.;  ; n++, f_2n += 2, f_x_2n += 2) {
	if(verbose)
	    REprintf("\n _OL_: n=%d",n);
	/* f_2n    === f + 2*n
	 * f_x_2n  === f - x + 2*n   > 0  <==> (f+2n)  >   x */
	if (f_x_2n > 0) {
	    /* find the error bound and check for convergence */

	    bound = t * x / f_x_2n;
	    if(verbose)
		REprintf("\n L10: n=%d; term= %g; bound= %g",n,term,bound);
	    is_r = is_it = FALSE;
	    /* convergence only if BOTH absolute and relative error < 'bnd' */
	    if (((is_b = (bound <= errmax)) &&
                 (is_r = (term <= reltol * ans))) || (is_it = (n > itrmax)))
            {
		if(verbose)
		    REprintf("BREAK n=%d %s; bound= %g %s, rel.err= %g %s\n",
			     n, (is_it ? "> itrmax" : ""),
			     bound, (is_b ? "<= errmax" : ""),
			     term/ans, (is_r ? "<= reltol" : ""));
		break; /* out completely */
            }
	}
	else if (n > itrmax) { // (not used ?  just in case)
	    REprintf("series not converged, n=%d > itrmax (while f_x_2n = %g <= 0)\n",
		     n, f_x_2n);
	    break;
	}

	// evaluate the next term of the expansion and then the partial sum :

        if(lamSml) {
            lu += l_lam - log(n); /* u = u* lam / n */
            if(lu >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
		if(verbose)
		    REprintf(" n=%d; nomore underflow in u = exp(lu) ==> change\n",
			     n);
                v = u = exp(lu); /* the first non-0 'u' */
                lamSml = FALSE;
            }
        } else {
	    u *= lam / n;
	    v += u;
	}
	if(tSml) {
            lt += l_x - log(f_2n);/* t <- t * (x / f2n) */
            if(lt >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
		if(verbose)
		    REprintf("  n=%d; nomore underflow in t = exp(lt) ==> change\n",
			     n);
                t = exp(lt); /* the first non-0 't' */
                tSml = FALSE;
            }
        } else {
	    t *= x / f_2n;
	}
        if(!lamSml && !tSml) {
	    term = v * t;
	    ans += term;
	} else {
	    i_[0]++; /* another initial 0-term */
	}
	terms[n] = term;

    } /* for(n ...) */

    if (is_it) {
	REprintf("pnchisq(x=%g, ..): not converged in %d iter.",
		 x, itrmax);
    }
    if(verbose)
	REprintf("\n == L_End: n=%d; term= %g; bound=%g\n",n,term,bound);
    *n_terms = n;
    return (ans);
}

/** Called from R's pnchisqIT()   as  .C(C_Pnchisq_it, ...)  --->> ../R/pnchisq.R
 */
// TODO: vectorize in x
void Pnchisq_it(double *x, double *f, double *theta,
		double *errmax, double *reltol, int *itrmax, int *verbose,
		int *i_0, int *n_terms, double *terms, double *prob)
{
    *prob = pnchisq_it(*x, *f, *theta, *errmax, *reltol, *itrmax, *verbose,
		       i_0, n_terms, terms);
}

//---------- current version of R's  pnchisq_raw() -----------------------------------------
// >>> ~/R/D/r-devel/R/src/nmath/pnchisq.c

static const double _dbl_min_exp = M_LN2 * DBL_MIN_EXP;
/*= -708.3964 for IEEE double precision */

double pnchisq_rawR(double x, double f, double theta /* = ncp */,
		    double cutoff_ncp, Rboolean small_logspace, int it_simple,
		    double errmax, double reltol, double epsS, int itrmax, int verbose,
		    Rboolean lower_tail, Rboolean log_p, LDOUBLE *sum, LDOUBLE *sum2);

// called in a loop over x from Pnchisq_R() which is called from R's pnchisqRC() :
double pnchisqR(double x, double df, double ncp, Rboolean lower_tail, Rboolean log_p,
		double cutoff_ncp, Rboolean small_logspace,
		//        80     || small_logspace was (lower_tail && f > 0. && log(x) < .....)
		Rboolean no_2nd_call, int it_simple,
		//         FALSE            110
		double errmax, double reltol, double epsS, int itrmax, int verbose)
		//     1e-12,  8*DBL_EPSILON,        1e-15    1000000,
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(df) || ISNAN(ncp))
	return x + df + ncp;
    if (!R_FINITE(df) || !R_FINITE(ncp))
	ML_WARN_return_NAN;
#endif
    if (df < 0. || ncp < 0.) ML_WARN_return_NAN;

    LDOUBLE sum, sum2;
    double
    ans = pnchisq_rawR(x, df, ncp, cutoff_ncp, small_logspace, it_simple,
		       errmax, reltol, epsS, itrmax,
		       verbose, lower_tail, log_p, &sum, &sum2);

    if (x <= 0. || x == ML_POSINF)
	return ans; // because it's perfect

    if(ncp >= cutoff_ncp) {
	if(lower_tail) {
	    ans = fmin2(ans, R_D__1);  /* e.g., pchisq(555, 1.01, ncp = 80) */
	} else { /* !lower_tail */
	    /* since we computed the other tail cancellation is likely */
	    // FIXME: There are cases where  ans == 0. if(!log_p) is perfect
	    if(ans < (log_p ? (-10. * M_LN10) : 1e-10)) {
		if(verbose)
		    REprintf(" ans := pnch.raw(*, ncp >= cutoff, <upper tail>)=%g \"too small\" -> precision warning\n",
			     ans);
		ML_WARNING(ME_PRECISION, "pnchisq");
	    }
	    if(!log_p && ans < 0.) ans = 0.;  /* Precaution PR#7099 */
	    // else if(log_p && ISNAN(ans)) ans = ML_NEGINF;
	}
    }
    /* MM: the following "hack" from c51179 (<--> PR#14216, by Jerry Lewis)
     * -- is "kind of ok" ... but potentially suboptimal: we do  log1p(- p(*, <other tail>, log=FALSE)),
     *    but that  p(*, log=FALSE) may already be an exp(.) or even expm1(..)
     *   <---> "in principle"  this check should happen there, not here  */
    if (no_2nd_call || !log_p || ans < -1e-8)
	return ans;
    else { // log_p  &&  ans >= -1e-8
	// prob. = exp(ans) is near one: we can do better using the other tail
	if(verbose)
	    REprintf("   pnchisq_raw(*, log_p): ans=%g => 2nd call, log1p(- <other tail no log>)\n", ans);
	ans = pnchisq_rawR(x, df, ncp, cutoff_ncp, small_logspace, it_simple,
			   errmax, reltol, epsS, itrmax,
			   verbose, !lower_tail, FALSE, &sum, &sum2);
	return log1p(-ans);
    }
}

double pnchisq_rawR(double x, double f, double theta /* = ncp */,
		    double cutoff_ncp, Rboolean small_logspace, int it_simple,
		    double errmax, double reltol, double epsS, int itrmax, int verbose,
		    Rboolean lower_tail, Rboolean log_p, LDOUBLE *sum, LDOUBLE *sum2)
{
    double lam, x2, f2, term, bound = -1., f_x_2n, f_2n;
    double l_lam = -1., l_x = -1.; /* initialized for -Wall */
    int n;
    Rboolean lamSml, tSml, is_r, is_b;
    LDOUBLE ans, u, v, t, lt, lu =-1;

    if (x <= 0.) {
	if(x == 0. && f == 0.) { // chi^2_0(.) has point mass at zero
#define _L  (-0.5 * theta) // = -lambda
	    return lower_tail ? R_D_exp(_L) : (log_p ? R_Log1_Exp(_L) : -expm1(_L));
	}
	/* x < 0  or {x==0, f > 0} */
	return R_DT_0;
    }
    if(!R_FINITE(x))	return R_DT_1;

    /* This is principally for use from qnchisq */
    R_CheckUserInterrupt();

    // cutoff_ncp was '80' hardcoded  //  it_simple was '110' hardcoded
    if(theta < cutoff_ncp) { /* use 110 for Inf, as ppois(110, 80/2, lower.tail=FALSE) is 2e-20 */
	int i;
	// Have  pgamma(x,s) < x^s / Gamma(s+1) (< and ~= for small x)
	// ==> pchisq(x, f) = pgamma(x, f/2, 2) = pgamma(x/2, f/2)
	//                  <  (x/2)^(f/2) / Gamma(f/2+1) < eps
	// <==>  f/2 * log(x/2) - log(Gamma(f/2+1)) < log(eps) ( ~= -708.3964 )
	// <==>        log(x/2) < 2/f*(log(Gamma(f/2+1)) + log(eps))
	// <==> log(x) < log(2) + 2/f*(log(Gamma(f/2+1)) + log(eps))
	if(small_logspace) { /* was
        if(lower_tail && f > 0. &&
	    log(x) < M_LN2 + 2/f*(lgamma(f/2. + 1) + _dbl_min_exp)) {
			     ==> is default via small.ncp.logspace.R2015(..) in ../R/pnchisq.R
			     */
	    // all  pchisq(x, f+2*i, lower_tail, FALSE), i=0,...,110 would underflow to 0.
	    // ==> work in log scale
	    double lambda = 0.5 * theta; // < cutoff_ncp/2
	    double pr = -lambda, log_lam = log(lambda);
	    *sum = *sum2 = (LDOUBLE) ML_NEGINF;
	    /* we need to renormalize here: the result could be very close to 1 */
	    if(verbose >= 2) REprintf("  logspace iterations: showing ' i, sum2;' :\n  ");
	    for(i = 0; i < it_simple;  pr += log_lam - log(++i)) {
		*sum2 = (LDOUBLE) logspace_add((double)*sum2, pr);
		*sum  = (LDOUBLE) logspace_add((double)*sum , pr + pchisq(x, f+2*i, lower_tail, TRUE));
		if(verbose >= 2)
		    REprintf(" %d: %" PR_g_ ";", i, *sum2);
		if (*sum2 >= -epsS) /*<=> EXP(sum2) >= 1-epsS */ break;
	    }
	    ans = *sum - *sum2;
	    if(verbose >= 2) REprintf(" final i=%d\n", i);
	    if(verbose)
		REprintf("pnchisq(x=%g, f=%g, th.=%g); th. < cutoff_ncp=%g, logspace: i=%d,"
			 " ans=(sum=%" PR_g_ ")-(sum2=%" PR_g_ ")\n",
			 x,f,theta, cutoff_ncp, i, *sum, *sum2);
	    if (i >= it_simple)
		MATHLIB_WARNING2(_("pnchisq(x=%g, ..): I: not converged in %d simple iterations"),
				 x, it_simple);
	    return (double)(lower_tail ? R_D_EXP(ans) : R_D_1EXP(ans));
	}
	else {
	    LDOUBLE lambda = 0.5 * theta; // < cutoff_ncp/2  ( = 40 )
	    LDOUBLE pr = EXP(-lambda); // does this need a feature test?
	    *sum = *sum2 = (LDOUBLE) 0;
	    /* we need to renormalize here: the result could be very close to 1 */
	    if(verbose >= 2) REprintf(" sum iterations: showing ' i, sum2;' :\n  ");
	    for(i = 0; i < it_simple;  pr *= lambda/++i) {
		// pr == exp(-lambda) lambda^i / i!  ==  dpois(i, lambda)
		*sum2 += pr;
		// pchisq(*, i, *) is  strictly decreasing to 0 for lower_tail=TRUE
		//                 and strictly increasing to 1 for lower_tail=FALSE
		*sum += pr * pchisq(x, f+2*i, lower_tail, FALSE);
		if(verbose >= 2) REprintf(" %d, %" PR_g_ ";", i, *sum2);
		if (*sum2 >= 1-epsS) break;
	    }
	    if(verbose >= 2) REprintf(" final i=%d\n", i);
	    ans = *sum / *sum2;
	    if(verbose)
		REprintf("pnchisq(x=%g, f=%g, theta=%g); theta < cutoff_ncp=%g: i=%d, "
			 "ans=(sum=%" PR_g_ ")/(sum2=%" PR_g_ ")\n",
			 x,f,theta, cutoff_ncp, i, *sum, *sum2);
	    if (i >= it_simple)
		MATHLIB_WARNING2(_("pnchisq(x=%g, ..): II: not converged in %d simple iterations"),
				 x, it_simple);
	    return (double) (log_p ? LOG(ans) : ans);
	}
    } // if(theta < cutoff_ncp)

    // else: theta == ncp >= cutoff_ncp --------------------------------------------
    if(verbose)
	REprintf("pnchisq(x=%g, f=%g, theta=%g >= cutoff_ncp = %g): ",
		 x,f,theta, cutoff_ncp);

    // Series expansion --- FIXME: log_p=TRUE, lower_tail=FALSE only applied at end ==> underflow

    lam = .5 * theta; // = lambda = ncp/2
    lamSml = (-lam < _dbl_min_exp);
    if(lamSml) {
	/* MATHLIB_ERROR(
	   "non centrality parameter (= %g) too large for current algorithm",
	   theta) */
        u = 0;
        lu = -lam;/* == ln(u) */
        l_lam = log(lam);
    } else {
	u = exp(-lam);
    }

    /* evaluate the first term */
    v = u;
    x2 = .5 * x;
    f2 = .5 * f;
    f_x_2n = f - x;

    if(verbose)
	REprintf("-- v=exp(-th/2)=%" PR_g_ ", x/2= %g, f/2= %g\n",v,x2,f2);

    if(f2 * DBL_EPSILON > 0.125 && /* very large f and x ~= f: probably needs */
       FABS(t = x2 - f2) <         /* another algorithm anyway */
       sqrt(DBL_EPSILON) * f2) {
	/* evade cancellation error */
	/* t = exp((1 - t)*(2 - t/(f2 + 1))) / sqrt(2*M_PI*(f2 + 1));*/
        lt = (1 - t)*(2 - t/(f2 + 1)) - M_LN_SQRT_2PI - 0.5 * log(f2 + 1);
	if(verbose)
	    REprintf(" (case I) ==> ");
    }
    else {
	/* Usual case 2: careful not to overflow .. : */
	lt = f2*log(x2) -x2 - lgammafn(f2 + 1);
    }
    if(verbose) REprintf(" lt= %" PR_g_ "", lt);
    tSml = (lt < _dbl_min_exp);
    if(tSml) {
	if(verbose)
	    REprintf(" is very small\n");

	if (x > f + theta +  5* sqrt( 2*(f + 2*theta))) {
	    /* x > E[X] + 5* sigma(X) */
	    return R_DT_1; /* FIXME: could be more accurate than 0. */
	} /* else */
	l_x = log(x);
	ans = term = 0.; t = 0;
    }
    else {
	t = EXP(lt);
	if(verbose)
	    REprintf(", t=exp(lt)= %" PR_g_ "\n", t);
	ans = term = (double) (v * t);
    }

    for (n = 1, f_2n = f + 2., f_x_2n += 2.; n <= itrmax ; n++, f_2n += 2, f_x_2n += 2) { // --------
	if(verbose >= 2) {
	    if(n % 1000 == 0)
		REprintf("\n _OL_: n=%d,  f_x_2n = %g", n, f_x_2n);
	    else
		REprintf(n % 100 == 0 ? ".\n" : ".");
	}
#ifndef MATHLIB_STANDALONE
	if(n % 1000 == 0) R_CheckUserInterrupt();
#endif
	/* f_2n    === f + 2*n
	 * f_x_2n  === f - x + 2*n   > 0  <==> (f+2n) > x  ==> have positive error bound
	 * <==> n > (x - f)/2 : when too large, MUST use different "algorithm": asymptotic approx ! */
	if (f_x_2n > 0) {
	    /* find the error bound and check for convergence */

	    bound = (double) (t * x / f_x_2n);
	    if(verbose >= 2 && n % 1000 == 0)
		REprintf("\n L10: n=%d; term, ans = %g, %" PR_g_ "; bound= %g",
			 n, term, ans, bound);
	    is_r = FALSE;
	    /* convergence only if BOTH absolute and relative error < 'bnd' */
	    if (((is_b = (bound <= errmax)) &&
                 (is_r = (term <= reltol * ans))))
            {
		if(verbose)
		    REprintf("BREAK from for(n = 1 ..): n=%d; bound= %g %s; term=%g, rel.err= %g %s\n",
			     n,
			     bound, (is_b ? "<= errmax" : ""), term,
			     (double)(term/ans), (is_r ? "<= reltol" : ""));
		break; /* out completely */
            }
	}

	/* evaluate the next term of the */
	/* expansion and then the partial sum */

        if(lamSml) {
            lu += l_lam - log(n); /* u = u* lam / n */
            if(lu >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
		if(verbose)
		    REprintf(" n=%d; nomore underflow in u = exp(lu) ==> change\n",
			     n);
                v = u = EXP(lu); /* the first non-0 'u' */
                lamSml = FALSE;
            }
        } else {
	    u *= lam / n;
	    v += u;
	}
	if(tSml) {
            lt += l_x - log(f_2n);/* t <- t * (x / f2n) */
            if(lt >= _dbl_min_exp) {
		/* no underflow anymore ==> change regime */
		if(verbose)
		    REprintf("  n=%d; nomore underflow in t = exp(lt) ==> change\n", n);
                t = EXP(lt); /* the first non-0 't' */
                tSml = FALSE;
            }
        } else {
	    t *= x / f_2n;
	}
        if(!lamSml && !tSml) {
	    term = (double) (v * t);
	    ans += term;
	}

    }// for(n ...) -----------------------

    if (n > itrmax) {
	MATHLIB_WARNING4(_("pnchisq(x=%g, f=%g, theta=%g, ..): not converged in %d iter."),
			 x, f, theta, itrmax);
    }
    double dans = (double) ans;
    if(verbose)
	// workaround MinGW Windows bug
#ifdef _WIN32
	REprintf("\n == L_End: n=%d; term= %g; bound=%g: [dbl]ans=%g\n",
		 n, term, bound, dans);
#else // working printf(<long double)
	REprintf("\n == L_End: n=%d; term= %g; bound=%g: ans=%" PR_g_ "\n",
		 n, term, bound, ans);
#endif
    return R_DT_val(dans);
}


/** Called from R's pnchisqRC()   as  .Call(C_R_Pnchisq_R, ...)  --->> ../R/pnchisq.R <<
 *
 * with slightly different name conventions:
 *
 *                        Name in R */
SEXP Pnchisq_R(SEXP x_, //	  q
	       SEXP f_, //	 df
	       SEXP theta_, //	ncp
	       SEXP lower_tail_, SEXP log_p_,
	       SEXP no_2nd_call_,
	       SEXP cutoff_ncp_, SEXP small_ncp_logspace_, SEXP it_simple_,
	       SEXP errmax_, SEXP reltol_, SEXP epsS_, SEXP itrmax_, SEXP verbose_)
{
    // vectorized in  (x,f,th)
    if(!isReal(x_) || !isReal(f_) || !isReal(theta_))
	error("'x', 'df', and 'ncp' must be \"double\" numeric vectors");
    if(!isLogical(small_ncp_logspace_)) error("'small.ncp.logspace' must be logical");
    double // of length 1 :
	errmax = asReal(errmax_),
	reltol = asReal(reltol_),
	epsS   = asReal(epsS_),
	cutoff_ncp = asReal(cutoff_ncp_);
    if(ISNAN(errmax) || errmax < 0) error("'errmax' must be numeric, >=0");
    if(ISNAN(reltol) || reltol < 0) error("'reltol' must be numeric, >=0");
    if(ISNAN(epsS)   || epsS  <= 0) error("'epsS' must be numeric, > 0");
    if(ISNAN(cutoff_ncp) || cutoff_ncp < 0) error("'cutoff_ncp' must be numeric, >=0");
    if(asLogical(lower_tail_) == NA_LOGICAL || asLogical(log_p_) == NA_LOGICAL)
	error("'lower.tail' and 'log.p' must be TRUE or FALSE (not NA!)");
    if(asLogical(no_2nd_call_) == NA_LOGICAL) error("'no2nd.call', must be TRUE or FALSE");
    Rboolean
	no_2nd_call = asLogical(no_2nd_call_),
	lower_tail = asLogical(lower_tail_),
	log_p      = asLogical(log_p_);
    int verbose   = asInteger(verbose_), // 0, 1, 2 ..
	itrmax    = asInteger(itrmax_),
	it_simple = asInteger(it_simple_);
    if(verbose == NA_INTEGER || verbose < 0)
	error("'verbose' must be TRUE, FALSE, or integer 0, 1,..");
    if(itrmax == NA_INTEGER || itrmax < 0)
	error("'itrmax' must be a non-negative integer");
    if(it_simple == NA_INTEGER || it_simple < 0)
	error("'it_simple' must be a non-negative integer");
    R_xlen_t
	n_x  = XLENGTH(x_),
	n_f  = XLENGTH(f_),
	n_th = XLENGTH(theta_),
	n_slg= XLENGTH(small_ncp_logspace_),
	n = n_x,
	i = (n_f >= n_th) ? n_f : n_th;
    if(!n_x || !n_f || !n_th || !n_slg) return allocVector(REALSXP, 0); // length 0
    // otherwise, recycle to common length n :
    if(n < i)     n = i;     // ==> n = max(length(x), length(f), length(theta))
    if(n < n_slg) n = n_slg; // ==> n = max(length(x), length(f), length(th.), length(sml_ncp.))
    if(verbose) {
	REprintf("Pnchisq_R(x, f, th, ... lower.tail=%d, log.p=%d, cut_ncp=%g, it_simple=%d,\n"
		 "  errmax=%g, reltol=%g, epsS=%g, itrmax=%d, verbose=%d)\n"
		 "  --> n:= max(length(.),..) = %lld\n",
		 lower_tail, log_p, cutoff_ncp, it_simple,
		 errmax, reltol, epsS, itrmax, verbose, (long long)n);
    }
    SEXP r_ = PROTECT(allocVector(REALSXP, n)); // result
    double *x = REAL(x_), *f = REAL(f_), *theta = REAL(theta_),
	*r = REAL(r_);
    int *small_ncp_logspace = LOGICAL(small_ncp_logspace_);
    for(i=0; i < n; i++) {
	r[i] = pnchisqR(x[i % n_x], f[i % n_f], theta[i % n_th],
			(Rboolean)lower_tail, (Rboolean)log_p,
			cutoff_ncp, (Rboolean)small_ncp_logspace[i % n_slg],
			no_2nd_call, it_simple,
			errmax, reltol, epsS, itrmax, verbose);
    }
    UNPROTECT(1);
    return(r_);
}
