#include <Rcpp.h>
#include <R_ext/RS.h>
#include <R_ext/Arith.h>
#include <R_ext/Rdynload.h>
#include <limits>

using namespace Rcpp;
using namespace std;

/// Wrapper for the objective function.  It is initialized to R's "c" function
static Function cf("c");

/** 
 * Fortran callable objective function evaluation.
 * 
 * @param n size of parameter vector
 * @param x parameter vector
 * @param ip print flag
 * 
 * @return objective function evaluation 
 */
extern "C"
double F77_NAME(calfun)(int const *n, double const x[], int const *ip) {
    Environment rho(cf.environment());
    IntegerVector cc(rho.get(".feval."));
    int nn = *n;
    cc[0]++;			// increment func eval count

    if (count_if(x, x + nn, R_finite) < nn)
	throw range_error("non-finite x values not allowed in calfun");

    SEXP pp = PROTECT(::Rf_allocVector(REALSXP, nn));
    copy(x, x + nn, REAL(pp));
    double f = ::Rf_asReal(::Rf_eval(PROTECT(::Rf_lang2(as<SEXP>(cf), pp)), as<SEXP>(rho)));
    UNPROTECT(2);

#if 0
    double f;
    try {
       f  = as<double>(cf(pp)); // evaluate objective
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
#endif
    if (!R_finite(f)) f = numeric_limits<double>::max();

    if (*ip == 3) {		// print eval info when very verbose
	Rprintf("%3d:%#14.8g:", cc[0], f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
    if (*ip > 3 && cc[0] % *ip == 0) { // print eval info every *ip if *ip >3
      Rprintf("%3d:%#14.8g:", cc[0], f);
      for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
      Rprintf("\n");
    }
    
    return f;
}

/** 
 * Construct the named and classed list to return from the optimizer
 * 
 * @param par parameter vector
 * @param cnm class name
 * 
 * @return an Rcpp::List object
 */
static SEXP rval(NumericVector par, string cnm, int ierr = 0) {
    Environment rho(cf.environment());
    SEXP feval = rho.get(".feval.");
    StringVector cl(2);
    cl[0] = cnm;
    cl[1] = "minqa";
    
    double f = ::Rf_asReal(::Rf_eval(PROTECT(::Rf_lang2(as<SEXP>(cf),
							as<SEXP>(par))),
				     as<SEXP>(rho)));
    UNPROTECT(1);

    List rr = List::create(_["par"] = par,
			   _["fval"] = f,
			   _["feval"] = feval,
			   _["ierr"] = ierr);
    rr.attr("class") = cl;
    return rr;
}

/// Declaration of Powell's bobyqa
extern "C" 
void F77_NAME(bobyqa)(const int *n, const int *npt, double X[],
		      const double xl[], const double xu[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun, double w[],
		      int *ierr);

/// Interface for bobyqa
extern "C"
SEXP bobyqa_cpp(SEXP parp, SEXP xlp, SEXP xup, SEXP ccp, SEXP fnp) {
    try {
	NumericVector par(parp), xl(xlp), xu(xup);
	Environment cc(ccp);
	cf = Function(fnp);	// install the objective function
	double rb = as<double>(cc.get("rhobeg")),
	    re = as<double>(cc.get("rhoend"));
	int ierr = 0,
	    ip = as<int>(cc.get("iprint")),
	    mxf = as<int>(cc.get("maxfun")),
	    n = par.size(), np = as<int>(cc.get("npt"));
	vector<double> w((np + 5) * (np + n) + (3 * n * (n + 5))/2);
	NumericVector pp = clone(par); // ensure that bobyqa doesn't modify the R object
	F77_NAME(bobyqa)(&n, &np, pp.begin(), xl.begin(), xu.begin(),
			 &rb, &re, &ip, &mxf, &w[0], &ierr);
	return rval(pp, "bobyqa", ierr);
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}

extern "C" 
void F77_NAME(uobyqa)(const int *n, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun,
		      double w[], int *ierr);

extern "C"
SEXP uobyqa_cpp(SEXP parp, SEXP ccp, SEXP fnp) {
    try {
	NumericVector par(parp);
	Environment cc(ccp);
	cf = Function(fnp);
	double rb = as<double>(cc.get("rhobeg")),
	    re = as<double>(cc.get("rhoend"));
	int ierr = 0, ip = as<int>(cc.get("iprint")),
	    mxf = as<int>(cc.get("maxfun")), n = par.size();
	Environment rho(cf.environment());
	vector<double>
	    w((n*(42+n*(23+n*(8+n))) + max(2*n*n + 4, 18*n)) / 4);
	NumericVector pp = clone(par); // ensure that uobyqa doesn't modify the R object
	
	F77_NAME(uobyqa)(&n, pp.begin(), &rb, &re, &ip, &mxf, &w[0], &ierr);
	return rval(pp, "uobyqa", ierr);
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
    
extern "C" 
void F77_NAME(newuoa)(const int *n, const int *npt, double X[],
		      const double *rhobeg, const double *rhoend,
		      const int *iprint, const int *maxfun,
		      double w[], int *ierr);

extern "C"
SEXP newuoa_cpp(SEXP parp, SEXP ccp, SEXP fnp) {
    try {
	NumericVector par(parp);
	Environment cc(ccp);
	cf = Function(fnp);
	double rb = as<double>(cc.get("rhobeg")),
	    re = as<double>(cc.get("rhoend"));
	int ierr = 0, ip = as<int>(cc.get("iprint")),
	    mxf = as<int>(cc.get("maxfun")),
	    n = par.size(), np = as<int>(cc.get("npt"));
	vector<double> w((np+13)*(np+n)+(3*n*(n+3))/2);
	NumericVector pp = clone(par); // ensure that newuoa doesn't modify the R object
	
	F77_NAME(newuoa)(&n, &np, pp.begin(), &rb, &re, &ip, &mxf, &w[0], &ierr);
	return rval(pp, "newuoa", ierr);
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}

/// Assorted error messages.
extern "C"
void F77_NAME(minqer)(const int *msgno) {
BEGIN_RCPP
    const char *msg = (char*)NULL;
    switch(*msgno) {
    case 10:
    case 101:
	msg = "NPT is not in the required interval";
	break;
    case 20:
	msg = "one of the differences XU(I)-XL(I) is less than 2*RHOBEG";
	break;
    case 320:
	msg = "bobyqa detected too much cancellation in denominator";
	break;
    case 390:
	msg = "maximum number of function evaluations exceeded";
	break;
    case 430:
    case 3701:
    case 2101:
	msg = "a trust region step failed to reduce q";
	break;
    default:
	throw range_error("minqer message number");
    }
    throw runtime_error(msg);
VOID_END_RCPP
}

/// Iteration output when rho changes and iprint >= 2
extern "C"
void F77_NAME(minqit)(const int *iprint, const double *rho, const int *nf,
		      const double *fopt, const int *n, const double xbase[],
		      const double xopt[]) {
    if (*iprint >= 2) {
	Rprintf("rho: %#8.2g eval: %3d fn: %#12g par:", *rho, *nf, *fopt);
	for(int i = 0; i < *n; i++) Rprintf("%#8g ", xbase[i] + xopt[i]);
	Rprintf("\n");
    }
}

/// Output at return (do we really need this - why not use the print method?)
extern "C" void
F77_NAME(minqir)(const int *iprint, const double *f, const int *nf,
		 const int *n, const double x[]) {
    if (*iprint > 0) {
	Rprintf("At return\n");
	Rprintf("eval: %3d fn: %#14.8g par:", *nf, *f);
	for (int i = 0; i < *n; i++) Rprintf(" %#8g", x[i]);
	Rprintf("\n");
    }
}

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(bobyqa_cpp, 5),
    CALLDEF(uobyqa_cpp, 3),
    CALLDEF(newuoa_cpp, 3),
    {NULL, NULL, 0}
};

/// Initializer for the package.  Registers the symbols for .Call.
extern "C" void R_init_minqa(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
