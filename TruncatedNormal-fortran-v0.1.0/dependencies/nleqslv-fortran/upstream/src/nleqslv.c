#include <R.h>
#include <Rinternals.h>

extern double getjacond(void);

typedef struct opt_struct {
    SEXP x;
    SEXP fcall; /* function */
    SEXP jcall; /* jacobian */
    SEXP env;   /* environment in which to evaluate calls */
    SEXP names; /* names of starting values */
    int  dsub;  /* number of subdiagonals if jacobian is banded */
    int  dsuper;/* number of superdiagonals if jacobian is banded */
} opt_struct, *OptStruct;

OptStruct OS;

/*
 * set to 1 if main function called
 * used to check for a recursive call
 * which is not allowed (for the time being)
 */
static int activeflag = 0;

void deactivatenleq(void);  /* called on.exit; see .R code */

void fcnval(double *xc, double *fc, int *n, int *flag);
void fcnjac(double *rjac, int *ldr, double *x, int *n);

void F77_NAME(nwnleq)(double *x, int *n, double *scalex, int *maxit, int *jacflg, double *xtol,
                      double *ftol, double *btol, double *cndtol,
                      int *method, int *global, int *xscalm, double *stepmx, double *delta, double *sigma,
                      double *rjac, int *ldr,
                      double *rwork, int *lrwork, double *rcdwrk, int *icdwrk, double *qrwork, int *qrwsiz,
                      void (*fcnjac)(double *rjac, int *ldr, double *x, int *n),
                      void (*fcnval)(double *xc, double *fc, int *n, int *flag),
                      int *outopt, double *xp, double *fp, double *gp, int *njcnt, int *nfcnt, int *iter,
                      int *termcd);

void F77_NAME(liqsiz)(int *n, int *wrksiz); /* returns size of optimal QR work memory */

void deactivatenleq(void)
{
    activeflag = 0;
}

static SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    return elmt;
}

static double *real_vector(int n)
{
    return (double*) R_alloc(n, sizeof(double));
}

static int *int_vector(int n)
{
    return (int*) R_alloc(n, sizeof(int));
}

static void trace_header(int method, int global, int xscalm, double sigma,
                         double delta, double stepmx, double ftol, double xtol, double btol, double cndtol)
{   /* print header for iteration tracing output for documentation purposes
     */

    Rprintf("  Algorithm parameters\n  --------------------\n");
    Rprintf("  Method: %s", method == 1 ? "Broyden" : "Newton");
    Rprintf("  Global strategy: ");
    switch(global)
    {
        case 0: Rprintf("none\n"); break;
        case 1: Rprintf("cubic linesearch\n"); break;
        case 2: Rprintf("quadratic linesearch\n"); break;
        case 3: Rprintf("geometric linesearch (reduction = %g)\n", sigma); break;
        case 4: Rprintf("double dogleg (initial trust region = %g)\n", delta); break;
        case 5: Rprintf("single dogleg (initial trust region = %g)\n", delta); break;
        case 6: Rprintf("More/Hebden/Lev/Marquardt (initial trust region = %g)\n", delta); break;
        default: error("Internal: invalid global value in trace_header\n");
    }

    Rprintf("  Maximum stepsize = %g\n", stepmx <= 0.0 ? DBL_MAX : stepmx);
    Rprintf("  Scaling: %s\n", xscalm == 0 ? "fixed" : "automatic");

    Rprintf("  ftol = %g xtol = %g btol = %g cndtol = %g\n\n", ftol,xtol,btol,cndtol);
    Rprintf("  Iteration report\n  ----------------\n");
}

static char *fcn_message(char *msg, size_t mbufsiz, int termcd)
{
    switch(termcd)
    {
        case 1: snprintf(msg, mbufsiz, "Function criterion near zero"); break;
        case 2: snprintf(msg, mbufsiz, "x-values within tolerance 'xtol'"); break;
        case 3: snprintf(msg, mbufsiz, "No better point found (algorithm has stalled)"); break;
        case 4: snprintf(msg, mbufsiz, "Iteration limit exceeded"); break;
        case 5: snprintf(msg, mbufsiz, "Jacobian is too ill-conditioned (1/condition=%7.1e) (see allowSingular option)",getjacond()); break;
        case 6: snprintf(msg, mbufsiz, "Jacobian is singular (1/condition=%7.1e) (see allowSingular option)",getjacond()); break;
        case 7: snprintf(msg, mbufsiz, "Jacobian is completely unusable (all zero entries?)"); break;
        case -10: snprintf(msg, mbufsiz, "User supplied Jacobian most likely incorrect"); break;
        default: snprintf(msg, mbufsiz, "'termcd' == %d should *NEVER* be returned! Please report bug where directed in the DESCRIPTION file.", termcd);
    }
    return msg;
}

#define max(a,b) ((a)>(b) ? (a):(b))
#define min(a,b) ((a)<(b) ? (a):(b))

static int  findcol(int row, int n, int k)
{
    int j, col = 0;

    for(j=k; j <= n; j += OS->dsub+OS->dsuper+1)
    {
        if( row >= max(j-OS->dsuper,1) && row <= min(j+OS->dsub,n) )
            col = j;
        break;
    }
    return col;
}

/*
 * interface to user supplied R function
 * (*flag) == 0 when function is called for function values only
 * (*flag) >  0 jacobian column number when function is called for numeric jacobian
 *         > *n (*flag-*n) is strip number for banded evaluation
 */

void fcnval(double *x, double *fc, int *n, int *flag)
{
    int i;
    SEXP sexp_fvec;

    for (i = 0; i < *n; i++) {
         if (!R_FINITE(x[i]))
             error("non-finite value for `x[%d]` supplied to function\n",i+1);
         REAL(OS->x)[i] = x[i];
    }

    SETCADR(OS->fcall, OS->x);
    PROTECT(sexp_fvec = eval(OS->fcall, OS->env));
    if(!isReal(sexp_fvec)) error("function must return a numeric vector");
    if(LENGTH(sexp_fvec) != *n) error("function return should be a vector of length %d but is of length %d\n",
                                       LENGTH(sexp_fvec), *n);
    for (i = 0; i < *n; i++) {
        fc[i] = REAL(sexp_fvec)[i];
        if( !R_FINITE(fc[i]) ) {
            fc[i] = sqrt(DBL_MAX / (double)(*n)); /* should force backtracking */
            if( *flag ) {
                if( *flag <= *n )
                    error("non-finite value(s) detected in jacobian (row=%d,col=%d)",i+1,*flag);
                else
                    error("non-finite value(s) detected in banded jacobian (row=%d,col=%d)",i+1,
                           findcol(i+1,*n,*flag-*n));
            }
        }
    }

    UNPROTECT(1);
}

void FCNJACDUM(double *rjac, int *ldr, double *x, int *n)
{
    error("FCNJACDUM should not have been called");
}

/*
 * interface to user supplied jacobian function
 */

void fcnjac(double *rjac, int *ldr, double *x, int *n)
{
    int i, j;
    SEXP sexp_fjac;
    SEXP jdims;

    for (i = 0; i < *n; i++) {
         if (!R_FINITE(x[i]))
             error("non-finite value for `x[%d]` supplied to jacobian function\n",i+1);
         REAL(OS->x)[i] = x[i];
    }

    SETCADR(OS->jcall, OS->x);
    PROTECT(sexp_fjac = eval(OS->jcall, OS->env));
    PROTECT(jdims = getAttrib(sexp_fjac,R_DimSymbol));

    /* test jacobian of function returning scalar.
     * Scalar jacobian allowed if *n==1
     */
    if(isReal(sexp_fjac) && IS_SCALAR(sexp_fjac,REALSXP) && *n == 1 )
        ;
    /* test for numerical matrix with correct dimensions */
    else if (!isReal(sexp_fjac) || !isMatrix(sexp_fjac) || INTEGER(jdims)[0]!=*n || INTEGER(jdims)[1]!=*n)
        error("The jacobian function must return a numerical matrix of dimension (%d,%d).",*n,*n);

    for (j = 0; j < *n; j++)
        for (i = 0; i < *n; i++) {
            if( !R_FINITE(REAL(sexp_fjac)[(*n)*j + i]) )
                error("non-finite value(s) returned by jacobian (row=%d,col=%d)",i+1,j+1);
            rjac[(*ldr)*j + i] = REAL(sexp_fjac)[(*n)*j + i];
        }

    UNPROTECT(2);
}

SEXP nleqslv(SEXP xstart, SEXP fn, SEXP jac, SEXP rmethod, SEXP rglobal, SEXP rxscalm,
             SEXP rjacobian, SEXP control, SEXP rho)
{
    double  *x, *rwork, *rcdwrk, *qrwork, *rjac;
    double  *xp, *fp, *gp, *scalex;
    double  *pjac;
    int     *icdwrk;
    const char *z;

    SEXP    eval_test;
    SEXP    sexp_x, sexp_diag, sexp_fvec, sexp_info, sexp_message, sexp_nfcnt, sexp_njcnt, sexp_iter;
    SEXP    sexp_jac;
    SEXP    out, out_names;
    SEXP    xnames;

    char    message[256];

    int     i, j, n, njcnt, nfcnt, iter, termcd, lrwork, qrwsiz, lrjac, ldr;
    int     maxit, method, global, xscalm;
    int     jactype, jacflg[4], dsub, dsuper;
    int     outopt[3];
    double  xtol, ftol, btol, stepmx, delta, sigma, cndtol;

    if( activeflag )
        error("Recursive call of nleqslv not possible");
    ++activeflag;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));

    if( isReal(xstart) )
        PROTECT(OS->x = duplicate(xstart)); // 1 protected
    else if(isInteger(xstart) || isLogical(xstart) )
        PROTECT(OS->x = coerceVector(xstart,REALSXP)); // 1 protected
    else
        error("Argument 'x' cannot be converted to numeric!");

    OS->names = getAttrib(xstart, R_NamesSymbol);

    /* paranoid check */
    n = length(OS->x);
    if( n < 1 ) error("(Internal) length OS->x should be >= 1");

    for (i = 0; i < n; i++)
        if( !R_FINITE(REAL(OS->x)[i]) )
            error("'x' contains a non-finite value at index=%d\n",i+1);

    if (!isFunction(fn)) error("fn is not a function!");
    PROTECT(OS->fcall = lang2(fn, OS->x)); // 2 protected

    if (!isEnvironment(rho)) error("rho is not an environment!");
    OS->env = rho;

    PROTECT(eval_test = eval(OS->fcall, OS->env)); // 3 protected
    if (!isReal(eval_test))
        error("evaluation of fn function returns non-numeric vector!");
    i = length(eval_test);
    if( i != n )
        error("Length of fn result <> length of x!");

    for (i = 0; i < n; i++)
        if( !R_FINITE(REAL(eval_test)[i]) )
            error("initial value of fn function contains non-finite values (starting at index=%d)\n"
                  "  Check initial x and/or correctness of function",i+1);

    UNPROTECT(1); // eval-test, 2 protected

    z = CHAR(STRING_ELT(rmethod, 0));
    if( strcmp(z,"Broyden") == 0 )
        method = 1;
    else
        method = 0;

    xtol    = asReal(getListElement(control, "xtol"));
    ftol    = asReal(getListElement(control, "ftol"));
    btol    = asReal(getListElement(control, "btol"));
    sigma   = asReal(getListElement(control, "sigma"));
    stepmx  = asReal(getListElement(control, "stepmax"));
    delta   = asReal(getListElement(control, "delta"));
    cndtol  = asReal(getListElement(control, "cndtol"));

    if(!R_FINITE(xtol)) error("'xtol' is not a valid finite number");
    if(!R_FINITE(ftol)) error("'ftol' is not a valid finite number");
    if(!R_FINITE(btol)) error("'btol' is not a valid finite number");
    if(!R_FINITE(sigma)) error("'sigma' is not a valid finite number");
    if(!R_FINITE(stepmx)) error("'stepmax' is not a valid finite number");
    if(!R_FINITE(delta)) error("'delta' is not a valid finite number");
    if(!R_FINITE(cndtol)) error("'cndtol' is not a valid finite number");

    maxit   = asInteger(getListElement(control, "maxit"));
    if(maxit == NA_INTEGER) error("'maxit' is not an integer");

    outopt[0] = asInteger(getListElement(control, "trace"));
    outopt[1] = asLogical(getListElement(control, "chkjac"));
    outopt[2] = asLogical(rjacobian) ? 1 : 0;

    z = CHAR(STRING_ELT(rglobal, 0));
    if( strcmp(z,"none") == 0 )
        global = 0;
    else if( strcmp(z,"cline") == 0 )
        global = 1;
    else if( strcmp(z,"qline") == 0 )
        global = 2;
    else if( strcmp(z,"gline") == 0 )
        global = 3;
    else if( strcmp(z,"dbldog") == 0 )
        global = 4;
    else if( strcmp(z,"pwldog") == 0 )
        global = 5;
    else if( strcmp(z,"hook") == 0 )
        global = 6;

    z = CHAR(STRING_ELT(rxscalm, 0));
    if( strcmp(z,"fixed") == 0 )
        xscalm = 0;
    else
        xscalm = 1;

    dsub   = asInteger(getListElement(control, "dsub"));
    dsuper = asInteger(getListElement(control, "dsuper"));
    /* -1 means not specified i.e. not active */
    if(dsub == NA_INTEGER || dsub < -1 || dsub > n-2) error("Invalid/impossible value for 'dsub'");
    if(dsuper == NA_INTEGER || dsuper < -1 || dsuper > n-2) error("Invalid/impossible value for 'dsuper'");
    if( (dsub < 0 && dsuper >= 0) || (dsuper < 0 && dsub >= 0) ) error("Both dsub and dsuper must be specified");
    if(method == 1 && dsub == 0 && dsuper == 0) error("method Broyden not implemented for dsub=dsuper=0!");

    /*
     * setup calculation type of jacobian
     */

    jactype = isNull(jac) ? 0 : 1;
    if( dsub >= 0 && dsuper >= 0 ) {
        /*
         * both dsub and dsuper specified and >= 0
         * banded jacobian
         */

        jactype  += 2;
        jacflg[1] = OS->dsub   = dsub;
        jacflg[2] = OS->dsuper = dsuper;
    }
    else {
        jacflg[1] = OS->dsub   = -1;
        jacflg[2] = OS->dsuper = -1;
    }

    jacflg[0] = jactype;

    /* for adjusting step when singular or illconditioned */
    jacflg[3] = asLogical(getListElement(control, "allowSingular")) ? 1 : 0;

    /*
     * query the optimal amount of memory Lapack needs
     * to execute blocked QR code
     * for largish n (500+) this speeds up significantly
     */

    F77_CALL(liqsiz)(&n, &qrwsiz);

    if( qrwsiz <= 0 )
        error("Error in querying amount of workspace for QR routines\n");

    /*
     * allocate memory for Fortran routines
     */

    qrwork  = real_vector(qrwsiz);

    lrwork = 9*n;
    ldr    = n;                             /* leading dimension of rjac */
    lrjac  = method == 1? 2*ldr*n : ldr*n;  /* Broyden needs 2*n columns; Newton needs n columns */

    x       = real_vector(n);
    xp      = real_vector(n);
    fp      = real_vector(n);
    gp      = real_vector(n);
    scalex  = real_vector(n);
    rwork   = real_vector(lrwork);
    rcdwrk  = real_vector(3*n);
    icdwrk  = int_vector(n);
    rjac    = real_vector(lrjac);

    /*
     * setup scaling if specified
     */

    /* copied from code in <Rsource>/src/library/stats/src/optim.c */
    sexp_diag = getListElement(control, "scalex");
    if( LENGTH(sexp_diag) != n )
        error("'scalex' is of the wrong length");
    PROTECT(sexp_diag = coerceVector(sexp_diag,REALSXP)); // 3 protected
    for (i = 0; i < n; i++) {
        scalex[i] = REAL(sexp_diag)[i];
        if( !R_FINITE(scalex[i]) )
            error("'scalex' contains a non-finite value at index=%d\n",i+1);
    }

    /*
     * initialize initial point
     */
    for (i = 0; i < n; i++)
        x[i] = REAL(OS->x)[i];

/*========================================================================*/

    if( outopt[0] == 1)
        trace_header(method, global, xscalm, sigma, delta, stepmx, ftol, xtol, btol, cndtol);

    if( isNull(jac) ) {
        /* numerical jacobian */
        F77_CALL(nwnleq)(x, &n, scalex, &maxit, jacflg, &xtol, &ftol, &btol, &cndtol,
                         &method, &global, &xscalm, &stepmx, &delta, &sigma,
                         rjac, &ldr,
                         rwork, &lrwork, rcdwrk, icdwrk, qrwork, &qrwsiz,
                         FCNJACDUM, &fcnval, outopt, xp, fp, gp, &njcnt, &nfcnt, &iter, &termcd);
    }
    else {
        /* user supplied jacobian */
        if (!isFunction(jac))
            error("jac is not a function!");
        PROTECT(OS->jcall = lang2(jac, OS->x)); // 4 protected
        F77_CALL(nwnleq)(x, &n, scalex, &maxit, jacflg, &xtol, &ftol, &btol, &cndtol,
                         &method, &global, &xscalm, &stepmx, &delta, &sigma,
                         rjac, &ldr,
                         rwork, &lrwork, rcdwrk, icdwrk, qrwork, &qrwsiz,
                         &fcnjac, &fcnval, outopt, xp, fp, gp, &njcnt, &nfcnt, &iter, &termcd);
        UNPROTECT(1); // OS->jcall, 3 protected
    }

/*========================================================================*/

    fcn_message(message, sizeof(message), termcd);

    PROTECT(sexp_x = allocVector(REALSXP,n)); // 4 protected
    for (i = 0; i < n; i++)
        REAL(sexp_x)[i] = xp[i];

    PROTECT(xnames = getAttrib(xstart,R_NamesSymbol)); // 5 protected
    if(!isNull(xnames)) setAttrib(sexp_x, R_NamesSymbol, xnames);

    PROTECT(sexp_fvec = allocVector(REALSXP,n)); // 6 protected
    for (i = 0; i < n; i++)
        REAL(sexp_fvec)[i] = fp[i];

    PROTECT(sexp_info = allocVector(INTSXP,1)); // 7 protected
    INTEGER(sexp_info)[0] = termcd;

    PROTECT(sexp_nfcnt = allocVector(INTSXP,1)); // 8 protected
    INTEGER(sexp_nfcnt)[0] = nfcnt;

    PROTECT(sexp_njcnt = allocVector(INTSXP,1)); // 9 protected
    INTEGER(sexp_njcnt)[0] = njcnt;

    PROTECT(sexp_iter = allocVector(INTSXP,1)); // 10 protected
    INTEGER(sexp_iter)[0] = iter;

    PROTECT(sexp_message = allocVector(STRSXP,1)); // 11 protected
    SET_STRING_ELT(sexp_message, 0, mkChar(message));

    for (i = 0; i < n; i++)
        REAL(sexp_diag)[i] = scalex[i];

    if( outopt[2] == 1 )
        PROTECT(out = allocVector(VECSXP,9)); // 12 protected
    else
        PROTECT(out = allocVector(VECSXP,8)); // 12 protected

    SET_VECTOR_ELT(out, 0, sexp_x);
    SET_VECTOR_ELT(out, 1, sexp_fvec);
    SET_VECTOR_ELT(out, 2, sexp_info);
    SET_VECTOR_ELT(out, 3, sexp_message);
    SET_VECTOR_ELT(out, 4, sexp_diag);
    SET_VECTOR_ELT(out, 5, sexp_nfcnt);
    SET_VECTOR_ELT(out, 6, sexp_njcnt);
    SET_VECTOR_ELT(out, 7, sexp_iter);
    if( outopt[2] == 1 ) {
        PROTECT(sexp_jac = allocMatrix(REALSXP, n, n)); // 13 protected
        SET_VECTOR_ELT(out, 8, sexp_jac);
        pjac = REAL(sexp_jac);
        for(j=0; j < n; j++)
            for(i=0; i < n; i++)
                pjac[j*n+i] = rjac[j*n+i];

        if(!isNull(xnames)) {
            SEXP rcnames;
            PROTECT(rcnames = allocVector(VECSXP, 2)); // 14 protected
            SET_VECTOR_ELT(rcnames, 0, duplicate(xnames));
            SET_VECTOR_ELT(rcnames, 1, duplicate(xnames));
            setAttrib(sexp_jac, R_DimNamesSymbol, rcnames);
            UNPROTECT(1); // rcnames, 13 protected
        }
        UNPROTECT(1); // sexp_jac (and pjac), 12 protected
    }

    if( outopt[2] == 1 )
        PROTECT(out_names = allocVector(STRSXP,9)); // 13 protected
    else
        PROTECT(out_names = allocVector(STRSXP,8)); // 13 protected

    SET_STRING_ELT(out_names, 0, mkChar("x"));
    SET_STRING_ELT(out_names, 1, mkChar("fvec"));
    SET_STRING_ELT(out_names, 2, mkChar("termcd"));
    SET_STRING_ELT(out_names, 3, mkChar("message"));
    SET_STRING_ELT(out_names, 4, mkChar("scalex"));
    SET_STRING_ELT(out_names, 5, mkChar("nfcnt"));
    SET_STRING_ELT(out_names, 6, mkChar("njcnt"));
    SET_STRING_ELT(out_names, 7, mkChar("iter"));
    if( outopt[2] == 1 )
        SET_STRING_ELT(out_names, 8, mkChar("jac"));

    setAttrib(out, R_NamesSymbol, out_names);

    UNPROTECT(13);

    return out;
}
