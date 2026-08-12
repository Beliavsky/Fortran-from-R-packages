/*  ===== File part of R package expm =====
 *
 *  Function to compute the matrix exponential
 *
 *     exp(M) = sum(n = 0:Inf; M^n / n!),
 *
 *  where M is an (n x n) matrix.
 *
 *  The functions therein use LAPACK and BLAS routines. Nicely
 *  formatted man pages for these can be found at
 *
 *    <http://www.mathkeisan.com/UsersGuide/E/>
 *
 *  AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Christophe
 *  Dutang, based on code in package Matrix.
 */

#include "expm.h"


/* For matrix exponential calculations. Pade constants
 *
 *   n_{pqj} = [(p + q - j)! p!]/[(p + q)! j! (p - j)!]
 *
 * and
 *
 *   d_{pqj} = [(p + q - j)! q!]/[(p + q)! j! (q - j)!]
 *
 * for p = q = 8 and j = 1, ..., 8.
 */
const static double padec88 [] =
{
    5.0000000000000000e-1,
    1.1666666666666667e-1,
    1.6666666666666667e-2,
    1.6025641025641026e-3,
    1.0683760683760684e-4,
    4.8562548562548563e-6,
    1.3875013875013875e-7,
    1.9270852604185938e-9,
};


/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. Mostly lifted from the core of fonction
 * expm() of package Matrix, which is itself based on the function of
 * the same name in Octave.
 */
void expm(double *x, int n, double *z, precond_type precond_kind)
{
    if (n == 1)
	z[0] = exp(x[0]);		/* scalar exponential */
    else
    {
	/* Constants */
	const double one = 1.0, zero = 0.0;
	const int i1 = 1, nsqr = n * n, np1 = n + 1;
	const Rboolean do_scale = (precond_kind == Ward_2 || precond_kind == Ward_buggy_octave);
	/* Variables */
	int i, j, is_uppertri = TRUE;;
	int ilo, ihi, iloscal, ihiscal, info, sqrpowscal;
	double infnorm, trshift, m1pj = -1;

	/* Arrays */
	int *pivot    = (int *) R_alloc(n, sizeof(int)); /* pivot vector */
	double *scale =
	    do_scale ?  (double *) R_alloc(n, sizeof(double)) : NULL; /* scale array */
	double *perm  = (double *) R_alloc(n, sizeof(double));/* permutation/sc array */

	double *work  = (double *) R_alloc(nsqr, sizeof(double)); /* workspace array */
	double *npp   = (double *) R_alloc(nsqr, sizeof(double)); /* num. power Pade */
	double *dpp   = (double *) R_alloc(nsqr, sizeof(double)); /* denom. power Pade */

	Memcpy(z, x, nsqr);

	/* Check if matrix x is upper triangular; stop checking as
	 * soon as a non-zero value is found below the diagonal. */
	for (i = 0; i < n - 1 && is_uppertri; i++)
	    for (j = i + 1; j < n; j++)
		if (!(is_uppertri = x[i * n + j] == 0.0))
		    break;


	/* Step 1 of preconditioning: shift diagonal by average diagonal. */
	trshift = 0.0;
	for (i = 0; i < n; i++)
	    trshift += x[i * np1];
	trshift /= n;		/* average diagonal element */
	if (trshift > 0.0)
	    for (i = 0; i < n; i++)
		z[i * np1] -= trshift;

	/* Step 2 of preconditioning: balancing with dgebal. */
	if(do_scale) {
	    if (is_uppertri) {
		/* no need to permute if x is upper triangular */
		ilo = 1;
		ihi = n;
	    }
	    else {
		F77_CALL(dgebal)("P", &n, z, &n, &ilo, &ihi, perm, &info FCONE);
		if (info)
		    error(_("LAPACK routine dgebal returned info code %d when permuting"), info);
	    }
	    F77_CALL(dgebal)("S", &n, z, &n, &iloscal, &ihiscal, scale, &info FCONE);
	    if (info)
		error(_("LAPACK routine dgebal returned info code %d when scaling"), info);
	}
	else if(precond_kind == Ward_1) {

	    F77_CALL(dgebal)("B", &n, z, &n, &ilo, &ihi, perm, &info FCONE);
	    if (info)
		error(_("LAPACK' dgebal(\"B\",.) returned info code %d"), info);

	}
	else {
	    error(_("invalid 'precond_kind: %d"), precond_kind);
	}


	/* Step 3 of preconditioning: Scaling according to infinity
	 * norm (a priori always needed). */
	infnorm = F77_CALL(dlange)("I", &n, &n, z, &n, work FCONE);
	sqrpowscal = (infnorm > 0) ? imax2((int) 1 + log(infnorm)/M_LN2, 0) : 0;
	if (sqrpowscal > 0) {
	    double scalefactor = R_pow_di(2, sqrpowscal);
	    for (i = 0; i < nsqr; i++)
		z[i] /= scalefactor;
	}

	/* Pade approximation (p = q = 8): compute x^8, x^7, x^6,
	 * ..., x^1 */
	for (i = 0; i < nsqr; i++)
	{
	    npp[i] = 0.0;
	    dpp[i] = 0.0;
	}
	for (j = 7; j >= 0; j--)
	{
	    /* npp = z * npp + padec88[j] * z */
	    F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, npp,
			     &n, &zero, work, &n FCONE FCONE);

	    /* npp <- work + padec88[j] * z */
	    for (i = 0; i < nsqr; i++)
		npp[i] = work[i] + padec88[j] * z[i];
	    /* dpp = z * dpp + (-1)^j * padec88[j] * z */
	    F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, dpp,
			     &n, &zero, work, &n FCONE FCONE);
	    for (i = 0; i < nsqr; i++)
		dpp[i] = work[i] + m1pj * padec88[j] * z[i];
	    m1pj *= -1; /* (-1)^j */
	}
	/* power 0 */
	for (i = 0; i < nsqr; i++)
	    dpp[i] *= -1.0;
	for (j = 0; j < n; j++)
	{
	    npp[j * np1] += 1.0;
	    dpp[j * np1] += 1.0;
	}

	/* Pade approximation is (dpp)^-1 * npp. */
	F77_CALL(dgetrf) (&n, &n, dpp, &n, pivot, &info);
	if (info)
	    error(_("LAPACK routine dgetrf returned info code %d"), info);
	F77_CALL(dgetrs) ("N", &n, &n, dpp, &n, pivot, npp, &n, &info FCONE);
	if (info)
	    error(_("LAPACK routine dgetrs returned info code %d"), info);

	Memcpy(z, npp, nsqr);

	/* Now undo all of the preconditioning */
	/* Preconditioning 3: square the result for every power of 2 */
	while (sqrpowscal--)
	{
	    F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, z, &n,
			    z, &n, &zero, work, &n FCONE FCONE);
	    Memcpy(z, work, nsqr);
	}


	/* Preconditioning 2: Inversion of 'dgebal()' :
	 * ------------------ Note that dgebak() seems *not* applicable */

	/* Step 2 a)  apply inverse scaling */
	if(do_scale) {
	    for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		    z[i + j * n] *= scale[i]/scale[j];
	}
	else if(precond_kind == Ward_1) { /* here, perm[ilo:ihi] contains scale[] */
	    for (j = 0; j < n; j++) {
		double sj = ((ilo-1 <= j && j < ihi)? perm[j] : 1.);
		for (i = 0; i < ilo-1; i++)	z[i + j * n] /= sj;
		for (i = ilo-1; i < ihi; i++)	z[i + j * n] *= perm[i]/sj;
		for (i = ihi+1; i < n; i++)	z[i + j * n] /= sj;
	    }
	}

	/* 2 b) Inverse permutation  (if not the identity permutation) */

	if (ilo != 1 || ihi != n) {

	    if(precond_kind == Ward_buggy_octave) {

		/* inverse permutation vector */
		int *invP  = (int *) R_alloc(n, sizeof(int));

		/* balancing permutation vector */
		for (i = 0; i < n; i++)
		    invP[i] = i;	/* identity permutation */

		/* leading permutations applied in forward order */
		for (i = 0; i < (ilo - 1); i++)
		{
		    int p_i = (int) (perm[i]) - 1;
		    int tmp = invP[i]; invP[i] = invP[p_i]; invP[p_i] = tmp;
		}

		/* trailing permutations applied in reverse order */
		for (i = n - 1; i >= ihi; i--)
		{
		    int p_i = (int) (perm[i]) - 1;
		    int tmp = invP[i]; invP[i] = invP[p_i]; invP[p_i] = tmp;
		}

		/* construct inverse balancing permutation vector */
		Memcpy(pivot, invP, n);
		for (i = 0; i < n; i++)
		    invP[pivot[i]] = i;

		/* apply inverse permutation */
		Memcpy(work, z, nsqr);
		for (j = 0; j < n; j++)
		    for (i = 0; i < n; i++)
			z[i + j * n] = work[invP[i] + invP[j] * n];

	    }
	    else if(precond_kind == Ward_2 || precond_kind == Ward_1) {

		/* ---- new code by Martin Maechler ----- */

#define SWAP_ROW(I,J) F77_CALL(dswap)(&n, &z[(I)], &n, &z[(J)], &n)

#define SWAP_COL(I,J) F77_CALL(dswap)(&n, &z[(I)*n], &i1, &z[(J)*n], &i1)

#define RE_PERMUTE(I)				\
		int p_I = (int) (perm[I]) - 1;	\
		SWAP_COL(I, p_I);		\
		SWAP_ROW(I, p_I)

		/* reversion of "leading permutations" : in reverse order */
		for (i = (ilo - 1) - 1; i >= 0; i--) {
		    RE_PERMUTE(i);
		}

		/* reversion of "trailing permutations" : applied in forward order */
		for (i = (ihi + 1) - 1; i < n; i++) {
		    RE_PERMUTE(i);
		}

	    }
/* 	    else if(precond_kind == Ward_1) { */

/* 	    } */

	}
	/* Preconditioning 1: Trace normalization */
	if (trshift > 0)
	{
	    double mult = exp(trshift);
	    for (i = 0; i < nsqr; i++)
		z[i] *= mult;
	}

    }
}

/* Main function, the only one used by .External(). */
SEXP do_expm(SEXP x, SEXP kind)
{
    SEXP dims, z;
    int n, nprot = 0;
    double *rx, *rz;
    const char *ch_kind = CHAR(asChar(kind));
    precond_type PC_kind = Ward_2 /* -Wall */;

    if (!isNumeric(x) || !isMatrix(x))
	error(_("invalid argument: not a numeric matrix"));
    if (isInteger(x)) {
	x = PROTECT(coerceVector(x, REALSXP)); 	nprot++;
    }
    rx = REAL(x);

    if(!strcmp(ch_kind, "Ward77")) {
	PC_kind = Ward_2;
    } else if(!strcmp(ch_kind, "buggy_Ward77")) {
	PC_kind = Ward_buggy_octave;
    } else if(!strcmp(ch_kind, "Ward77_1")) {
	PC_kind = Ward_1;
    } else
	error(_("invalid 'kind' argument: %s\n"), ch_kind);

    dims = getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    if (n != INTEGER(dims)[1])
	error(_("non-square matrix"));
    if (n == 0) {
	UNPROTECT(nprot);
	return(allocMatrix(REALSXP, 0, 0));
    }

    PROTECT(z = allocMatrix(REALSXP, n, n)); nprot++;
    rz = REAL(z);

    expm(rx, n, rz, PC_kind);
/*  ---- */
    setAttrib(z, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(nprot);
    return z;
}
