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
 *  AUTHORS: Christophe Dutang, based on code eigen,
 *
 *  i.e., function 'modLa_rg' and 'modLa_dgesv' in R's
 *  <Rsource>/src/modules/lapack/lapack.c, used in eigen()
 */

#include "expm-eigen.h"


void expm_eigen(double *x, int n, double *z, double tol)
{
    if (n == 1)
        z[0] = exp(x[0]);		/* scalar exponential */
    else
    {
        const int nsqr = n * n;

        const Rcomplex cone = { .r = 1., .i = 0.}, czero = { .r = 0., .i = 0.};
        int i, j;
        int info, lwork, is_conjug, is_diag;
        double onenorm, rcond, tmp;

        /* Arrays */
        int *ipiv = (int *) R_alloc(n, sizeof(int)); /* permutation vector */
        double *left, *right, *workdiag; /* left and right eigenvectors and workspace for diagonalisation */
        double *wR = (double *) R_alloc(n, sizeof(double)); /* real part of eigenvalues */
        double *wI = (double *) R_alloc(n, sizeof(double)); /* imaginary part of eigenvalues */
        double *rworksing = (double *) R_alloc(2*n, sizeof(double)); /* working vector to test the singularity */
        Rcomplex *eigvect = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* (right) eigenvectors matrix */
        Rcomplex *eigvectinv = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* its inverse */
        Rcomplex *expeigval; /* complex matrix diag(exp(eigenvalues)) */
        Rcomplex *ctmp = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* temp working variable */
        Rcomplex *worksing = (Rcomplex *) R_alloc(2*n, sizeof(Rcomplex)); /* workspace to test the singularity */

        Memcpy(z, x, nsqr);

        /* Test if x is diagonalisable by computing its eigenvalues and (right) eigenvectors */
        /* code based on modLa_rg in lapack.c, used in eigen.R */
        left = (double *) 0;
        right = (double *) R_alloc(nsqr, sizeof(double));

        /* 1 - ask for optimal size of work array */
        lwork = -1;
        F77_CALL(dgeev)("N", "V", &n, z, &n, wR, wI,
			left, &n, right, &n, &tmp, &lwork, &info FCONE FCONE);
        if (info != 0)
            error(_("error code %d from Lapack routine dgeev"), info);
        lwork = (int) tmp;
        workdiag = (double *) R_alloc(lwork, sizeof(double));

        /* 2 - compute eigenvalues and (right) eigenvectors */
        F77_CALL(dgeev)("N", "V", &n, z, &n, wR, wI,
			left, &n, right, &n, workdiag, &lwork, &info FCONE FCONE);
        if (info != 0)
            error(_("error code %d from Lapack routine dgeev"), info);

        /* try to invert the eigenvectors matrix */
        /* 1 - build the Rcomplex matrix with eigenvectors */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                is_conjug = 0;
                if(i < n-1)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i+1] && wI[i] == -wI[i+1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[i * n + j];
                        eigvect[i * n + j].i = right[(i+1) * n + j];
                    }
                }
                if(i > 0)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i-1] && wI[i] == -wI[i-1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[(i-1) * n + j];
                        eigvect[i * n + j].i = -right[i * n + j];
                    }
                }
                /* real eigenvalues */
                if(!is_conjug)
                {
                    eigvect[i * n + j].r = right[i * n + j];
                    eigvect[i * n + j].i = 0.0;
                }
                /* eigvectinv initialise with the identity matrix */
                eigvectinv[i * n +j].r = (i == j) ? 1.0 : 0.0;
                eigvectinv[i * n +j].i = 0.0;
            }
        }

        /* 2 - store the matrix eigvect (because function zgesv will change it) */
        Memcpy(ctmp, eigvect, nsqr);

        /* 3 - solve a linear complex equation system with eigvectinv equals
         * to matrix identity. hence, on exit eigvectinv contains the
         * inverse of complex matrix eigvect. code base on solve.R */
        F77_CALL(zgesv)(&n, &n, eigvect, &n, ipiv, eigvectinv, &n, &info);

        if (info > 0)
            is_diag = 0; //matrix eigvect is exactly singular.
        if (info < 0)
            error(_("argument %d of Lapack routine dgesv had invalid value"), -info);
        if (info == 0)
            is_diag = 1;

        /* check if matrix eigvectinv is numerically singular */
        if (is_diag)
        {
            /* compute the reciprocal condition number of eigvectinv. */

            /* 1 - compute the one norm of the matrix eigvectinv */
            onenorm = F77_CALL(zlange)("1", &n, &n, eigvectinv, &n, (double*) NULL FCONE);

            /* 2 - estimates the reciprocal of the condition number
             * when the one norm is used. */
            F77_CALL(zgecon)("1", &n, eigvectinv, &n, &onenorm, &rcond,
			     worksing, rworksing, &info FCONE);

            if (rcond < tol)
                is_diag=0;
        }


        if (is_diag)
        {

            /* x is diagonalisable so
             * compute complex matrix operations :
             * eigvect %*% diag(exp(eigenvalues)) %*% eigvectinv */

            /* 1 - expeigval is the complex matrix diag(exp(eigenvalues)) */
            /* code based on z_exp in complex.c */
            expeigval = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex));
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if(i == j)
                    {
                        expeigval[i * n +j].r = exp(wR[i]) * cos(wI[i]);
                        expeigval[i * n +j].i = exp(wR[i]) * sin(wI[i]);
                    }
                    else
                    {
                        expeigval[i * n +j].r = 0.0;
                        expeigval[i * n +j].i = 0.0;
                    }
                }
            }

            /* 2 - restore the matrix eigvect */
            Memcpy(eigvect, ctmp, nsqr);

            /* 3 - compute (complex) matrix product: ctmp <- eigvect * expeigval */
            F77_CALL(zgemm)("N", "N", &n, &n, &n, &cone, eigvect, &n,
			    expeigval, &n, &czero, ctmp, &n FCONE FCONE);

            /* 4 - compute (complex) matrix product: expeigval <- ctmp * eigvectinv */
            F77_CALL(zgemm)("N", "N", &n, &n, &n, &cone, ctmp, &n,
			    eigvectinv, &n, &czero, expeigval, &n FCONE FCONE);

            /* store the real part in z */
            /* the matrix exponential is always real,
             * even if x has complex eigen values. */
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    z[i * n + j] = expeigval[i * n + j].r;

        }
	else
	    expm(x, n, z, Ward_2);
    }
}

/* Main function, the only one used by .Call(). */
SEXP do_expm_eigen(SEXP x, SEXP tolin)
{
    SEXP dims, z;
    int n, nprot = 0;
    double *rx, *rz;
    double tol = asReal(tolin);

    if (!isNumeric(x) || !isMatrix(x))
	error(_("invalid argument: not a numeric matrix"));
    if (isInteger(x)) {
	x = PROTECT(coerceVector(x, REALSXP)); nprot++;
    }
    rx = REAL(x);

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

    expm_eigen(rx, n, rz, tol);

    setAttrib(z, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(nprot);
    return z;
}
