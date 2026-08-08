/* calmaj.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
/* Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) INRIA */

/* Copyright (C) 2012 - 2016 - Scilab Enterprises */

/* This file is hereby licensed under the terms of the GNU GPL v2.0, */
/* pursuant to article 5.3.4 of the CeCILL v.2.1. */
/* This file was originally licensed under the terms of the CeCILL v2.1, */
/* and continues to be available under such terms. */
/* For more information, see the COPYING file which you should have received */
/* along with this program. */

/* Subroutine */ int calmaj_(double *dh, int *n, double *g1, 
	double *sig, double *w, int *ir, int *mk, double *
	epsmc, int *nfac)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int i__, j, k, nfac1, n2fac, nnfac;
    extern /* Subroutine */ int majour_(double *, double *, 
	    double *, int *, double *, int *, int *, 
	    double *);


/*     subroutine de qnbd */
    /* Parameter adjustments */
    --w;
    --g1;
    --dh;

    /* Function Body */
    if (*nfac == *n) {
	goto L50;
    }
    nfac1 = *nfac + 1;
    nnfac = *n - *nfac;
    n2fac = *nfac * nfac1 / 2;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = g1[i__] * *sig;
    }
    k = n2fac;
    if (*nfac == 0) {
	goto L25;
    }
    i__1 = *nfac;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = nfac1; i__ <= i__2; ++i__) {
	    ++k;
	    dh[k] += g1[i__] * w[j];
	}
    }
L25:
    k = n2fac + *nfac * nnfac;
    i__1 = *n;
    for (j = nfac1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    ++k;
	    dh[k] += g1[i__] * w[j];
	}
    }
L50:
    *ir = *nfac;
    if (*nfac == 0) {
	return 0;
    }
    majour_(&dh[1], &g1[1], &w[1], nfac, sig, ir, mk, epsmc);
    return 0;
} /* calmaj_ */

