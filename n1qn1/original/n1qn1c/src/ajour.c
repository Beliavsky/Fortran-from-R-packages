/* ajour.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#define DOUBLE_EPS     DBL_EPSILON
//#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
#define safe_zero(a) (a)


#include <stdio.h>
#include <float.h>

#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

#include "n1qn1.h"

/* Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) INRIA */

/* Copyright (C) 2012 - 2016 - Scilab Enterprises */

/* This file is hereby licensed under the terms of the GNU GPL v2.0, */
/* pursuant to article 5.3.4 of the CeCILL v.2.1. */
/* This file was originally licensed under the terms of the CeCILL v2.1, */
/* and continues to be available under such terms. */
/* For more information, see the COPYING file which you should have received */
/* along with this program. */

/* Subroutine */ int ajour_(int *mode, int *n, int *nc, int *
	nr, double *h__, double *w, int *indi)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Local variables */
    double a, b, c__;
    int i__, j, k;
    double u, v, h1;
    int i1;
    double h2, ai, di;
    int ii, ij, ik, nh, nj, nk, nl, ko;
    double wi;
    int nw;
    double di1;
    int nh1, nr1, nr2, inc;
    double hij;
    int nii, nkk, nrr, inc1, nsaut;



/*     mode = +1 factorise la ligne nc (indices de depart) */
/*          = -1 defactorise   ' */
/*     nr nbre de lignes factorisees */
/*     h mat de dim n */
/*     w,d vect de travail */
/*     indi(i) ligne ou est stockee la ligne i de depart */

    /* Parameter adjustments */
    --indi;
    --w;
    --h__;

    /* Function Body */
    inc = indi[*nc];
    nr1 = *nr + 1;
    nr2 = *nr - 1;
    nrr = *n - *nr;
    nii = *n - inc;
    nkk = *nr - inc;
    if (*mode == -1) {
	goto L240;
    }

/*      addition d'une ligne a l */

/*          stockage des elements de la colonne inc dans w */
    nsaut = nii + 1;
    nh = inc * (*n + 1) - inc * (inc + 1) / 2;
    nw = *n;
    if (inc == *n) {
	goto L20;
    }
    i__1 = nii;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[nw] = h__[nh];
	--nw;
	--nh;
    }
L20:
    w[nr1] = h__[nh];
    --nh;
    if (inc == nr1) {
	goto L60;
    }
    i__1 = inc - nr1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nl = nii + i__ - 1;
	if (nl == 0) {
	    goto L35;
	}
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh + nsaut] = h__[nh];
	    --nh;
	}
L35:
	w[nw] = h__[nh];
	--nw;
	--nh;
	++nsaut;
    }
    i__1 = inc - nr1;
    for (j = 1; j <= i__1; ++j) {
	h__[nh + nsaut] = h__[nh];
	--nh;
    }

L60:
    --nw;
    nsaut = 1;
    if (*nr == 0) {
	goto L125;
    }
    if (inc == *n) {
	goto L80;
    }
    i__1 = nii;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh + nsaut] = h__[nh];
	--nh;
    }
L80:
    if (*nr == 1) {
	goto L110;
    }
    i__1 = nr2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[nw] = h__[nh];
	--nw;
	--nh;
	++nsaut;
	if (*n == nr1) {
	    goto L100;
	}
	i__2 = *n - nr1;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh + nsaut] = h__[nh];
	    --nh;
	}
    }
L100:
L110:
    w[nw] = h__[nh];
    --nh;
    ++nsaut;
    if (inc == nr1) {
	goto L125;
    }
    i__1 = inc - nr1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh + nsaut] = h__[nh];
	--nh;
    }
/*         mise a jour de l */
L125:
    if (*nr != 0) {
	goto L130;
    }
    if (w[1] > 0.) {
	goto L220;
    }
    *mode = -1;
    return 0;
L130:
    if (*nr == 1) {
	goto L160;
    }
    i__1 = *nr;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ij = i__;
	i1 = i__ - 1;
	v = w[i__];
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    v -= h__[ij] * w[j];
	    ij = ij + *nr - j;
	}
	w[i__] = v;
    }
L160:
    ij = 1;
    v = w[nr1];
    i__1 = *nr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__];
	hij = safe_zero(h__[ij]);
/* Computing 2nd power */
	d__1 = wi;
	v -= d__1 * d__1 / hij;
	w[i__] = wi / hij;
	ij = ij + nr1 - i__;
    }
    if (v > 0.) {
	goto L180;
    }
    *mode = -1;
    return 0;
L180:
    w[nr1] = v;
/*          stockage de w dans h */
    nh = *nr * (*nr + 1) / 2;
    nw = nr1;
    nsaut = nw;
    h__[nh + nsaut] = w[nw];
    --nw;
    --nsaut;
    if (*nr == 1) {
	goto L220;
    }
    i__1 = nr2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh + nsaut] = w[nw];
	--nw;
	--nsaut;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh + nsaut] = h__[nh];
	    --nh;
	}
    }
L220:
    h__[nr1] = w[1];
    if (*n == nr1) {
	goto L233;
    }
    nh1 = *nr * (*n + 1) - *nr * (*nr + 1) / 2 + 1;
    nw = nr1;
    i__1 = *n - nr1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh1 + i__] = w[nw + i__];
    }
/*          mise a jour de indi */
L233:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = indi[i__];
	if (ii <= *nr || ii >= inc) {
	    goto L235;
	}
	indi[i__] = ii + 1;
    }
L235:
    ++(*nr);
    indi[*nc] = *nr;
    *mode = 0;
    return 0;

/*      soustraction d'une ligne a l */

/*          recherche des composantes de h */
L240:
    i__1 = *nr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ik = i__;
	ij = inc;
	ii = 1;
	ko = min(ik,inc);
	v = 0.;
	if (ko == 1) {
	    goto L252;
	}
	i__2 = ko - 1;
	for (k = 1; k <= i__2; ++k) {
	    nk = nr1 - k;
	    v += h__[ij] * h__[ik] * h__[ii];
	    ij = ij + nk - 1;
	    ii += nk;
	    ik = ik + nk - 1;
	}
L252:
	a = 1.;
	b = 1.;
	if (ko == i__) {
	    goto L253;
	}
	a = h__[ik];
L253:
	if (ko == inc) {
	    goto L260;
	}
	b = h__[ij];
	w[i__] = v + a * b * h__[ii];
    }
L260:
/*          mise a jour de l */
    if (inc == *nr) {
	goto L315;
    }
    inc1 = inc - 1;
    nh = inc1 * nr1 - inc1 * inc / 2 + 2;
    nh1 = nh + nkk;
    di = h__[nh - 1];
    i__1 = nkk;
    for (j = 1; j <= i__1; ++j) {
	di1 = h__[nh1];
	++nh1;
	a = h__[nh];
	ai = a * di;
/* Computing 2nd power */
	d__1 = a;
	c__ = safe_zero(d__1 * d__1 * di + di1);
	h__[nh] = c__;
	++nh;
	if (j == nkk) {
	    goto L315;
	}
	i__2 = nkk - j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    h1 = h__[nh];
	    h2 = h__[nh1];
	    u = ai * h1 + h2 * di1;
	    h__[nh] = u / c__;
	    h__[nh1] = -h1 + a * h2;
	    ++nh;
	    ++nh1;
	}
	++nh;
	di = di * di1 / c__;
    }
/*          condensation de la matrice l */
L315:
    nh = inc + 1;
    nsaut = 1;
    nj = *nr - 2;
    if (inc == 1) {
	++nj;
    }
    if (*nr == 1) {
	goto L440;
    }
    i__1 = nr2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nj;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh - nsaut] = h__[nh];
	    ++nh;
	}
	++nsaut;
	++nh;
	if (i__ == inc - 1) {
	    goto L430;
	}
	--nj;
	if (nj == 0) {
	    goto L440;
	}
    }
L430:
/*          mise a jour de la matrice h */
L440:
    nh = *nr * nr2 / 2 + 1;
    nw = 1;
    nsaut = *nr;
    if (inc == 1) {
	goto L470;
    }
    i__1 = inc - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh] = w[nw];
	++nw;
	--nsaut;
	if (*n == *nr) {
	    goto L455;
	}
	i__2 = nrr;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh + j] = h__[nh + nsaut + j];
	}
L455:
	nh = nh + nrr + 1;
    }
L470:
    ++nw;
    if (*nr == *n) {
	goto L485;
    }
    i__1 = nrr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[*nr + i__] = h__[nh + nsaut + i__ - 1];
    }
    nsaut += nrr;
L485:
    if (inc == *nr) {
	goto L510;
    }
    i__1 = nkk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--nsaut;
	h__[nh] = w[nw];
	++nw;
	if (*nr == *n) {
	    goto L495;
	}
	i__2 = nrr;
	for (j = 1; j <= i__2; ++j) {
	    h__[nh + j] = h__[nh + nsaut + j];
	}
L495:
	nh = nh + nrr + 1;
    }
L510:
    h__[nh] = w[inc];
    if (*nr == *n) {
	goto L540;
    }
    i__1 = nrr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[nh + i__] = w[*nr + i__];
    }
/*          mise a jour de indi */
L540:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = indi[i__];
	if (ii <= inc || ii > *nr) {
	    goto L550;
	}
	indi[i__] = ii - 1;
    }
L550:
    indi[*nc] = *nr;
    --(*nr);
    *mode = 0;
    return 0;
} /* ajour_ */

