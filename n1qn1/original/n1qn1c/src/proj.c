/* proj.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
#include <float.h>
#include <math.h>

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

/* Subroutine */ int proj_(int *n, double *binf, double *bsup, 
	double *x)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int i__;

    /* Parameter adjustments */
    --x;
    --bsup;
    --binf;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
/* Computing MIN */
	d__3 = x[i__], d__4 = bsup[i__];
	d__1 = binf[i__], d__2 = min(d__3,d__4);
	x[i__] = max(d__1,d__2);
    }
    return 0;
} /* proj_ */

