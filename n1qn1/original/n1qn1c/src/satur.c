/* satur.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int satur_(int *n, double *x, double *binf, 
	double *bsup, double *d__, double *ttmin, double *
	ttsup, double *topt, double *tg, double *td, double *
	tmi, int *icoi, int *icos, int *iproj)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Local variables */
    double e;
    int i__;
    double ep, tb;
    int inf;
    double cres;


/*      subroutine calculant ,dans un intervalle donne, un pas proche */
/*      de tmi saturant une contrainte */
/*         topt:pas calculer (=0 s'il n'existe pas un tel pas         (s) */
/*        ttmin,ttsup:bornes de l'intervalle dans lequel doit */
/*         etre topt                                                  (e) */
/*        tmi:pas au voisinnage duquel on calcul topt                 (e) */
/*        iproj:indicateur de projection                              (e) */
/*             =0:on cherche un pas saturant une contrainte dans */
/*                 l'intervalle ttmin,ttsup */
/*             =1:on cherche un pas dans l'intervalle tg,td et on */
/*                le ramene dans l'intervalle ttmin,ttsup */
/*                par projection */
/*       icos:indice de la variable saturee par la borne superieure */
/*       icoi:indice de la variable saturee par la borne inferieure */
/*       inf:indicateur pour l initialisation de icoi et icos */
/*            =0:la borne superieure est atteinte */
/*            =1:la borne superieure est atteinte */
/*            =2:le pas est obtenu par projection sur ttmin ttsup */


    /* Parameter adjustments */
    --d__;
    --bsup;
    --binf;
    --x;

    /* Function Body */
    *icoi = 0;
    *icos = 0;
    ep = *tmi;

/*        boucle */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	inf = 0;
/*        calcul du pas saturant la ieme contrainte:tb */
	cres = d__[i__];
	if (cres < 0.) {
	    goto L61;
	} else if (cres == 0.) {
	    goto L70;
	} else {
	    goto L62;
	}
L61:
	tb = (binf[i__] - x[i__]) / d__[i__];
	inf = 1;
	goto L63;
L62:
	tb = (bsup[i__] - x[i__]) / d__[i__];
L63:
	if (tb > *ttsup || tb < *ttmin) {
/*        projection de tb sur l intervalle ttmin,ttsup */
	    if (*iproj == 0 || tb < *tg || tb > *td) {
		goto L70;
	    }
	    tb = max(tb,*ttmin);
	    tb = min(tb,*ttsup);
	    inf = 2;
	}
/*        recherche du pas le plus proche de tmi */
	e = (d__1 = tb - *tmi, fabs(d__1));
	if (e >= ep) {
	    goto L70;
	}
	*topt = tb;
	ep = e;
/*        mise a jour de icoi,icos */
	*icoi = 0;
	*icos = 0;
	if (inf == 0) {
	    *icos = i__;
	}
	if (inf == 1) {
	    *icoi = i__;
	}
L70:
	;
    }
    return 0;
} /* satur_ */

