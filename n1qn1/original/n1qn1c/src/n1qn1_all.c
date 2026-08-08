/* n1qn1_all.f -- translated by f2c (version 20160102).
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
#include <R.h>

#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

#define DOUBLE_EPS     DBL_EPSILON
//#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
#define safe_zero(a) (a)

#include "n1qn1.h"

/* Table of constant values */

static int c__1 = 1;
static double c_b32 = 0.;

/*     Modified by Matthew Fidler in 2017 for different outputs to the R console */
int vff_(int *n, double *g)
{
    /* System generated locals */
    int ret_val, i__1;

    /* Local variables */
    int i__;
    int ret;

    /* Parameter adjustments */
    --g;

    /* Function Body */
    ret = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (!R_FINITE(g[i__])) {
	    ret = 1;
	    goto L7710;
	}
    }
L7710:
    ret_val = ret;
    return ret_val;
} /* vff_ */

/*     Note this should be thread safe since it doesn't use global variables */
/*     Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) 1987 - INRIA - Claude LEMARECHAL */

/* This file must be used under the terms of the CeCILL. */
/* This source file is licensed as described in the file COPYING, which */
/* you should have received as part of this distribution.  The terms */
/* are also available at */
/* http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt */

/* Subroutine */ int n1qn1_(U_fp simul, int *n, double *x, double 
	*f, double *g, double *var, double *eps, int *mode, 
	int *niter, int *nsim, int *imp, double *zm, int *
        izs, float *rzs, double *dzs, int *id)
{
    int nd, nw, nga, ngb, nxa, nxb;
    extern /* Subroutine */ int n1qn1a_(U_fp, int *, double *, 
	    double *, double *, double *, double *, int *,
	     int *, int *, int *, double *, double *, 
	    double *, double *, double *, double *, 
           double *, int *, float *, double *, int *);


/* !but */
/*     minimisation d une fonction reguliere sans contraintes */
/* !origine */
/*     c. lemarechal, inria, 1987 */
/*     Copyright INRIA */
/* !methode */
/*     direction de descente calculee par une methode de quasi-newton */
/*     recherche lineaire de type wolfe */
/* !liste d appel */
/*     simul    : point d'entree au module de simulation (cf normes modulopt i) */
/*     n1qn1 appelle toujours simul avec indic = 4 ; le module de */
/*     simulation doit se presenter sous la forme subroutine simul */
/*     (n,x, f, g, izs, rzs, dzs) et e^tre declare en external dans le */
/*     programme appelant n1qn1. */
/*     n (e)    : nombre de variables dont depend f. */
/*     x (e-s)   : vecteur de dimension n ; en entree le point initial ; */
/*                 en sortie : le point final calcule par n1qn1. */
/*     f (e-s)   : scalaire ; en entree valeur de f en x (initial), en sortie */
/*                 valeur de f en x (final). */
/*     g (e-s)   : vecteur de dimension n : en entree valeur du gradient en x */
/*                 (initial), en sortie valeur du gradient en x (final). */
/*     var (e)   : vecteur strictement positif de dimension n. amplitude de la */
/*                 modif souhaitee a la premiere iteration sur x(i).une bonne */
/*                 valeur est 10% de la difference (en valeur absolue) avec la */
/*                 coordonee x(i) optimale */
/*     eps (e-s) : en entree scalaire definit la precision du test d'arret. */
/*      le programme considere que la convergence est obtenue lorque il lui */
/*      est impossible de diminuer f en attribuant a au moins une coordonnee */
/*      x(i) une variation superieure a eps*var(i). */
/*      en sortie, eps contient le carre de la norme du gradient en x (final). */
/*     mode (e)     : definit l approximation initiale du hessien */
/*                  =1 n1qn1 l initialise lui-meme */
/*                  =2 le hessien est fourni dans zm sous forme compressee (zm */
/*                     contient les colonnes de la partie inferieure du hessien) */
/*     niter (e-s)  : en entree nombre maximal d'iterations : en sortie nombre */
/*                    d'iterations reellement effectuees. */
/*     nsim (e-s)  : en entree nombre maximal d'appels a simul (c'est a dire */
/*         avec indic = 4). en sortie le nombre de tels appels reellement faits. */
/*      imp (e)   : contro^le les messages d'impression : */
/*                  0 rien n'est imprime */
/*                  = 1 impressions initiales et finales */
/*                  = 2 une impression par iteration (nombre d'iterations, */
/*                      nombre d'appels a simul, valeur courante de f). */
/*                  >=3 informations supplementaires sur les recherches */
/*                      lineaires ; */
/*                      tres utile pour detecter les erreurs dans le gradient. */
/*      lp (e)    : le numero du canal de sortie, i.e. les impressions */
/*                  commandees par imp sont faites par write (lp, format). */
/*     zm     : memoire de travail pour n1qn1 de   dimension n*(n+13)/2. */
/*     izs,rzs,dzs memoires reservees au simulateur (cf doc) */

/* ! */
    /* Parameter adjustments */
    --var;
    --g;
    --x;
    --zm;
    --izs;
    --rzs;
    --dzs;

    /* Function Body */
    nd = *n * (*n + 1) / 2 + 1;
    nw = nd + *n;
    nxa = nw + *n;
    nga = nxa + *n;
    nxb = nga + *n;
    ngb = nxb + *n;
    n1qn1a_((U_fp)simul, n, &x[1], f, &g[1], &var[1], eps, mode, niter, nsim, 
	    imp, &zm[1], &zm[nd], &zm[nw], &zm[nxa], &zm[nga], &zm[nxb], &zm[ngb],
	    &izs[1], &rzs[1], &dzs[1], id);
    return 0;
} /* n1qn1_ */

/* Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) INRIA */

/* This file must be used under the terms of the CeCILL. */
/* This source file is licensed as described in the file COPYING, which */
/* you should have received as part of this distribution.  The terms */
/* are also available at */
/* http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt */

/* Subroutine */ int n1qn1a_(S_fp simul, int *n, double *x, 
	double *f, double *g, double *scale, double *acc, 
	int *mode, int *niter, int *nsim, int *iprint, 
	double *h__, double *d__, double *w, double *xa, 
	double *ga, double *xb, double *gb, int *izs, float *rzs,
        double *dzs, int *id)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    double c__;
    int i__, j, k;
    double v;
    int i1;
    double cc, fa, fb, hh;
    int ii, ij, ik, jk, ni, ip, ir, np;
    double gl1, gl2, dga, dgb, dff;
    int ial;
    extern int vff_(int *, double *);
    int nip, itr;
    double fmin, gmin;
    int nfun, isfv;
    double step;
    int indic, iecri;
    double stmin, stepbd, steplb;
    extern /* Subroutine */ int majour_(double *, double *, 
	    double *, int *, double *, int *, int *, 
	    double *);


/*     A (very) few modifs by Bruno (14 March 2005): I have translated some output */
/*     information in english (but they don't use format instruction */
/*     which is put in the second arg of write). Also for the linear */
/*     search output information I divide by the direction vector norm */
/*     to get the "normalized" directionnal derivative. Note that this is */
/*     just for output (the computing code is normally not modified). */
/* (blas routine) added by Bruno to get */
/* a better information concerning directionnal derivative */

/*              calcul initial de fonction-gradient */

    /* Parameter adjustments */
    --gb;
    --xb;
    --ga;
    --xa;
    --w;
    --d__;
    --scale;
    --g;
    --x;
    --h__;
    --izs;
    --rzs;
    --dzs;

    /* Function Body */
    indic = 4;
    (*simul)(&indic, n, &x[1], f, &g[1], &izs[1], &rzs[1], &dzs[1], id);
/*     next line added by Serge to avoid Inf and Nan's (04/2007) */
    if (!R_FINITE(*f) && vff_(n, &g[1]) != 1) {
	indic = -1;
    }
    if (indic > 0) {
	goto L13;
    }
    *acc = 0.;
    *niter = 1;
    *nsim = 1;
    return 0;
L13:
    nfun = 1;
    iecri = 0;
    itr = 0;
    np = *n + 1;
/*                  initialisation du hessien, en fonction de var */
    if (*mode >= 2) {
	goto L60;
    }
L20:
    c__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = c__, d__3 = (d__1 = g[i__] * scale[i__], fabs(d__1));
	c__ = max(d__2,d__3);
    }
    if (c__ <= 0.) {
	c__ = 1.;
    }
    k = *n * np / 2;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__[i__] = 0.;
    }
    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      double tmp = scale[i__] * scale[i__];
      tmp = safe_zero(tmp);
	h__[k] = c__ * .01 / tmp;
	k = k + np - i__;
    }
    goto L100;
/*               factorisation du hessien */
L60:
    if (*mode >= 3) {
	goto L80;
    }
    k = *n;
    if (*n > 1) {
	goto L300;
    }
    if (h__[1] > 0.) {
	goto L305;
    }
    h__[1] = 0.;
    k = 0;
    goto L305;
L300:
    np = *n + 1;
    ii = 1;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
      hh = safe_zero(h__[ii]);
	ni = ii + np - i__;
	if (hh > 0.) {
	    goto L301;
	}
	h__[ii] = 0.;
	--k;
	ii = ni + 1;
	goto L304;
L301:
	ip = ii + 1;
	ii = ni + 1;
	jk = ii;
	i__2 = ni;
	for (ij = ip; ij <= i__2; ++ij) {
	    v = h__[ij] / hh;
	    i__3 = ni;
	    for (ik = ij; ik <= i__3; ++ik) {
		h__[jk] -= h__[ik] * v;
		++jk;
	    }
	    h__[ij] = v;
	}
    }
L304:
    if (h__[ii] > 0.) {
	goto L305;
    }
    h__[ii] = 0.;
    --k;
L305:

    if (k >= *n) {
	goto L100;
    }
L70:
    goto L20;
/*                verification que la diagonale est positive */
L80:
    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (h__[k] <= 0.) {
	    goto L70;
	}
	k = k + np - i__;
    }
/*                quelques initialisations */
L100:
    dff = 0.;
L110:
    fa = *f;
    isfv = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xa[i__] = x[i__];
	ga[i__] = g[i__];
    }
/*                   iteration */
L130:
    ++itr;
    ial = 0;
    if (itr > *niter) {
	goto L250;
    }
    ++iecri;
    if (iecri != -(*iprint)) {
	goto L140;
    }
    iecri = 0;
    indic = 1;
    (*simul)(&indic, n, &x[1], f, &g[1], &izs[1], &rzs[1], &dzs[1], id);
/*     error in user function */
    if (indic == 0) {
	goto L250;
    }
/*               calcul de la direction de recherche */
L140:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = -ga[i__];
    }
    w[1] = d__[1];
    if (*n > 1) {
	goto L400;
    }
    d__[1] /= h__[1];
    goto L412;
L400:
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ij = i__;
	i1 = i__ - 1;
	v = d__[i__];
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    v -= h__[ij] * d__[j];
	    ij = ij + *n - j;
	}
	w[i__] = v;
	d__[i__] = v;
    }
    d__[*n] /= h__[ij];
    np = *n + 1;
    i__1 = *n;
    for (nip = 2; nip <= i__1; ++nip) {
	i__ = np - nip;
	ii = ij - nip;
	v = d__[i__] / safe_zero(h__[ii]);
	ip = i__ + 1;
	ij = ii;
	i__2 = *n;
	for (j = ip; j <= i__2; ++j) {
	    ++ii;
	    v -= h__[ii] * d__[j];
	}
	d__[i__] = v;
    }
L412:
/*               calcul du pas minimum */
/*               et de la derivee directionnelle initiale */
    c__ = 0.;
    dga = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = c__, d__3 = (d__1 = d__[i__] / safe_zero(scale[i__]), fabs(d__1));
	c__ = max(d__2,d__3);
	dga += ga[i__] * d__[i__];
    }
/*               test si la direction est de descente */
    if (dga >= 0.) {
	goto L240;
    }
/*               initialisation du pas */
    stmin = 0.;
    stepbd = 0.;
    steplb = *acc / safe_zero(c__);
    fmin = fa;
    gmin = dga;
    step = 1.;
    if (dff <= 0.) {
/* Computing MIN */
	d__1 = step, d__2 = 1. / safe_zero(c__);
	step = min(d__1,d__2);
    }
    if (dff > 0.) {
/* Computing MIN */
	d__1 = step, d__2 = (dff + dff) / safe_zero(-dga);
	step = min(d__1,d__2);
    }
L170:
    c__ = stmin + step;
    if (nfun >= *nsim) {
	goto L250;
    }
    ++nfun;
/*              calcul de fonction-gradient */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xb[i__] = xa[i__] + c__ * d__[i__];
    }
    indic = 4;
    (*simul)(&indic, n, &xb[1], &fb, &gb[1], &izs[1], &rzs[1], &dzs[1], id);
/*     next line added by Serge to avoid Inf and Nan's (04/2007) */
    if (!R_FINITE(fb) && vff_(n, &gb[1]) != 1) {
	indic = -1;
    }
/*              test sur indic */
    if (indic > 0) {
	goto L185;
    }
    if (indic < 0) {
	goto L183;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = xb[i__];
	g[i__] = gb[i__];
    }
    goto L250;
L183:
    stepbd = step;
    ial = 1;
    step /= 10.;
    if (stepbd > steplb) {
	goto L170;
    }
    goto L240;
/*             stockage si c'est la plus petite valeur */
L185:
    isfv = min(2,isfv);
    if (fb > *f) {
	goto L220;
    }
    if (fb < *f) {
	goto L200;
    }
    gl1 = 0.;
    gl2 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = scale[i__] * g[i__];
	gl1 += d__1 * d__1;
/* Computing 2nd power */
	d__1 = scale[i__] * gb[i__];
	gl2 += d__1 * d__1;
    }
    if (gl2 >= gl1) {
	goto L220;
    }
L200:
    isfv = 3;
    *f = fb;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = xb[i__];
	g[i__] = gb[i__];
    }
/*               calcul de la derivee directionnelle */
L220:
    dgb = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dgb += gb[i__] * d__[i__];
    }
    if (*iprint < 3) {
	goto L231;
    }
    /* s = fb - fa; */
L231:
    if (fb - fa <= c__ * .1 * dga) {
	goto L280;
    }
    ial = 0;
/*               iteration terminee si le pas est minimum */
    if (step > steplb) {
	goto L270;
    }
L240:
    if (isfv >= 2) {
	goto L110;
    }
/*               ici, tout est termine */
L250:
    *acc = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*acc += g[i__] * g[i__];
    }
    *niter = itr;
    *nsim = nfun;
    return 0;
/*               interpolation cubique */
L270:
    stepbd = step;
    c__ = gmin + dgb - (fb - fmin) * 3. / safe_zero(step);
    if (c__ == 0.) {
	goto L250;
    }
    cc = fabs(c__) - gmin * (dgb / safe_zero(fabs(c__)));
    cc = sqrt((fabs(c__))) * sqrt((max(0.,cc)));
    c__ = (c__ - gmin + cc) / safe_zero(dgb - gmin + cc + cc);
    step *= max(.1,c__);
    goto L170;
/*               ceci est un pas de descente */
L280:
    if (ial == 0) {
	goto L285;
    }
    if (stepbd > steplb) {
	goto L285;
    }
    goto L240;
L285:
    stepbd -= step;
    stmin = c__;
    fmin = fb;
    gmin = dgb;
/*               extrapolation */
    step = stmin * 9.;
    if (stepbd > 0.) {
	step = stepbd * .5;
    }
    c__ = dga + dgb * 3. - (fb - fa) * 4. / safe_zero(stmin);
    if (c__ > 0.) {
/* Computing MIN */
/* Computing MAX */
	d__3 = 1., d__4 = -dgb / safe_zero(c__);
	d__1 = step, d__2 = stmin * max(d__3,d__4);
	step = min(d__1,d__2);
    }
    if (dgb < dga * .7) {
	goto L170;
    }
/*                 recherche lineaire terminee, test de convergence */
    isfv = 4 - isfv;
    if (stmin + step <= steplb) {
	goto L240;
    }
/*                 formule de bfgs */
    ir = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xa[i__] = xb[i__];
	xb[i__] = ga[i__];
	d__[i__] = gb[i__] - ga[i__];
	ga[i__] = gb[i__];
    }
    d__1 = 1. / safe_zero(dga);
    majour_(&h__[1], &xb[1], &w[1], n, &d__1, &ir, &c__1, &c_b32);
    ir = -ir;
    d__1 = 1. / safe_zero(stmin * (dgb - dga));
    majour_(&h__[1], &d__[1], &d__[1], n, &d__1, &ir, &c__1, &c_b32);
/* ww edits */
/*     write(*,*) (h(kk), kk=1,(n*(n+1))/2) */
/*                  test du rang de la nouvelle matrice */
    if (ir < *n) {
	goto L250;
    }
/*                  nouvelle iteration */
    dff = fa - fb;
    fa = fb;
    goto L130;
} /* n1qn1a_ */

/* Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) INRIA */

/* This file must be used under the terms of the CeCILL. */
/* This source file is licensed as described in the file COPYING, which */
/* you should have received as part of this distribution.  The terms */
/* are also available at */
/* http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt */

/* Subroutine */ int majour_(double *hm, double *hd, double *dd, 
	int *n, double *hno, int *ir, int *indic, double *
	eps)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Local variables */
    double b;
    int i__, j;
    double r__, y, gm;
    int ll, mm, np;
    double del, hml, hon, honm;
    int iplus;


    /* Parameter adjustments */
    --hm;
    --dd;
    --hd;

    /* Function Body */
    if (*n == 1) {
	goto L100;
    }

    np = *n + 1;
    if (*hno > 0.) {
	goto L99;
    }

    if (*hno == 0.) {
	goto L999;
    }
    if (*ir == 0) {
	goto L999;
    }
    hon = 1. / safe_zero(*hno);
    ll = 1;
    if (*indic == 0) {
	goto L1;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (hm[ll] == 0.) {
	    goto L2;
	}
/* Computing 2nd power */
	d__1 = dd[i__];
	hon += d__1 * d__1 / safe_zero(hm[ll]);
	ll = ll + np - i__;
    }
L2:
    goto L3;

L1:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd[i__] = hd[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iplus = i__ + 1;
	del = dd[i__];
	if (hm[ll] > 0.) {
	    goto L6;
	}
	dd[i__] = 0.;
	ll = ll + np - i__;
	goto L5;
L6:
/* Computing 2nd power */
	d__1 = del;
	hon += d__1 * d__1 / safe_zero(hm[ll]);
	if (i__ == *n) {
	    goto L7;
	}
	i__2 = *n;
	for (j = iplus; j <= i__2; ++j) {
	    ++ll;
	    dd[j] -= del * hm[ll];
	}
L7:
	++ll;
    }
L5:

L3:
    if (*ir <= 0) {
	goto L9;
    }
    if (hon > 0.) {
	goto L10;
    }
    if (*indic - 1 <= 0) {
	goto L99;
    } else {
	goto L11;
    }
L9:
    hon = 0.;
    *ir = -(*ir) - 1;
    goto L11;
L10:
    hon = *eps / safe_zero(*hno);
    if (*eps == 0.) {
	--(*ir);
    }
L11:
    mm = 1;
    honm = hon;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = np - i__;
	ll -= i__;
	if (hm[ll] != 0.) {
/* Computing 2nd power */
	    d__1 = dd[j];
	    honm = hon - d__1 * d__1 / safe_zero(hm[ll]);
	}
	dd[j] = hon;
	hon = honm;
    }
    goto L13;

L99:
    mm = 0;
    honm = 1. / safe_zero(*hno);
L13:
    ll = 1;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iplus = i__ + 1;
	del = hd[i__];
	if (hm[ll] > 0.) {
	    goto L14;
	}
	if (*ir > 0) {
	    goto L15;
	}
	if (*hno < 0.) {
	    goto L15;
	}
	if (del == 0.) {
	    goto L15;
	}
	*ir = 1 - *ir;
/* Computing 2nd power */
	d__1 = del;
	hm[ll] = d__1 * d__1 / safe_zero(honm);
	if (i__ == *n) {
	    goto L999;
	}
	i__2 = *n;
	for (j = iplus; j <= i__2; ++j) {
	    ++ll;
	    hm[ll] = hd[j] / safe_zero(del);
	}
	goto L999;
L15:
	hon = honm;
	ll = ll + np - i__;
	goto L98;
L14:
	hml = del / safe_zero(hm[ll]);
	if (mm <= 0) {
	    goto L17;
	} else {
	    goto L18;
	}
L17:
	hon = honm + del * hml;
	goto L19;
L18:
	hon = dd[i__];
L19:
	r__ = hon / safe_zero(honm);
	hm[ll] *= r__;
	if (r__ == 0.) {
	    goto L20;
	}
	if (i__ == *n) {
	    goto L20;
	}
	b = hml / safe_zero(hon);
	if (r__ > 4.) {
	    goto L21;
	}
	i__2 = *n;
	for (j = iplus; j <= i__2; ++j) {
	    ++ll;
	    hd[j] -= del * hm[ll];
	    hm[ll] += b * hd[j];
	}
	goto L23;
L21:
	gm = honm / safe_zero(hon);
	i__2 = *n;
	for (j = iplus; j <= i__2; ++j) {
	    ++ll;
	    y = hm[ll];
	    hm[ll] = b * hd[j] + y * gm;
	    hd[j] -= del * y;
	}
L23:
	honm = hon;
	++ll;
L98:
	;
    }

L20:
    if (*ir < 0) {
	*ir = -(*ir);
    }
    goto L999;
L100:
/* Computing 2nd power */
    d__1 = hd[1];
    hm[1] += *hno * (d__1 * d__1);
    *ir = 1;
    if (hm[1] > 0.) {
	goto L999;
    }
    hm[1] = 0.;
    *ir = 0;
L999:
    return 0;
} /* majour_ */

