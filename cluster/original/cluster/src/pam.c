/*
 * PAM := Partitioning Around Medoids
 *
 * original Id: pam.f,v 1.16 2003/06/03 13:40:56 maechler translated by
 * f2c (version 20031025) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <float.h>

#include <Rmath.h>
#include <Rinternals.h>

#include <R_ext/Print.h>/* for diagnostics */
#include <R_ext/Utils.h>/* for interrupting */

#include "cluster.h"
#include "ind_2.h"

// carries out a clustering using the k-medoid approach

SEXP cl_Pam(SEXP k_, SEXP n_,
	    SEXP do_diss_, /* == !diss;  if true, compute distances from x (= x_or_diss);
			      otherwise distances provided by x_or_diss */
	    SEXP x_or_diss,// this "is"  if(do_diss) "x[]" (n x p) else "dys[]"
	    SEXP all_stats_, // all_stats == !cluster.only
	    SEXP medoids, // NULL or integer(k) subset {1:n}
	    SEXP do_swap_, SEXP trace_lev_,
	    SEXP keep_diss_, SEXP pam_once_,

	    // the next 3 are only needed  if(do_diss)
	    SEXP val_md, SEXP j_md, // "md" := [m]issing [d]ata
	    SEXP diss_kind_) // = 1 ("euclidean")  or 2 ("manhattan")
{
    const int kk = asInteger(k_), n = asInteger(n_),
	pam_once = asInteger(pam_once_),
	trace_lev = asInteger(trace_lev_);
    const Rboolean all_stats = asLogical(all_stats_)
	, med_given = (medoids != R_NilValue) // if true, med[] contain initial medoids
	, do_diss = asLogical(do_diss_)
	, do_swap = asLogical(do_swap_)
	, keep_diss = asLogical(keep_diss_) // only  if(keep_diss)  return dys[] ..
	, do_syl = all_stats && (1 < kk && kk < n);


#ifdef once_we_get_n_from_args
    int n, p = NA_INTEGER;
    if (do_diss) { // <-- was 'jdyss != 1' i.e.  jdyss == 0
	SEXP dims = getAttrib(x_or_diss, R_DimSymbol);
	n = INTEGER(dims)[0];
	p = INTEGER(dims)[1];
    } else {
	n = asInteger(getAttrib(x_or_diss, install("Size")));
    }
#endif

    int i, nhalf; // nhalf := #{distances}+1 = length(dys)
    double s;
    if (n % 2 == 0) { // avoid overflow of n * (n - 1)
	nhalf = n / 2 * (n - 1) + 1;
    } else {
	nhalf = (n - 1) / 2 * n + 1;
    }

    int   *nsend = (int*) R_alloc(n, sizeof(int))
	, *nelem = (int*) R_alloc(all_stats ? n : 1, sizeof(int)) /* Rboolean */
	, *nrepr = (int*) R_alloc(n, sizeof(int))
	, *med
	;
    double
	*radus = (double*) R_alloc( n, sizeof(double)),
	*damer = (double*) R_alloc( n, sizeof(double)),
        *separ = (double*) R_alloc(kk, sizeof(double));
    int clusinf_dim1 = kk;

    if(med_given) {
	if(TYPEOF(medoids) != INTSXP || LENGTH(medoids) != kk)
	    error(_("Invalid 'medoids'"));
	PROTECT(medoids = duplicate(medoids));
    } else {
	PROTECT(medoids = allocVector(INTSXP, kk));
    }
    med = INTEGER(medoids);

    SEXP nms,
	ans = PROTECT(allocVector(VECSXP, keep_diss ? 9 : 9-1));
    setAttrib(ans, R_NamesSymbol,
	      nms = allocVector(STRSXP, keep_diss ? 9 : 9-1));
    int nprot = 2; // <- ++ for each PROTECT() below
    SEXP dys_, avsyl_, obj_, clu_, clusinf_, sylinf_, nisol_,ttsyl_;

    // these are only used  if(do_diss) :
    double *valmd; int *jtmd; int *diss_kind;
    if (do_diss) { // <-- was 'jdyss != 1' i.e.  jdyss == 0
	PROTECT(dys_ = allocVector(REALSXP, nhalf)); nprot++;
	valmd = REAL(val_md);
	jtmd = INTEGER(j_md);
	diss_kind = INTEGER(diss_kind_); // = 1 ("euclidean")  or 2 ("manhattan")
    } else {
	dys_ = x_or_diss; // a pointer to the same thing
    }
    // Creating the SEXPs as list components, so they are auto-PROTECTed:
    SET_STRING_ELT(nms, 0, mkChar("clu"));
    SET_VECTOR_ELT(ans, 0,         clu_ = allocVector(INTSXP, n));
    SET_STRING_ELT(nms, 1, mkChar("med")); SET_VECTOR_ELT(ans, 1, medoids);
    SET_STRING_ELT(nms, 2, mkChar("silinf"));
    if(do_syl)
    SET_VECTOR_ELT(ans, 2,         sylinf_ = all_stats ? allocMatrix(REALSXP, n, 4)
		   : allocVector(REALSXP, 1));
    SET_STRING_ELT(nms, 3, mkChar("obj"));
    SET_VECTOR_ELT(ans, 3,         obj_ = allocVector(REALSXP, 2));
    SET_STRING_ELT(nms, 4, mkChar("isol")); // 'isolation' {0,1,2}:  L=1, L*=2 - clusters
    SET_VECTOR_ELT(ans, 4,         nisol_ = allocVector(INTSXP, all_stats ? kk : 1));
    SET_STRING_ELT(nms, 5, mkChar("clusinf"));
    SET_VECTOR_ELT(ans, 5,         clusinf_ = all_stats ? allocMatrix(REALSXP, kk, 5)
		   : allocVector(REALSXP, 1));
    SET_STRING_ELT(nms, 6, mkChar("avsil"));
    SET_VECTOR_ELT(ans, 6,         avsyl_ = allocVector(REALSXP, n));
    SET_STRING_ELT(nms, 7, mkChar("ttsil"));
    if(do_syl)
    SET_VECTOR_ELT(ans, 7,         ttsyl_ = allocVector(REALSXP, 1));
    if(keep_diss) {
	SET_STRING_ELT(nms, 8, mkChar("dys"));	SET_VECTOR_ELT(ans, 8, dys_);
    }

    int *ncluv = INTEGER(clu_),
	*nisol = INTEGER(nisol_);
    double
	*dys    = REAL(dys_),
	*avsyl  = REAL(avsyl_),
	*obj    = REAL(obj_),
	*clusinf= REAL(clusinf_);

    if (do_diss) { // <-- was 'jdyss != 1' i.e.  jdyss == 0
	double *x = REAL(x_or_diss);
	SEXP dims = getAttrib(x_or_diss, R_DimSymbol);
	int p = INTEGER(dims)[1];
	if(trace_lev)
	    Rprintf("C pam(): computing %d dissimilarities from  %d x %d  matrix: ",
		    nhalf, n, p);
	int jhalt = dysta((int*)&n, &p, x, dys, diss_kind, jtmd, valmd); /* --> ./dysta.c */
	if (jhalt != 0) {
	    if(trace_lev) Rprintf(" dysta()-error: jhalt=%d\n", jhalt);
	    UNPROTECT(nprot);
	    return ScalarInteger(jhalt); // i.e., integer error code instead of a named list
	}
	// else
	if(trace_lev) Rprintf("[Ok]\n");
    }

    /* s := max( dys[.] ), the largest distance */
    for (i = 1, s = 0.; i < nhalf; ++i) /* dys[0] == 0. not used here */
	if (s < dys[i])
	    s = dys[i];

    /* FIXME: work with med[] = (i_1, i_2, ..., i_k)
     * ----- instead nrepr[] = (b_1, ... b_n)   b_i in {0,1} */
    for (i = 0; i < n; ++i)
	nrepr[i] = 0;
    if(med_given) { /* if true, med[] contain initial medoids */

	/* for the moment, translate these to nrepr[] 0/1 :
	 * not assuming that the med[] indices are sorted */
	for (int k = 0; k < kk; k++)
	    nrepr[med[k] - 1] = 1;
    }

/*     Build + Swap [but no build if(med_given); swap only if(do_swap) : */

    bswap(kk, n, nrepr,			// <-  3
	  med_given, do_swap, trace_lev,// <-  6
	  radus, damer, avsyl, 		// <-  9
	  dys, s, obj,			// <- 12
	  pam_once);

    if(trace_lev) Rprintf("end{bswap()}, ");
/*     Compute Clustering & STATs if(all_stats): */
    cstat(kk, n, nsend, nrepr, all_stats,
	  radus, damer, avsyl, separ, &s, dys, ncluv, nelem, med, nisol);
    if(trace_lev) Rprintf("end{cstat()}\n");
    if(all_stats) {
	for (int k = 0; k < kk; ++k) {
	    clusinf[k]=		(double)       nrepr[k];
	    clusinf[k + clusinf_dim1]	     = radus[k];
	    clusinf[k + (clusinf_dim1 << 1)] = avsyl[k];
	    clusinf[k + clusinf_dim1 * 3]    = damer[k];
	    clusinf[k + (clusinf_dim1 << 2)] = separ[k];
	}
	if (do_syl) { // Compute Silhouette info :
	    double
		*ttsyl  = REAL(ttsyl_),
		*sylinf = REAL(sylinf_);

	    dark(kk, n, ncluv, dys, s,
		 // -->
		 nsend, nelem, nrepr, radus, damer, avsyl, ttsyl, sylinf);
	}
    }
    UNPROTECT(nprot);
    return ans;
} /* cl_Pam */



/* -----------------------------------------------------------
     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
*/
void bswap(int kk, int n, int *nrepr,
	   Rboolean med_given, Rboolean do_swap, int trace_lev,
	   /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
	   double *dysma, double *dysmb, double *beter,
	   const double dys[], double s, double *obj, int pamonce)
{
    int i, j, ij, k,h, dig_n;
    double sky;

    /* Parameter adjustments */
    --nrepr;
    --beter;

    --dysma; --dysmb;

    if(trace_lev) Rprintf("pam()'s bswap(*, s=%g, pamonce=%d): ", s, pamonce);

    s = s * 1.1 + 1.;// larger than all dys[]  (but DBL_MAX is too large)


/* IDEA: when n is large compared to k (= kk),
 * ----  rather use a "sparse" representation:
 * instead of boolean vector nrepr[] , use  ind_repr <- which(nrepr) !!
 */
    for (i = 1; i <= n; ++i)
	dysma[i] = s;

    if(med_given) {
	if(trace_lev) Rprintf("medoids given; ");

	/* compute dysma[] : dysma[j] = D(j, nearest_representative) */
	for (i = 1; i <= n; ++i) {
	    if (nrepr[i] == 1)
		for (j = 1; j <= n; ++j) {
		    ij = ind_2(i, j);
		    if (dysma[j] > dys[ij])
			dysma[j] = dys[ij];
		}
	}
    }
    else { // ====== first algorithm: BUILD. ======================

	if(trace_lev) Rprintf("build %d medoids:\n", kk);

	/* find  kk  representatives  aka medoids :  */

	for (k = 1; k <= kk; ++k) {

	    R_CheckUserInterrupt();

	    /* compute beter[i] for all non-representatives:
	     * also find ammax := max_{..} and nmax := argmax_i{beter[i]} ... */
	    int nmax = -1; /* -Wall */
	    double ammax, cmd;
	    ammax = 0.;
	    for (i = 1; i <= n; ++i) {
		if (nrepr[i] == 0) {
		    beter[i] = 0.;
		    for (j = 1; j <= n; ++j) {
			cmd = dysma[j] - dys[ind_2(i, j)];
			if (cmd > 0.)
			    beter[i] += cmd;
		    }
		    if (ammax <= beter[i]) {
			/*  does < (instead of <= ) work too? -- NO! */
			ammax = beter[i];
			nmax = i;
		    }
		}
	    }

	    nrepr[nmax] = 1;/* = .true. : found new representative */
	    if (trace_lev >= 2)
		Rprintf("    new repr. %d\n", nmax);

	    /* update dysma[] : dysma[j] = D(j, nearest_representative) */
	    for (j = 1; j <= n; ++j) {
		ij = ind_2(nmax, j);
		if (dysma[j] > dys[ij])
		    dysma[j] = dys[ij];
	    }
	}
	/* output of the above loop:  nrepr[], dysma[], ... */
    }

    if(trace_lev) /* >= 2 (?) */ {
	dig_n = 1+(int) log10(n); // = 1 + floor(log10(n))
	Rprintf("  after build: medoids are");
	for (i = 1; i <= n; ++i)
	    if(nrepr[i] == 1) Rprintf(" %*d", dig_n, i);
	if(trace_lev >= 3) {
	    Rprintf("\n  and min.dist dysma[1:n] are\n");
	    for (i = 1; i <= n; ++i) {
		Rprintf(" %6.3g", dysma[i]);
		if(i % 10 == 0) Rprintf("\n");
	    }
	    if(n % 10 != 0) Rprintf("\n");
	} else Rprintf("\n");
    } else dig_n = 1;// -Wall

    sky = 0.;
    for (j = 1; j <= n; ++j)
	sky += dysma[j];
    obj[0] = sky / n;

    if (do_swap && (kk > 1 || med_given)) {

	double dzsky;
	int hbest = -1, nbest = -1, kbest= -1; // -Wall
	int *medoids, *clustmembership, *clustmembershipb;
	double *fvect;
	if(pamonce) {
	    // +1 --> use 1-based indices (as R)
	    medoids = (int*) R_alloc(kk+1, sizeof(int));
	    clustmembership = (int*) R_alloc(n+1, sizeof(int));
	    fvect = (double*) R_alloc(n+1, sizeof(double));
	    for (int k = 1, i = 1; i <= n; ++i) {
		if (nrepr[i]) {
		    medoids[k] = i;
		    k++;
		}
	    }
	} else { // -Wall :
	    clustmembership = medoids = (int*) NULL;
	    clustmembershipb = medoids = (int*) NULL;
	    fvect = (double*) NULL;
	}
	int *best_h = NULL; double *best_d = NULL;
	if (pamonce == 4 || pamonce == 5) { // Schubert and Rousseeuw 2019 FastPAM2
	    best_h = (int*) R_alloc(kk+1, sizeof(int));
	    best_d = (double*) R_alloc(kk+1, sizeof(double));
	}
	// Improvement part of Schubert and Rousseeuw 2021 FasterPAM
	double *remove_loss = NULL;
	if (pamonce >= 3 && pamonce <= 6) {
	    remove_loss = (double*) R_alloc(kk+1, sizeof(double));
	    clustmembershipb = (int*) R_alloc(n+1, sizeof(int));
	}

/* ====== second algorithm: SWAP. ====== */

	/* Hmm: In the following, we RE-compute dysma[];
	 *      don't need it first time; then only need *update* after swap */

	// Used by FasterPAM only (Schubert and Rousseeuw 2021, pamonce = 6)
	int lastswap = 0;

/*--   Loop : */
    L60:
	if(pamonce == 0) { // original algorithm
	    for (j = 1; j <= n; ++j) {
		/*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
		 *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
		dysma[j] = s;
		dysmb[j] = s;
		for (i = 1; i <= n; ++i) {
		    if (nrepr[i]) {
			ij = ind_2(i, j);
			if (dysma[j] > dys[ij]) {
			    dysmb[j] = dysma[j];
			    dysma[j] = dys[ij];
			} else if (dysmb[j] > dys[ij]) {
			    dysmb[j] = dys[ij];
			}
		    }
		}
	    }
	} else if (pamonce >= 3 && pamonce <= 6) { // Schubert and Rousseeuw 2019 FastPAM1/2
	    for (j = 1; j <= n; ++j) {
		/*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
		 *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
		double dysmaj = s, dysmbj = s;
		for(k = 1; k <= kk; k++) {
		    double dysij = dys[ind_2(medoids[k], j)];
		    if (dysmaj > dysij) {
			clustmembershipb[j] = clustmembership[j];
			dysmbj = dysmaj;
			clustmembership[j] = k; // Use medoid k, not i
			dysmaj = dysij;
		    } else if (dysmbj > dysij) {
			clustmembershipb[j] = k;
			dysmbj = dysij;
		    }
		}
		dysma[j] = dysmaj;
		dysmb[j] = dysmbj;
	    }
	} else { // pamonce == 1 or == 2 :
	    for (j = 1; j <= n; ++j) {
		/*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
		 *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
		dysma[j] = s;
		dysmb[j] = s;
		for(k = 1; k <= kk; k++) {
		    i = medoids[k];
		    ij = ind_2(i, j);
		    if (dysma[j] > dys[ij]) {
			//store cluster membership
			clustmembership[j] = i;
			dysmb[j] = dysma[j];
			dysma[j] = dys[ij];
		    } else if (dysmb[j] > dys[ij]) {
			dysmb[j] = dys[ij];
		    }
		}
	    }
	}

// FasterPAM improvement: do not recompute the nearest, but update it in the end as necessary
L60b:
	dzsky = 1.; /* 1 is arbitrary > 0; only dzsky < 0 matters in the end */

	if(pamonce == 0) { // original algorithm
	    for (h = 1; h <= n; ++h) if (!nrepr[h]) {
		    R_CheckUserInterrupt();
		    for (i = 1; i <= n; ++i) if (nrepr[i]) {
			    double dz = 0.;
			    /* dz := T_{ih} := sum_j C_{jih}  [p.104] : */
			    for (j = 1; j <= n; ++j) { /* if (!nrepr[j]) { */
				int hj = ind_2(h, j);
				ij = ind_2(i, j);
				if (dys[ij] == dysma[j]) {
				    double small = dysmb[j] > dys[hj] ? dys[hj] : dysmb[j];
				    dz += (- dysma[j] + small);
				} else if (dys[hj] < dysma[j]) /* 1c. */
				    dz += (- dysma[j] + dys[hj]);
			    }
			    if (dzsky > dz) {
				dzsky = dz; /* dzsky := min_{i,h} T_{i,h} */
				hbest = h;
				nbest = i;
			    }
			}
		}
	} else if (pamonce == 1 || pamonce == 2) {
	    for(k = 1; k <= kk; k++) {
		R_CheckUserInterrupt();
		i=medoids[k];
		double removeCost = 0.;
		//Compute cost for removing the medoid
		for (j = 1; j <= n; ++j) {
		    if(clustmembership[j] == i) {
			removeCost+=(dysmb[j]-dysma[j]);
			fvect[j]=dysmb[j];
		    }
		    else{
			fvect[j]=dysma[j];
		    }
		}

		if (pamonce == 1) {
		    // Now check possible new medoids h
		    for (h = 1; h <= n; ++h) if (!nrepr[h]) {
			    double addGain = removeCost;
			    // Compute gain of adding h as a medoid:
			    for (j = 1; j <= n; ++j) {
				int hj = ind_2(h, j);
				if(dys[hj] < fvect[j])
				    addGain += (dys[hj]-fvect[j]);
			    }
			    if (dzsky > addGain) {
				dzsky = addGain; /* dzsky := min_{i,h} T_{i,h} */
				hbest = h;
				nbest = i;
				kbest = k;
			    }
			}

		} else { // pamonce == 2 :

		    // Now check possible new medoids h
		    for (h = 1; h <= n; ++h) if (!nrepr[h]) {
			    double addGain = removeCost - fvect[h]; // - fvect[h] since dys[h,h]=0;
			    // Compute gain of adding h as a medoid:
			    int ijbase = (h-2)*(h-1)/2;
			    for (j = 1; j < h; ++j) {
				int hj = ijbase+j;
				if(dys[hj] < fvect[j])
				    addGain += (dys[hj]-fvect[j]);
			    }
			    ijbase += h;// = (h-2)*(h-1)/2 + h
			    for (j = h+1; j <= n; ++j) {
				ijbase += j-2;
				if(dys[ijbase] < fvect[j])
				    addGain += (dys[ijbase]-fvect[j]);
			    }
			    if (dzsky > addGain) {
				dzsky = addGain; /* dzsky := min_{i,h} T_{i,h} */
				hbest = h;
				nbest = i;
				kbest = k;
			    }
			}
		}
	    }
	} else if (pamonce >= 3 && pamonce <= 6) { // Schubert and Rousseeuw variants
	    // cost of removing current medoids
	    for(k = 1; k <= kk; k++)
		remove_loss[k] = 0;
	    for (j = 1; j <= n; ++j)
		remove_loss[clustmembership[j]] += dysmb[j] - dysma[j];

	    if (pamonce == 3) { // Schubert and Rousseeuw 2019 FastPAM1 -- improved 2021 version
		// Now check possible new medoids h
		for (h = 1; h <= n; ++h) if (!nrepr[h]) {
		    R_CheckUserInterrupt();
		    for(k = 1; k <= kk; k++)
			beter[k] = 0;
		    double acc = 0;
		    // Compute gain of substituting h for each other medoid:
		    for (j = 1; j <= n; ++j) {
			double dist_h = dys[ind_2(h, j)]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    for(k = 1; k <= kk; k++) {
			double loss = beter[k] + acc + remove_loss[k];
			if (loss < dzsky) {
			    dzsky = loss;
			    hbest = h;
			    nbest = medoids[k];
			    kbest = k;
			}
		    }
		}
	    } else if (pamonce == 4) { // Schubert and Rousseeuw 2019 FastPAM2
		for(k = 1; k <= kk; k++) {
		    best_d[k] = 1.; // arbitrary > 0
		}
		// Now check possible new medoids h
		for (h = 1; h <= n; ++h) if (!nrepr[h]) {
		    R_CheckUserInterrupt();
		    for(k = 1; k <= kk; k++)
			beter[k] = 0;
		    // Compute gain of substituting h for each other medoid:
		    double acc = 0;
		    for (j = 1; j <= n; ++j) {
			double dist_h = dys[ind_2(h, j)]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    for(k = 1; k <= kk; k++) {
			double loss = beter[k] + acc + remove_loss[k];
			if (loss < best_d[k]) {
			    best_d[k] = loss;
			    best_h[k] = h;
			}
		    }
		}
		// Pick one best medoid at the end (additional swaps come later).
		for(k = 1; k <= kk; k++) {
		    if (best_d[k] < dzsky) {
			dzsky = best_d[k];
			kbest = k;
		    }
		}
		hbest = best_h[kbest];
		nbest = medoids[kbest];
	    } else if (pamonce == 5) { // Schubert and Rousseeuw 2019 FastPAM2,
		// with linearized memory access as in Reynolds 2
		for(k = 1; k <= kk; k++) {
		    best_d[k] = 1.; // arbitrary > 0
		}
		// Now check possible new medoids h
		for (h = 1; h <= n; ++h) if (!nrepr[h]) {
		    R_CheckUserInterrupt();
		    for(k = 1; k <= kk; k++)
			beter[k] = 0;
		    double acc = -dysma[h];
		    beter[clustmembership[h]] = dysma[h] - dysmb[h];
		    // Compute gain of substituting h for each other medoid:
		    int ijbase = (h-2)*(h-1)/2;
		    for (j = 1; j < h; ++j) {
			double dist_h = dys[ijbase+j]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    ijbase += h;// = (h-2)*(h-1)/2 + h
		    for (j = h+1; j <= n; ++j) {
			ijbase += j-2;
			double dist_h = dys[ijbase]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    for(k = 1; k <= kk; k++) {
			double loss = beter[k] + acc + remove_loss[k];
			if (loss < best_d[k]) {
			    best_d[k] = loss;
			    best_h[k] = h;
			}
		    }
		}
		// Pick one best medoid at the end (additional swaps come later).
		for(k = 1; k <= kk; k++) {
		    if (best_d[k] < dzsky) {
			dzsky = best_d[k];
			kbest = k;
		    }
		}
		hbest = best_h[kbest];
		nbest = medoids[kbest];
	    } else { // Schubert and Rousseeuw, 2021, FasterPAM
		// Check possible new medoids h
		// Beginning at the last swap performed
		int h2 = 1;
		for (h = lastswap + 1; h2 <= n; ++h, ++h2) {
		    if (h == n+1) h = 1; // wrap around
		    if (nrepr[h]) continue;
		    R_CheckUserInterrupt();
		    for(k = 1; k <= kk; k++)
			beter[k] = 0;
		    double acc = -dysma[h];
		    beter[clustmembership[h]] += dysma[h] - dysmb[h];
		    // Compute gain of substituting h for each other medoid:
		    int ijbase = (h-2)*(h-1)/2;
		    for (j = 1; j < h; ++j) {
			double dist_h = dys[ijbase+j]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    ijbase += h;// = (h-2)*(h-1)/2 + h
		    for (j = h+1; j <= n; ++j) {
			ijbase += j-2;
			double dist_h = dys[ijbase]; // New medoid
			int memb = clustmembership[j]; // Medoid nr of nearest
			double distcur = dysma[j]; // Nearest
			double distsec = dysmb[j]; // Second nearest
			// New medoid is new best:
			if (dist_h < distcur) {
			    acc += dist_h - distcur;
			    beter[memb] += distcur - distsec;
			} else if (dist_h < distsec) {
			    beter[memb] += dist_h - distsec;
			}
		    }
		    for(k = 1; k <= kk; k++) {
			double loss = beter[k] + acc + remove_loss[k];
			if (loss < dzsky) {
			    dzsky = loss;
			    hbest = h;
			    nbest = medoids[k];
			    kbest = k;
			}
		    }
		    // Execute the first swap found.
		    if (dzsky < 0) {
			break;
		    }
		}
	    }
	} else {
	    Rprintf("Invalid pamonce value: %d\n", pamonce);
	    return;
	}

	if (dzsky < - 16*DBL_EPSILON * fabs(sky)) { // basically " < 0 ",
	    // but ' < 0 ' gave infinite loop, swapping the identical objects
	    // found an improving swap

	    if(trace_lev >= 2)
		Rprintf( "   swp new %*d <-> %*d old; decreasing diss. %7g by %g\n",
			 dig_n, hbest, dig_n, nbest, sky, dzsky);
	    nrepr[hbest] = 1;
	    nrepr[nbest] = 0;
	    if(pamonce)
		medoids[kbest]=hbest;
	    if(pamonce == 6)
		lastswap = hbest;

	    sky += dzsky;
	    if (pamonce == 3 || pamonce == 6) { // FasterPAM improvement, Schubert and Rousseeuw 2021
		// Update dysma, dysmb rather than recomputing it
		double newloss = 0;
		for (j = 1; j <= n; ++j) {
		    /*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
		     *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
		    double dnew = dys[ind_2(hbest, j)];
		    if (clustmembership[j] == kbest) {
			if (dnew < dysmb[j]) {
			    // Clustmembership already kbest
			    dysma[j] = dnew;
			} else {
			    double dysmaj = s, dysmbj = s;
			    for(k = 1; k <= kk; k++) {
				double dysij = dys[ind_2(medoids[k], j)];
				if (dysmaj > dysij) {
				    clustmembershipb[j] = clustmembership[j];
				    dysmbj = dysmaj;
				    //store cluster membership
				    clustmembership[j] = k; // Use medoid k, not i
				    dysmaj = dysij;
				} else if (dysmbj > dysij) {
				    clustmembershipb[j] = k; // Use medoid k, not i
				    dysmbj = dysij;
				}
			    }
			    dysma[j] = dysmaj;
			    dysmb[j] = dysmbj;
			}
		    } else {
			if (dnew < dysma[j]) {
			    clustmembershipb[j] = clustmembership[j]; // not kbest
			    dysmb[j] = dysma[j];
			    clustmembership[j] = kbest; // Use medoid k, not i
			    dysma[j] = dnew;
			} else if (clustmembershipb[j] == kbest) {
			    double dysmaj = s, dysmbj = s;
			    for(k = 1; k <= kk; k++) {
				double dysij = dys[ind_2(medoids[k], j)];
				if (dysmaj > dysij) {
				    clustmembershipb[j] = clustmembership[j];
				    dysmbj = dysmaj;
				    //store cluster membership
				    clustmembership[j] = k; // Use medoid k, not i
				    dysmaj = dysij;
				} else if (dysmbj > dysij) {
				    clustmembershipb[j] = k; // Use medoid k, not i
				    dysmbj = dysij;
				}
			    }
			    dysma[j] = dysmaj;
			    dysmb[j] = dysmbj;
			} else if (dnew < dysmb[j]) {
			    clustmembershipb[j] = kbest; // Use medoid k, not i
			    dysmb[j] = dnew;
			}
		    }
		    newloss += dysma[j];
		}
		sky = newloss; // more accurate than accumulating changes
		goto L60b;
	    }
	    // FastPAM2, Schubert and Rousseeuw 2019
	    if (pamonce == 4 || pamonce == 5)
		do {
		    best_d[kbest] = 1; // Deactivate
		    // Find next best:
		    dzsky = 1;
		    for(k = 1; k <= kk; k++) {
			if (best_d[k] < dzsky) {
			    dzsky = best_d[k];
			    kbest = k;
			}
		    }
		    hbest = best_h[kbest];
		    nbest = medoids[kbest];
		    if(dzsky >= - 16*DBL_EPSILON * fabs(sky))
			break;
		    // else try more :

		    // FIXME: duplicated code from above - update stats.
		    for (j = 1; j <= n; ++j) {
			/*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
			 *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
			dysma[j] = s;
			dysmb[j] = s;
			for(k = 1; k <= kk; k++) {
			    i = medoids[k];
			    ij = ind_2(i, j);
			    if (dysma[j] > dys[ij]) {
				//store cluster membership
				clustmembership[j] = k; // Use medoid k, not i
				dysmb[j] = dysma[j];
				dysma[j] = dys[ij];
			    } else if (dysmb[j] > dys[ij]) {
				dysmb[j] = dys[ij];
			    }
			}
		    }
		    // Recompute dzsky value, again.
		    dzsky = 0;
		    for (j = 1; j <= n; ++j) { /* if (!nrepr[j]) { */
			int hj = ind_2(hbest, j);
			ij = ind_2(nbest, j);
			if (dys[ij] == dysma[j]) {
			    dzsky += (- dysma[j] + fmin2(dysmb[j], dys[hj]) );
			} else if (dys[hj] < dysma[j]) /* 1c. */
			    dzsky += (- dysma[j] + dys[hj]);
		    }
		    if(dzsky >= - 16*DBL_EPSILON * fabs(sky))
			break;
		    // else try more :
		    if(trace_lev >= 2)
			Rprintf( "   fswp new %*d <-> %*d old; decreasing diss. %7g by %g\n",
				 dig_n, hbest, dig_n, nbest, sky, dzsky);
		    nrepr[hbest] = 1;
		    nrepr[nbest] = 0;
		    medoids[kbest]=hbest;
		    sky += dzsky;

		} while(1); // (pamonce = 4 or 5)

	    goto L60;
	}
    }
    obj[1] = sky / n;
} /* bswap */


/* -----------------------------------------------------------
 cstat(): Compute STATistics (numerical output) concerning each partition
*/
void cstat(int kk, int nn, int *nsend, int *nrepr, Rboolean all_stats,
	   double *radus, double *damer, double *avsyl, double *separ, double *s,
	   double *dys, int *ncluv, int *nelem, int *med, int *nisol)
{
    int j, k, ja, jk, nplac, ksmal = -1/* -Wall */;
    double ss = *s * 1.1 + 1.;

    /* Parameter adjustments */
    --ncluv;
    --nrepr;
    --nsend;

    /* nsend[j] := i,  where x[i,] is the medoid to which x[j,] belongs */
    for (j = 1; j <= nn; ++j) {
	if (nrepr[j] == 0) {
	    double dsmal = ss;
	    for (k = 1; k <= nn; ++k) {
		if (nrepr[k] == 1) {
		    int kj_ = ind_2(k, j);
		    if (dsmal > dys[kj_]) {
			dsmal = dys[kj_];
			ksmal = k;
		    }
		}
	    }
	    nsend[j] = ksmal;
	} else {
	    nsend[j] = j;
	}
    }
    /* ncluv[j] := k , the cluster number  (k = 1..kk) */
    jk = 1;
    nplac = nsend[1];
    for (j = 1; j <= nn; ++j) {
	ncluv[j] = 0;
	if (nsend[j] == nplac)
	    ncluv[j] = 1;
    }
    for (ja = 2; ja <= nn; ++ja) {
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    for (j = 2; j <= nn; ++j) {
		if (nsend[j] == nplac)
		    ncluv[j] = jk;
	    }
	    if (jk == kk)
		break;
	}
    }

    if(all_stats) { /*	   analysis of the clustering. */
	/* int numl; */
	--avsyl; // <-> [1]-indexing
	--damer;
	--med;
	--nelem;
	--nisol;
	--radus;
	--separ;
	for (k = 1; k <= kk; ++k) {
	    int ntt = 0, m = -1/* -Wall */;
	    double ttt = 0.;
	    radus[k] = -1.;
	    R_CheckUserInterrupt();
	    for (j = 1; j <= nn; ++j) {
		if (ncluv[j] == k) {
		    double djm;
		    ++ntt;
		    m = nsend[j];
		    nelem[ntt] = j;
		    djm = dys[ind_2(j, m)];
		    ttt += djm;
		    if (radus[k] < djm)
			radus[k] = djm;
		}
	    }
	    if(ntt == 0) error(_("pam(): Bug in C level cstat(), k=%d: ntt=0"), k);
	    avsyl[k] = ttt / ntt;
	    med[k] = m;
	}

	if (kk == 1) {
	    damer[1] = *s;
	    nrepr[1] = nn;
	    nisol[1] = 0;
	    separ[1] = 0.;
	    return;
	}
	/*  ELSE	  kk > 1 : */

	/* numl = number of L-clusters. */
	/* numl = 0; */
	for (k = 1; k <= kk; ++k) {
	    /*
	      identification of cluster k:
	      nelem= vector of object indices,
	      nel  = number of objects
	    */
	    int nel = 0;
	    R_CheckUserInterrupt();

	    for (j = 1; j <= nn; ++j) {
		if (ncluv[j] == k) {
		    ++nel;
		    nelem[nel] = j;
		}
	    }
	    nrepr[k] = nel;
	    if (nel == 1) {
		int nvn = nelem[1];
		damer[k] = 0.;
		separ[k] = ss;
		for (j = 1; j <= nn; ++j) {
		    if (j != nvn) {
			int mevj = ind_2(nvn, j);
			if (separ[k] > dys[mevj])
			    separ[k] = dys[mevj];
		    }
		}

		/* Is cluster k
		   1) an L-cluster	 or
		   2) an L*-cluster ? */
		/* if (separ[k] == 0.) */
		/*     ++numl; */

	    }
	    else { /*	       nel != 1 : */
		double dam = -1., sep = ss;
		Rboolean kand = TRUE;
		for (ja = 1; ja <= nel; ++ja) {
		    int jb, nvna = nelem[ja];
		    double aja = -1., ajb = ss;
		    for (jb = 1; jb <= nn; ++jb) {
			int jndz = ind_2(nvna, jb);
			if (ncluv[jb] == k) {
			    if (aja < dys[jndz])
				aja = dys[jndz];
			} else {
			    if (ajb > dys[jndz])
				ajb = dys[jndz];
			}
		    }
		    if (kand && aja >= ajb)
			kand = FALSE;
		    if (dam < aja)
			dam = aja;
		    if (sep > ajb)
			sep = ajb;
		}
		separ[k] = sep;
		damer[k] = dam;
		if (kand) {
		    /* ++numl; */
		    if (dam >= sep) /*	L-cluster */
			nisol[k] = 1;
		    else/*		L*-cluster */
			nisol[k] = 2;
		    continue /* k */;
		}
	    }
	    /* nel = 1 or (!kand) : */
	    nisol[k] = 0;

	}/* for(k) */

    } /* all_stats */

} /* cstat */

/* -----------------------------------------------------------
     Compute Silhouette Information :
 */
void dark(
    // input:
    int kk, int nn, const int ncluv[], const double dys[], double s,
    // output:
    int *nsend, int *nelem, int *negbr,
    double *syl, double *srank, double *avsyl, double *ttsyl,
    double *sylinf)
{
    int k, nsylr;
    /* pointers to sylinf[] columns -- sylinf[nn, 4] : */
    double *sylinf_2, *sylinf_3, *sylinf_4;
    sylinf_2 = sylinf	+ nn;
    sylinf_3 = sylinf_2 + nn;
    sylinf_4 = sylinf_3 + nn;

    /* Parameter adjustments */
    --avsyl;
    --ncluv;

    nsylr = 0;
    *ttsyl = 0.;
    for (k = 1; k <= kk; ++k) {
	/* nelem[0:(ntt-1)] := indices (1-based) of obs. in cluster k : */
	int j,l, ntt = 0;
	for (j = 1; j <= nn; ++j) {
	    if (ncluv[j] == k) {
		nelem[ntt] = j;
		++ntt;
	    }
	}

	for (j = 0; j < ntt; ++j) {/* (j+1)-th obs. in cluster k */
	    int k_, nj = nelem[j];
	    double dysb = s * 1.1 + 1.;
	    negbr[j] = -1;
	    /* for all clusters  k_ != k : */
	    for (k_ = 1; k_ <= kk; ++k_) if (k_ != k) {
		double db = 0.;
		int nbb = 0;
		for (l = 1; l <= nn; ++l) if (ncluv[l] == k_) {
		    ++nbb;
		    if (l != nj)
			db += dys[ind_2(nj, l)];
		}
		db /= nbb; /* now  db(k_) := mean( d[j, l]; l in C_{k_} ) */
		if (dysb > db) {
		    dysb = db;
		    negbr[j] = k_;
		}
	    }/* negbr[j] := arg max_{k_} db(k_) */
	    if (ntt > 1) {
		double dysa = 0.;
		for (l = 0; l < ntt; ++l) {
		    int nl = nelem[l];
		    if (nj != nl)
			dysa += dys[ind_2(nj, nl)];
		}
		dysa /= ntt - 1;
		if (dysa > 0.) {
		    if (dysb > 0.) {
			if (dysb > dysa)
			    syl[j] = 1. - dysa / dysb;
			else if (dysb < dysa)
			    syl[j] = dysb / dysa - 1.;
			else /* dysb == dysa: */
			    syl[j] = 0.;

			if (syl[j] < -1.)
			    syl[j] = -1.;
			else if (syl[j] > 1.)
			    syl[j] = 1.;

		    } else {
			syl[j] = -1.;
		    }
		}
		else /* dysa == 0 */ if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;
	    }
	    else { /*	  ntt == 1: */
		syl[j] = 0.;
	    }
	} /* for( j ) */
	avsyl[k] = 0.;
	if (ntt == 0) /* this can happen when medoids are user-specified !*/
	    continue; /* next k */

	for (j = 0; j < ntt; ++j) {
	    int lang=-1 /*Wall*/;
	    double symax = -2.;
	    for (l = 0; l < ntt; ++l) {
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = symax; /* = syl[lang] */
	    avsyl[k] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[k];
	avsyl[k] /= ntt;
	if (ntt == 1) {
	    sylinf  [nsylr] = (double) k;
	    sylinf_2[nsylr] = (double) negbr[0];
	    sylinf_3[nsylr] = 0.;
	    sylinf_4[nsylr] = (double) nelem[0];
	    ++nsylr;
	} else {
	    for (j = 0; j < ntt; ++j) {
		int lplac = nsend[j];
		sylinf	[nsylr] = (double) k;
		sylinf_2[nsylr] = (double) negbr[lplac];
		sylinf_3[nsylr] = srank[j];
		sylinf_4[nsylr] = (double) nelem[lplac];
		++nsylr;
	    }
	}
    } /* for (k) */
    *ttsyl /= nn;
} /* dark */
