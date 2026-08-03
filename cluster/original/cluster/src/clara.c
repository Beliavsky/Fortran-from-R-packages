
/*   Clustering LARge Applications
     ~		~~~   ~
     Clustering program based upon the k-medoid approach,
     and suitable for data sets of at least 100 objects.
     (for smaller data sets, please use program pam.)
 */

/* original Id: clara.f,v 1.10 2002/08/27 15:43:58 maechler translated by
 * f2c (version 20010821) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <math.h>

#include <R_ext/Print.h>/* for diagnostics */
#include <R_ext/Random.h>/* when R's RNG is used */
#include <R_ext/Utils.h>/* for interrupting */

#include "cluster.h"
#include "ind_2.h"

static void
dysta2(int nsam, int jpp, int *nsel,
       double *x, int n, double *dys, DISS_KIND diss_kind,
       int *jtmd, double *valmd, Rboolean has_NA, Rboolean *toomany_NA);

void cl_clara(int *n,  /* = number of objects */
	      int *jpp,/* = number of variables */
	      int *kk, /* = number of clusters, 1 <= kk <= n-1 */
	      double *x,	/* Input:  the data x[n, jpp] _rowwise_ (transposed)
				 * Output: the first `n' values are the `clustering'
				 *	   (integers in 1,2,..,kk) */
	      int *nran,	/* = #{random samples} drawn	   (= `samples' in R)*/
	      int *nsam,	/* = #{objects} drawn from data set (`sampsize' in R) */
	      double *dys,/* [1:(1 + (nsam * (nsam - 1))/2)]
			   * Output: to contain the distances,  -1. instead of NA */
	      int *mdata,	/*= {0,1}; 1: min(x) is missing value (NA);  0: no NA */
	      double *valmd,/*[j]= missing value code (instead of NA) for x[,j]*/
	      int *jtmd,	/* [j]= {-1,1};	 -1: x[,j] has NA; 1: no NAs in x[,j] */
	      DISS_KIND *diss_kind, // aka 'metric' \in {EUCLIDEAN, MANHATTAN, JACCARD, GOWER}
// NB: GOWER is *buggy* (<--> ../tests/clara-gower_valgrind.Rout ; ---*AND* disabled from R 2023-11-30
	      int/*logical*/ *rng_R,/*= {0,1};  0 : use clara's internal weak RNG;
				     *	        1 : use R's RNG (and seed) */
	      int/*logical*/ *pam_like,/* if (1), we do "swap()" as in pam(), otherwise
					  use the code as it was in clara() "forever"
					  upto 2011-04 */
	      int *correct_d,/* option for dist.computation: if (0), use the "fishy"
			     formula to update distances in the NA-case,
			     if (1), use a dysta2()-compatible formula */
 	      int *nrepr, /* logical (0/1): 1 = "is representative object"  */
	      int *nsel,
	      int *nbest,/* x[nbest[j],] : the j-th obs in the final sample */
	      int *nr, int *nrx,/* prov. and final "medoids" aka representatives */
	      double *radus, double *ttd, double *ratt,
	      double *ttbes, double *rdbes, double *rabes,
	      int *mtt, double *obj,
	      double *avsyl, double *ttsyl, double *sylinf,
	      int *jstop, int *trace_lev,
	      double *tmp, /* = double [ 3 * nsam ] */
	      int *itmp	/* = integer[ 6 * nsam ] */
    )
{

#define tmp1 tmp
#define tmp2 &tmp[*nsam]

#define ntmp1 itmp
#define ntmp2 &itmp[*nsam]
#define ntmp3 &itmp[nsamb]
#define ntmp4 &itmp[nsamb+ *nsam]
#define ntmp5 &itmp[2*nsamb]
#define ntmp6 &itmp[2*nsamb+ *nsam]

    /* Local variables */

    Rboolean has_NA = mdata[0];
    int j, jk, js, l, nrun;

    *jstop = 0;
    double rnn = (double) (*n), zba = -1., sx = -1.;/* Wall */

    /* n_dys := size of distance array dys[] (-1): dys[0] remains unused */
    int n_dys = *nsam * (*nsam - 1) / 2;/* >= 0 */
    int nsamb = *nsam * 2, n_sam;
    Rboolean
	full_sample = (*n == *nsam), /* only one sub sample == full data */
	lrg_sam     = (*n < nsamb);  /* sample more than *n/2 */
    if (lrg_sam)/* generate indices for the other, smaller half */
	n_sam = *n - *nsam;
    else
	n_sam = *nsam;

    if(*trace_lev) Rprintf("C clara(): (nsam,nran,n) = (%d,%d,%d);%s\n",
			   *nsam, *nran, *n,
			   full_sample ? " 'full_sample',":
			   (lrg_sam ? " 'large_sample',": ""));
    if(*rng_R && !full_sample)
	GetRNGstate();
    else /* << initialize `random seed' of the very simple randm() below */
	nrun = 0;

#define NEW_rand_k_trace_print(_nr_)					\
        rand_k= 1+ (int)(rnn* ((*rng_R)? unif_rand(): randm(&nrun)));	\
	if (rand_k > *n) {/* should never happen */			\
	    warning(_("C level clara(): random k=%d > n **\n"), rand_k); \
	    rand_k = *n;						\
	}								\
	if(*trace_lev >= 4) {						\
	    Rprintf("... {" #_nr_ "}");					\
	    if(*rng_R) Rprintf("R unif_rand()");			\
	    else       Rprintf("nrun=%5d", nrun);			\
	    Rprintf(" -> k{ran}=%d\n", rand_k);				\
	}

/* __LOOP__ :  random subsamples are drawn and partitioned into kk clusters */

    Rboolean kall = FALSE; /* kall becomes TRUE iff we've found a "valid sample",
		     i.e. one for which all d(j,k) can be computed */
    Rboolean dyst_toomany_NA = FALSE;
    int nunfs = 0;
    for (int jran = 1; jran <= *nran; ++jran) {
	if (!full_sample) {/* `real' case: sample size < n */
	    if(*trace_lev) Rprintf("- clara sample %2d ", jran);
	    int rand_k, ntt = 0;
	    if (kall && nunfs+1 != jran && !lrg_sam) {
		/* Have had (at least) one valid sample; use its representatives
		 * nrx[] :  nsel[] := sort(nrx[])  for the first j=1:k */
		if(*trace_lev >= 2) Rprintf(" if (kall && nunfs...): \n");

		for (jk = 0; jk < *kk; ++jk)
		    nsel[jk] = nrx[jk];
		for (jk = 0; jk < *kk-1; ++jk) { /* sort(nsel[0:(kk-1)] */
		    /* FIXME: nsel[] is 0-indexed, but *contains* 1-indices*/
		    int nsm = nsel[jk],
			jsm = jk;
		    for (int jkk = jk + 1; jkk < *kk; ++jkk) {
			if (nsm > nsel[jkk]) {
			    nsm = nsel[jkk];
			    jsm = jkk;
			}
		    }
		    nsel[jsm] = nsel[jk]; nsel[jk] = nsm;
		}
		ntt = *kk;
	    }
	    else { /* no valid sample  _OR_  lrg_sam */
		if(*trace_lev >= 2) Rprintf(" finding 1st... new k{ran}:\n");

		/* Loop finding random index `rand_k' not yet in nrx[0:(*kk-1)] : */
	    L180:
		NEW_rand_k_trace_print(180)

		if (kall) {
		    for (jk = 0; jk < *kk; ++jk)
			if (rand_k == nrx[jk])
			    goto L180;
		}
		/* end Loop */

		nsel[ntt] = rand_k;
		if (++ntt == n_sam)
		    goto L295;
	    }

	    if(*trace_lev >= 2) {
		Rprintf(".. kall: %s, ", (kall) ? "T" : "FALSE");
		if(*trace_lev == 2) {
		    Rprintf("nsel[ntt=%d] = %d\n", ntt, nsel[ntt]);
		} else { /* trace_lev >= 3 */
		    Rprintf("\n... nrx [0:%d]= ",*kk-1);
		    for (jk = 0; jk < *kk; jk++) Rprintf("%d ",nrx[jk]);
		    Rprintf("\n... nsel[0:%d]= ",ntt-1);
		    for (jk = 0; jk < ntt; jk++) Rprintf("%d ",nsel[jk]);
		    Rprintf("\n");
		}
	    }

	    do { /*------------------------------------------------------------
		 * Loop finding random index 'rand_k' in {1:n},
		 * not in nrx[0:(k-1)] nor nsel[1:ntt] : */
	    L210:
		NEW_rand_k_trace_print(210)

		if (kall && lrg_sam) {
		    for (jk = 0; jk < *kk; ++jk) {
			if (rand_k == nrx[jk])
			    goto L210;
		    }
		}
		/* insert rand_k into nsel[1:ntt] or after  and increase ntt : */
		for (int ka = 0; ka < ntt; ++ka)
		    if (nsel[ka] >= rand_k) {
			if (nsel[ka] == rand_k)
			    goto L210;
			else {// nsel[ka] > rand_k :
			    for (int na = ntt; na > ka; --na)
				nsel[na] = nsel[na-1];
			    nsel[ka] = rand_k;
			    /* continue _outer_ loop */ goto L290;
			}
		    }
		// else: rand_k > nsel[ka]  for all ka = 0:(ntt-1) :
		nsel[ntt] = rand_k;

	    L290:
		++ntt;
	    } while (ntt < n_sam);

	L295:
	    if(*trace_lev) Rprintf(" [ntt=%d, nunfs=%d] ", ntt, nunfs);
	    if (lrg_sam) {
		/* have indices for smaller _nonsampled_ half; revert this: */
		for (j = 1, jk = 0, js = 0; j <= *n; j++) {
		    if (jk < n_sam && nsel[jk] == j)
			++jk;
		    else
			nrepr[js++] = j;
		}
		for (j = 0; j < *nsam; ++j)
		    nsel[j] = nrepr[j];
	    }
	    if(*trace_lev >= 3) {
		Rprintf(".. nsel[1:%d]= ", *nsam);
		for (jk = 0; jk < *nsam; jk++) Rprintf("%d ",nsel[jk]);
	    }
	}
	else { /* full_sample : *n = *nsam -- one sample is enough ! */
	    for (j = 0; j < *nsam; ++j)
		nsel[j] = j+1;/* <- uses 1-indices for its *values*! */
	}
	char nl[] = " ";
	if(*trace_lev) {
	    Rprintf(" -> dysta2()");
	    if(*trace_lev >= 2) nl[0] = '\n';
	}

	dysta2(*nsam, *jpp, nsel, x, *n, dys, *diss_kind,
	       jtmd, valmd, has_NA, &dyst_toomany_NA);
	if(dyst_toomany_NA) {
	    if(*trace_lev)
		Rprintf(" gave dyst_toomany_NA --> new sample.%s", nl);
	    dyst_toomany_NA = FALSE;
	    ++nunfs;
	    continue;/* random sample*/
	} else if(*trace_lev) Rprintf(";%s", nl);

	double s = 0., sky;
	for(l = 1; l <= n_dys; l++) /* dys[0] is not used here */
	    if (s < dys[l])
		s = dys[l];
	if(*trace_lev >= 2)
	    Rprintf(" clara -> s:= max{dys[1..%d]} = %g;", n_dys,s);

	bswap2(*kk, *nsam, s, dys, *pam_like, *trace_lev,
	       /* --> */ &sky, nrepr,
	       /* dysma */tmp1, /*dysmb*/tmp2,
	       /* beter[], only used here */&tmp[nsamb]);

	if(*trace_lev >= 2)
	    Rprintf(" end{bswap2}: sky = %g\n", sky);

	Rboolean nafs;
	double zb;
	selec(*kk, *n, *jpp, *diss_kind, &zb, *nsam, has_NA, jtmd, valmd,
	      *trace_lev, nrepr, nsel, dys, x, nr, &nafs, ttd, radus, ratt,
	      ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, tmp1, tmp2, *correct_d);

	if (nafs) { /* couldn't assign some observation (to a cluster)
		     * because of too many NA s */
	    ++nunfs;
	    if(*trace_lev >= 2)
		Rprintf(" selec() -> 'NAfs'");
	}
	else if(!kall || zba > zb) { /* 1st proper sample  or  new best */
	    if(*trace_lev >= 2) Rprintf(kall ? "new best" : "1st proper sample");
	    kall = TRUE;
	    zba = zb;
	    for (jk = 0; jk < *kk; ++jk) {
		ttbes[jk] = ttd	 [jk];
		rdbes[jk] = radus[jk];
		rabes[jk] = ratt [jk];
		nrx  [jk] = nr	 [jk];
	    }
	    for (js = 0; js < *nsam; ++js)
		nbest[js] = nsel[js];
	    sx = s;
	}
	/* if(*trace_lev) Rprintf("%sobj= %g\n", (*trace_lev < 2) ?", " :" ",  zb/rnn); */
	if(*trace_lev) Rprintf(" obj= %g\n", zb/rnn);

	if(full_sample) break; /* out of resampling */
    }
/* --- end random sampling loop */
    if(*rng_R && !full_sample)
	PutRNGstate();

    if(*trace_lev && !full_sample) {
	Rprintf("C clara() -> best sample _found_ ");
	if(*trace_lev >= 2) {
	    Rprintf("; nbest[1:%d] =\n c(", *nsam);
	    for (js = 0; js < *nsam; ) {
		Rprintf("%d", nbest[js]);
		if(++js < *nsam) Rprintf((js % 20) ? "," : ",\n  ");
	    }
	    Rprintf(")\n");
	}
    }

    *obj = zba / rnn;

    if (nunfs >= *nran) { *jstop = 1; return; }
    /* else */
    if (!kall) { *jstop = 2; return; }

/* for the best subsample, the objects of the entire data set
   are assigned to their clusters */
    dysta2(*nsam, *jpp, nbest, x, *n, dys, *diss_kind, jtmd, valmd,
	   has_NA, &dyst_toomany_NA);
    if(dyst_toomany_NA) {
	error(_("clara()'s C level dysta2(nsam=%d, p=%d, n=%d) gave 'toomany_NA'"),
	      *nsam, *jpp, *n);
    }
    if(*trace_lev) Rprintf(" resul(), ");
    resul(*kk, *n, *jpp, *diss_kind, has_NA, jtmd, valmd, x, nrx, mtt, *correct_d);

    if (*kk > 1) {
	if(*trace_lev) Rprintf(" black()");
	black(*kk, *jpp, *nsam, nbest, dys, sx, x,
	      /* compute --> */
	      avsyl, ttsyl, sylinf,
	      ntmp1, ntmp2, ntmp3, ntmp4, /* syl[] */ tmp1, tmp2);
    }
    if(*trace_lev) Rprintf(" and return() from C.\n");
    return;
} /* End clara() ---------------------------------------------------*/
#undef tmp1
#undef tmp2

#undef ntmp1
#undef ntmp2
#undef ntmp3
#undef ntmp4
#undef ntmp5
#undef ntmp6



/**
 * Compute Dissimilarities for the selected sub-sample ---> dys[,]
 */
void dysta2(int nsam, int jpp, int *nsel,
	    double *x, int n, double *dys, DISS_KIND diss_kind,
	    int *jtmd, double *valmd, Rboolean has_NA, Rboolean *toomany_NA)
{
    int nlk = 0;
    int current_group = 0; // for GOWER only: Counter for overall index
    dys[0] = 0.;/* very first index; *is* used because ind_2(i,i) |-> 0 ! */
    for (int l = 1; l < nsam; ++l) {
	int lsel = nsel[l];
	if(lsel <= 0 || lsel > n)
	    error(_("C level dysta2(): nsel[%s= %d] = %d is outside 0..n, n=%d"),
		  "l", l, lsel, n);
	for (int k = 0; k < l; ++k) { /* compute d(nsel[l], nsel[k]) {if possible}*/
	    int ksel = nsel[k];
	    if(ksel <= 0 || ksel > n)
		error(_("C level dysta2(): nsel[%s= %d] = %d is outside 0..n, n=%d"),
		      "k", k, ksel, n);
	    ++nlk;
	    int npres = 0, j, lj, kj, N_ones = 0;
	    double clk = 0.;
	    for (j = 0, lj = lsel-1, kj = ksel-1; j < jpp;
		 ++j, lj += n, kj += n) {
		if (has_NA && jtmd[j] < 0) { /* x[,j] has some Missing (NA) */
		    /* in the following line (Fortran!), x[-2] ==> seg.fault
		       {BDR to R-core, Sat, 3 Aug 2002} */
		    if (x[lj] == valmd[j] || x[kj] == valmd[j]) {
			continue /* next j */;
		    }
		}
		++npres; // number of present pairs
		if (diss_kind == EUCLIDEAN)
		    clk += (x[lj] - x[kj]) * (x[lj] - x[kj]);
		else if (diss_kind == JACCARD) {
		    if( x[lj] > 0.9 && x[kj] > 0.9) { // both "are 1" - increment numerator
			clk++ ; N_ones++ ;
		    } else if( x[lj] > 0.9 || x[kj] > 0.9)// any is 1 - increment N_ones
			N_ones++ ;
		}
                else if (diss_kind == MANHATTAN)
                    clk += fabs(x[lj] - x[kj]);
// NB: GOWER is *buggy* (<--> ../tests/clara-gower_valgrind.Rout ; ---*AND* disabled from R 2023-11-30
                else if (diss_kind == GOWER) {
                    if (x[lj] == x[kj]) {
                        continue /* next j */;
                    }
                    double diff;
                    if (jtmd[j] >= 0) { // Numerical feature

                        diff = fabs(x[lj] - x[kj]);

                        // Find the maximum and minimum values of this feature
                        double max_val = -1.0 / 0.0;
                        double min_val = 1.0 / 0.0;
                        for (int i = 0; i < nsam; ++i){

                            double val = x[i + nsam * current_group];
                            if (val > max_val) {
                                max_val = val;
                            }
                            if (val < min_val) {
                                min_val = val;
                            }

                        }
                        current_group = (current_group + 1) % 2; // Toggle between 0 and 1
                        double range = max_val - min_val;
                        diff /= range;

                    } else { // Categorical feature
                        diff = (x[lj] != x[kj]);
                    }
                    clk += diff;
                }
	    } /* for( j ..) */
	    if (npres == 0) {/* cannot compute d(.,.) because of too many NA */
		*toomany_NA = TRUE;
		dys[nlk] = -1.;
	    } else {
		double d1 = clk * (jpp / (double) npres);
		dys[nlk] =
		      (diss_kind == EUCLIDEAN) ? sqrt(d1)
		    : (diss_kind ==  JACCARD)  ? 1 - clk / (double) N_ones
                    : (diss_kind == MANHATTAN) ? d1
                    : /* diss_kind == GOWER */ d1 / jpp;
	    }
	} /* for( k ) */
    } /* for( l ) */
    return;
} /* End dysta2() -----------------------------------------------------------*/

double randm(int *nrun)
{
/* we programmed this generator ourselves because we wanted it
   to be machine independent. it should run on most computers
   because the largest int used is less than 2^30 . the period
   is 2^16=65536, which is good enough for our purposes. */
    /* MM: improved the original speed-wise only: */
    *nrun = (*nrun * 5761 + 999) & 0177777;
    /* Masking off all but the last 16 bits is equivalent to  % 65536 */
    return ((double) (*nrun) / 65536.);
} /* randm() */

/* bswap2() : called once [per random sample] from clara() : */
void bswap2(int kk, int n, /* == nsam == 'sampsize', here in clara */
	    double s, const double dys[],
	    Rboolean pam_like, int trace_lev,
	    // result (though *only* nrepr[] is really used in caller:)
	    double *sky, int *nrepr,
	    double *dysma, double *dysmb, double *beter)
{
    int i, j, ij, k,h, hbest = -1, nbest = -1;/* init for -Wall */
    double dzsky;

    /* Parameter adjustments */
    --nrepr;
    --beter;

    --dysma;	--dysmb;

    if(trace_lev >= 2) {
	if(trace_lev == 2)
	    Rprintf("\n bswap2(), 1. BUILD:");
	else
	    Rprintf("\n bswap2(*, s=%g), 1. BUILD: ", s);
    }

    s = s * 1.1 + 1.;/* value larger than all dissimilarities */

/* ====== first algorithm: BUILD. ====== */

    for (i = 1; i <= n; ++i) {
	nrepr[i] = 0;
	dysma[i] = s;
    }

    for(k = 0; k < kk; k++) {
	int nmax = -1; /* -Wall */
	double ammax = 0.;
	for (i = 1; i <= n; ++i) {
	    if (nrepr[i] == 0) {
		beter[i] = 0.;
		for (j = 1; j <= n; ++j) {
		    double cmd = dysma[j] - dys[ ind_2(i, j)];
		    if (cmd > 0.)
			beter[i] += cmd;
		}
		if (ammax <= beter[i]) {
		    /*	does < (instead of <= ) work too? -- NO! */
		    ammax = beter[i];
		    nmax = i;
		}
	    }
	}

	nrepr[nmax] = 1;/* = .true. : found new representative */
	if(trace_lev >= 2) {
	    if(trace_lev == 2)
		Rprintf(" %d", nmax);
	    else
		Rprintf("    new repr. %d\n", nmax);
	}

	/* update dysma[] : dysma[j] = D(j, nearest_representative) */
	for (j = 1; j <= n; ++j) {
	    ij = ind_2(nmax, j);
	    if (dysma[j] > dys[ij])
		dysma[j] = dys[ij];
	}
    }
    // output of the above loop:  nrepr[], dysma[], ...

    *sky = 0.;
    for (j = 1; j <= n; ++j)
	*sky += dysma[j];

    if(trace_lev >= 2) /* >= 2 (?) */ {
	Rprintf("  after build: medoids are");
	for (i = 1; i <= n; ++i)
	    if(nrepr[i] == 1) Rprintf(" %2d", i);
	if(trace_lev >= 3) {
	    Rprintf("\n  and min.dist dysma[1:n] are\n");
	    for (i = 1; i <= n; ++i) {
		Rprintf(" %6.3g", dysma[i]);
		if(i % 10 == 0) Rprintf("\n");
	    }
	    if(n % 10 != 0) Rprintf("\n");
	} else Rprintf("\n");
	Rprintf(" --> sky = sum_j D_j= %g\n", *sky);
    }

    if (kk == 1)
	return;

// asky = *sky / ((double) n);

/* ====== second algorithm: SWAP. ====== */

/* Big LOOP : */
L60:

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

    dzsky = 1.; /* 1 is arbitrary > 0; only dzsky < 0 matters in the end */
    for (h = 1; h <= n; ++h) if (!nrepr[h]) {
        for (i = 1; i <= n; ++i) if (nrepr[i]) {
	    double dz = 0.;
	    /* dz := T_{ih} := sum_j C_{jih}  [p.104] : */
	    for (j = 1; j <= n; ++j) {
		int ij = ind_2(i, j),
		    hj = ind_2(h, j);
		if (dys[ij] == dysma[j]) {
		    double small;
		    if(pam_like)
			small = dysmb[j] > dys[hj] ? dys[hj] : dysmb[j];
		    else // old clara code which differs from pam()'s
			// and seems a bit illogical:
			small = dysmb[j] > dys[ij] ? dys[hj] : dysmb[j];
		    dz += (- dysma[j] + small);
		}
		else if (dys[hj] < dysma[j])
		    dz += (- dysma[j] + dys[hj]);
	    }
	    if (dzsky > dz) {
		dzsky = dz; // dzsky := min_{i,h} T_{i,h}
		hbest = h;
		nbest = i;
	    }
	}
    }

    /* once had some 64-bit compiler / data configuration that looped forever*/
    R_CheckUserInterrupt();

    if (dzsky < 0.) { /* found an improving swap */
	if(trace_lev >= 3)
	    Rprintf( "   swp new %d <-> %d old; decreasing diss. by %g\n",
		     hbest, nbest, dzsky);
	nrepr[hbest] = 1;
	nrepr[nbest] = 0;
	*sky += dzsky;
	goto L60;
    }
    if(trace_lev >= 2 && hbest != -1) // in my examples  hbest == -1 and it does not print:
	Rprintf( "  Last swap: new %d <-> %d old; decreasing diss. by %g\n",
		 hbest, nbest, dzsky);

} /* End of bswap2() -------------------------------------------------- */

/* selec() : called once [per random sample] from clara() */
void selec(int kk, int n, int jpp, DISS_KIND diss_kind,
	   double *zb, int nsam, Rboolean has_NA, int *jtmd, double *valmd,
	   int trace_lev, // currently unused -- FIXME !?
	   int *nrepr, int *nsel, double *dys, double *x, int *nr,
	   Rboolean *nafs, /* := TRUE if a distance cannot be calculated */
	   double *ttd, double *radus, double *ratt,
	   // [i]tmp* for clara(), i.e. not used later!
	   int *nrnew, int *nsnew, int *npnew, int *ns, int *np, int *new,
	   double *ttnew, double *rdnew, int correct_d)
{
/* Parameter adjustments */
    --nsel;    --nrepr;

    --ratt;
    --radus; --ttd;    --np;	--nr;	 --ns;

    --rdnew; --ttnew; --npnew; --nrnew; --nsnew;
    --new;

    /* nafs := TRUE  if a distance cannot be calculated (because of NA s)*/
    *nafs = FALSE;

    /* identification of representative objects, and initializations */
    for (int j = 1, jk = 0; j <= nsam; ++j) {
	if (nrepr[j] != 0) {
	    ++jk;
	    nr	 [jk] = nsel[j];
	    ns	 [jk] = 0;
	    ttd	 [jk] = 0.;
	    radus[jk] = -1.;
	    np	 [jk] = j;
	}
    }

/* - assignment of the objects of the entire data set to a cluster,
 * - computation of some statistics,
 * - determination of the new ordering of the clusters */

    double pp = (double) (jpp);
    *zb = 0.;
    int newf = 0;
    for(int i = 1; i <= n; i++) { /* for  obj  i,  find its cluster  jk in 1:kk */
	int jkabc = -1/* -Wall */;
	double dsum,
	    dnull = -9./* -Wall */;
	if (!has_NA) {
	    for (int jk = 1; jk <= kk; ++jk) {
		dsum = 0.;
		int nr_k = nr[jk];
		if (nr_k != i) {
		    int N_ones = 0;
		    double tra = 0.; // init only for JACCARD
		    for (int jp = 0; jp < jpp; ++jp) {
			int na = (nr_k - 1) + jp * n,
			    nb = (i    - 1) + jp * n;
			// d_j(x[na], x[nb]) :
			if (diss_kind == JACCARD) {
			    if(x[na] > 0.9 && x[nb] > 0.9) {
				// both "are 1" - increment numerator (and denom.)
				tra += 1; N_ones ++;
			    } else if( x[na] > 0.9 || x[nb] > 0.9) {
				// any is 1 - increment denominator N_ones
				N_ones ++;
			    }
			} else { // Euclidean or Manhattan
			    tra = fabs(x[na] - x[nb]);
			    if (diss_kind == EUCLIDEAN)
				tra *= tra;
			    dsum += tra;
			}
		    }
		    if (diss_kind == JACCARD)
			dsum = 1 - tra / (double)N_ones;

		    if (jk != 1 && dsum >= dnull)
			continue /* next jk */;

		} /* (nr_k != i) */

		// new best: dsum < "previous" dnull
		//---------------- FIXME for JACCARD -------------------------
		dnull = dsum;
		jkabc = jk;
	    } /* for (jk  1:kk) */
	}
	else { // _has_ missing data
	    Rboolean first = TRUE;
	    for (int jk = 1; jk <= kk; ++jk) {
		dsum = 0.;
		int nr_k = nr[jk];
		if (nr_k != i) {
		    int nobs = 0, N_ones = 0;
		    double tra = 0.; // init only for JACCARD
		    for (int jp = 0; jp < jpp; ++jp) {
			int na = (nr_k - 1) + jp * n,
			    nb = (i    - 1) + jp * n;
			if (jtmd[jp] < 0) {
			    if (x[na] == valmd[jp] || x[nb] == valmd[jp])
				continue /* next jp */;
			}
			nobs++;
			if (diss_kind == JACCARD) {
			    if(x[na] > 0.9 && x[nb] > 0.9) {
				// both "are 1" - increment numerator (and denom.)
				tra += 1; N_ones ++;
			    } else if( x[na] > 0.9 || x[nb] > 0.9) {
				// any is 1 - increment denominator N_ones
				N_ones ++;
			    }
			} else { // Euclidean or Manhattan
			    tra = fabs(x[na] - x[nb]);
			    if (diss_kind == EUCLIDEAN)
				tra *= tra;
			    dsum += tra;
			}
		    }
		    if (nobs == 0) /* all pairs partially missing */
			continue /* next jk */;
		    if (diss_kind == JACCARD)
			dsum = 1 - tra / (double)N_ones;
		    if(correct_d) // correct -- only since 2017-06
			dsum *= (pp / nobs);
		    else
			dsum *= (nobs / pp);
		}
		if (first)
		    first = FALSE;
		else if (dnull <= dsum)
		    continue /* next jk */;
		/* here : first was TRUE {i.e. 1st time} or
		 *	  dnull > dsum	 {i.e. new best} */
		dnull = dsum;
		jkabc = jk;
	    }/* for(jk ..) */

	    if (first) { /* found nothing */
		*nafs = TRUE; return;
	    }
	} /* else: has_NA */

	if (diss_kind == EUCLIDEAN)
	    dnull = sqrt(dnull);

	*zb += dnull;
	ttd[jkabc] += dnull;
	if (radus[jkabc] < dnull)
	    radus[jkabc] = dnull;

	++ns[jkabc];
	if (newf < kk) {
	    if (newf != 0) {
		for (int jnew = 1; jnew <= newf; ++jnew) {
		    if (jkabc == new[jnew])
			goto L90;/* next i */
		}
	    }
	    ++newf;
	    new[newf] = jkabc;
	}
    L90:
	;
    } /* for( i = 1..n ) -------------------------------------------------*/


/*     a permutation is carried out on vectors nr,ns,np,ttd,radus
     using the information in vector new. */

    for (int j = 1; j <= kk; ++j) {
	int nj = new[j];
	nrnew[j] = nr[nj];
	nsnew[j] = ns[nj];
	npnew[j] = np[nj];
	ttnew[j] = ttd[nj];
	rdnew[j] = radus[nj];
    }
    for (int j = 1; j <= kk; ++j) {
	nr[j] = nrnew[j];
	ns[j] = nsnew[j];
	np[j] = npnew[j];
	ttd[j] = ttnew[j] / (double) ns[j];
	radus[j] = rdnew[j];
    }

    if (kk > 1) {
	/* computation of ratt[ka] := minimal distance of medoid ka to any
	   other medoid for comparison with the radius of cluster ka. */

	for (int ka = 1; ka <= kk; ++ka) {
	    Rboolean first = TRUE;
	    int npa = np[ka];
	    for (int kb = 1; kb <= kk; ++kb) {
		if (kb == ka)
		    continue /* next kb */;

		int npb = np[kb],
		    npab = ind_2(npa, npb);
		if (first)
		    first = FALSE;
		else if (dys[npab] >= ratt[ka])
		    continue /* next kb */;

		ratt[ka] = dys[npab];
		if (ratt[ka] == 0.)
		    ratt[ka] = -1.;
	    }
	    if (ratt[ka] > -0.5)
		ratt[ka] = radus[ka] / ratt[ka];
	}
    }
    return;
} /* End selec() -----------------------------------------------------------*/

void resul(int kk, int n, int jpp, DISS_KIND diss_kind, Rboolean has_NA,
	   int *jtmd, double *valmd, double *x, int *nrx, int *mtt,
	   int correct_d) /* option for dist.computation:
              if (0), use the "fishy" formula to update distances in the NA-case,
 	      if (1), use a dysta2()-compatible formula */
{
// __FIXME__  "Jaccard" not yet supported ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<


/* clustering vector is incorporated into x, and ``printed''. */

    for(int i = 0; i < n; i++) {

	for (int jk = 0; jk < kk; ++jk) {
	    if (nrx[jk] == i + 1)/* 1-indexing */
		goto L220; /* continue next i (i.e., outer loop) */
	}
	int njnb = i, nrjk, /* = nrx[jk -1 */
	    jksky = -1/* Wall */;
	double dnull = -9./* Wall */;
	if (!has_NA) {
	    for (int jk = 0; jk < kk; ++jk) {
		double dsum = 0.;
		nrjk = (nrx[jk] - 1);
		for (int j = 0; j < jpp; ++j) {
		    double tra = fabs(x[nrjk + j * n] -
				      x[njnb + j * n]);
		    if (diss_kind == EUCLIDEAN)
			tra *= tra;
		    dsum += tra;
		}
		if (diss_kind == EUCLIDEAN)
		    dsum = sqrt(dsum);
		if (jk == 0 || dnull > dsum) { // have new best
		    dnull = dsum;
		    jksky = jk;
		}
	    }
	}
	else { /* _has_ missing data */
	    double pp = (double) jpp;
	    for (int jk = 0; jk < kk; ++jk) {
		double dsum = 0.;
		nrjk = (nrx[jk] - 1);
		int nobs = 0;
		for (int j = 0; j < jpp; ++j) {
		    int na = nrjk + j * n,
			nb = njnb + j * n;
		    if (jtmd[j] < 0) {
			if (x[na] == valmd[j] || x[nb] == valmd[j])
			    continue /* next j */;
		    }
		    nobs++;
		    double tra = fabs(x[na] - x[nb]);
		    if (diss_kind == EUCLIDEAN)
			tra *= tra;
		    dsum += tra;
		}
		if (diss_kind == EUCLIDEAN)
		    dsum = sqrt(dsum);
		if(correct_d) // correct -- only since 2016-04
		    dsum *= (pp / nobs);
		else
		    dsum *= (nobs / pp); // MM: "fishy" (had note since r4321, 2007-05-01 !)

		if (jk == 0 || dnull > dsum) { // have new best
		    dnull = dsum;
		    jksky = jk;
		}
	    }
	}
	x[njnb] = (double) jksky + 1;/* 1-indexing */

    L220:
	;
    } /* for(i = 0; i < n ..)*/

    for (int jk = 0; jk < kk; ++jk)
	x[nrx[jk] - 1] = (double) jk + 1;/* 1-indexing */

    /* mtt[k] := size(k-th cluster) : */
    for (int ka = 0; ka < kk; ++ka) {
	mtt[ka] = 0;
	for(int i = 0; i < n; i++) {
	    if (((int) x[i]) == ka + 1)/* 1-indexing */
		++mtt[ka];
	}
    }
    return;
} /* end resul() -----------------------------------------------------------*/


// called 'dark()' in  ./pam.c
void black(int kk, int jpp, int nsam, int *nbest,
	   double *dys, double s, double *x,
	   /* --> Output : */
	   double *avsyl, double *ttsyl, double *sylinf,
	   /* but the following output vectors are never used by clara() : */
	   int *ncluv, int *nsend, int *nelem, int *negbr,
	   double *syl, double *srank)
{
/* Silhouettes computation and "drawing"  --> syl[] and sylinf[] */

    /* System generated locals */
    int sylinf_dim1 = nsam,
	sylinf_offset = 1 + sylinf_dim1 * 1;


/* Parameter adjustments */
    --avsyl;

    --srank;    --syl;
    --negbr; --nelem; --nsend;
    --ncluv;	--nbest;

    sylinf -= sylinf_offset;

/*
     construction of clustering vector (ncluv)
     of selected sample (nbest).
*/

    for (int l = 1; l <= nsam; ++l) {
	ncluv[l] = (int) x[nbest[l] - 1];
    }

/*     "drawing" of the silhouettes */

    int nsylr = 0;
    *ttsyl = 0.;
    for (int numcl = 1; numcl <= kk; ++numcl) {
	int ntt = 0, j, l;
	double dysa, dysb;
	for (j = 1; j <= nsam; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) {
	    int nj = nelem[j];
	    dysb = s * 1.1 + 1.;
	    negbr[j] = -1;

	    for (int nclu = 1; nclu <= kk; ++nclu) {
		if (nclu != numcl) {
		    int nbb = 0;
		    double db = 0.;
		    for (l = 1; l <= nsam; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    db += dys[ind_2(nj, l)];
			}
		    }
		    db /= (double) nbb;
		    if (db < dysb) {
			dysb = db;
			negbr[j] = nclu;
		    }
		}
	    }

	    if (ntt == 1) {
		syl[j] = 0.;	    continue /* j */;
	    }
	    dysa = 0.;
	    for (l = 1; l <= ntt; ++l) {
		int nl = nelem[l];
		dysa += dys[ind_2(nj, nl)];
	    }
	    dysa /= (double) (ntt - 1);
	    if (dysa <= 0.) {
		if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;

		continue /* j */;
	    }

	    if (dysb > 0.) {
		if (dysb > dysa)
		    syl[j] = 1. - dysa / dysb;
		else if (dysb < dysa)
		    syl[j] = dysb / dysa - 1.;
		else /* (dysb == dysa) */
		    syl[j] = 0.;

		if (syl[j] < -1.)
		    syl[j] = -1.;
		else if (syl[j] > 1.)
		    syl[j] = 1.;
	    }
	    else {
		syl[j] = -1.;
	    }

	} /* for(j ..) */

	avsyl[numcl] = 0.;
	for (j = 1; j <= ntt; ++j) {
	    int lang = -1/* -Wall */;
	    double symax = -2.;
	    for (l = 1; l <= ntt; ++l) {
		if (syl[l] > symax) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[numcl];
	avsyl[numcl] /= ntt;

	if (ntt >= 2) {
	    for (l = 1; l <= ntt; ++l) {
		int lplac = nsend[l];
		++nsylr;
		sylinf[nsylr + sylinf_dim1       ] = (double) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac];
		sylinf[nsylr + sylinf_dim1 * 3   ] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[nelem[lplac]];
	    }
	}
	else {
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1       ] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3   ] = 0.;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[nelem[1]];
	}

    } /*  for( numcl ) */
    *ttsyl /= (double) (nsam);
    return;
} /* black */
