/* Declare everything, Fortran & C -- so we can register them */

#include <R.h>
#include <Rinternals.h>
/* -> Rconfig.h, but also Boolean.h RS.h */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("cluster", String)
#else
#define _(String) (String)
#endif

// These codes must match those in ../R/clara.q  <==>  'diss_kind'
typedef enum {
    EUCLIDEAN = 1,
    MANHATTAN = 2,
    JACCARD = 3
  , ZJACCARD  = 4
  , SMATCHING = 5
  , GOWER = 6
} DISS_KIND;

/* --------- ./clara.c ------------------*/

double randm(int *nrun);

void cl_clara(int *n,  /* = number of objects */
	      int *jpp,/* = number of variables */
	      int *kk, /* = number of clusters, 1 <= kk <= n-1 */
	      double *x,	/* Input:  the data x[n, jpp] _rowwise_ (transposed)
				 * Output: the first `n' values are the `clustering'
				 *	   (integers in 1,2,..,kk) */
	      int *nran,	/* = #{random samples} drawn	   (= `samples' in R)*/
	      int *nsam,	/* = #{objects} drawn from data set (`sampsize' in R) */
	      double *dys,/* [1:(1 + (nsam * (nsam - 1))/2)]
			   * Output: to contain the distances */
	      int *mdata,	/*= {0,1}; 1: min(x) is missing value (NA);  0: no NA */
	      double *valmd,/*[j]= missing value code (instead of NA) for x[,j]*/
	      int *jtmd,	/* [j]= {-1,1};	 -1: x[,j] has NA; 1: no NAs in x[,j] */
	      DISS_KIND *diss_kind, // = {EUCLIDEAN, MANHATTAN, ....}
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
    );

void bswap2(int kk, int nsam, double s, const double dys[],
	    Rboolean pam_like, int trace_lev,
	    // result:
	    double *sky, int *nrepr,
	    double *dysma, double *dysmb, double *beter);

void selec(int kk, int n, int jpp, DISS_KIND diss_kind,
	   double *zb, int nsam, Rboolean has_NA, int *jtmd, double *valmd,
	   int trace_lev,
	   int *nrepr, int *nsel, double *dys, double *x, int *nr,
	   Rboolean *nafs, double *ttd, double *radus, double *ratt,
	   int *nrnew, int *nsnew, int *npnew, int *ns, int *np, int *new,
	   double *ttnew, double *rdnew, int correct_d);

void resul(int kk, int n, int jpp, DISS_KIND diss_kind, Rboolean has_NA,
	   int *jtmd, double *valmd, double *x, int *nrx, int *mtt, int correct_d);

void black(int kk, int jpp, int nsam, int *nbest,
	   double *dys, double s, double *x,
	   /* --> Output : */
	   double *avsyl, double *ttsyl, double *sylinf,
	   int *ncluv, int *nsend, int *nelem, int *negbr,
	   double *syl, double *srank);

/* -------- ./dysta.c --- (dysta.f was in pam.f) -------------------- */
int dysta(int *nn, int *jpp, double *x, double *dys, int *ndyst,
	  int *jtmd, double *valmd);
/* --------- ./pam.c ------------------*/

#ifdef _UNUSED_C_pam
void cl_pam(int *nn, int *jpp, int *kk, double *x, double *dys,
	    int *jdyss, /* jdyss = 0 : compute distances from x
			 *	      = 1 : distances provided	in x */
	    double *valmd, int *jtmd,
	    int *ndyst, int *nsend, int *nrepr, int *nelem,
	    double *radus, double *damer, double *ttd, double *separ,
	    double *ttsyl, double *obj, int *med, int *ncluv,
	    double *clusinf, double *sylinf, int *nisol, int* optim);
#endif

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
	    SEXP dist_kind); // = 1 ("euclidean")  or 2 ("manhattan")

void bswap(int kk, int nsam, int *nrepr,
	   /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
	   Rboolean med_given, Rboolean do_swap, int trace_lev,
	   double *dysma, double *dysmb, double *beter,
	   const double dys[], double s, double *obj, int pamonce);

void cstat(int kk, int nn, int *nsend, int *nrepr, Rboolean all_stats,
	   double *radus, double *damer, double *ttd, double *separ, double *s,
	   double *dys, int *ncluv, int *nelem, int *med, int *nisol);

void dark(int kk, int nn, const int ncluv[], const double dys[], double s,
	  int *nsend, int *nelem, int *negbr,
	  double *syl, double *srank, double *avsyl, double *ttsyl,
	  double *sylinf);


/* --------- ./spannel.c ------------------*/

void cl_sweep(double *, int *, int *, int *, double *);

void spannel(int *ncas, /* = number of objects */
	     int *ndep, /* = number of variables */
	     double *dat,/* [ncas, 0:ndep] */
	     double *dstopt, /* = squared distances [1:ncas] */
	     double *cov,/* matrix [0:ndep, 0:ndep] */
	     double *varsum,	/* [1:ndep] */
	     double *varss,	/* [1:ndep] */
	     double *prob,	/* [1:ncas] */
	     double *work,	/* [0:ndep] */
	     double *eps,
	     int *maxit, /* = maximal # iterations (and returns #{iter.})*/
	     int *ierr);

void sildist(double *d,		/* distance : in matrix or dist format; i.e.,
				   of length n^2 or n*(n-1)/2; see 'ismat' */
	     int    *n,		/* number of Subjects (attr(d,'Size')) */
	     int    *clustering,/* clustering vector, values from {1..k} */
	     int    *k,		/* number of clusters */
	     double *diC,	/* diC */
	     int    *counts,	/* counts[k] :=  #{cluster k} */
	     double *si,	/* (a_i - b_i) / max(ai,bi) */
	     int    *neighbor,	/* neighbor */
	     int    *ismat);	/* boolean : is 'd' a matrix or 'dist' ? */

void cl_fanny(int *nn, int *jpp, int *kk,
	      double *x, double *dss, int *jdyss, double *valmd,
	      int *jtmd, int *ndyst, int *nsend, int *nelem,
	      int *negbr, double *syl, double *p, double *dp,
	      double *pt, int *nfuzz, double *esp, double *ef,
	      double *dvec, double *ttsyl, double *obj,
	      int *ncluv, double *sylinf, double *r, double *tol, int *maxit);


/* ================= Fortran things (remainder) ======================== */

/* -------- ./daisy.c ---------------------------------- */
void cldaisy(int *nn, int *jpp, double *x,
		       double *valmd, double *weights,
		       int *jtmd, int *jdat, int *vtype,
		       int *ndyst, int *mdata, double *disv);

/* -------- ./mona.c ---------------------------------- */
void clmona(int *nn, int *pp, int *x, int *jerr,
		      int *nban, int *ner, int *kwan, int *lava, int *jlack);

/* -------- ./twins.c ---------------------------------- */
void R_bncoef(int *nn, double *ban, double *cf);
double bncoef(int  nn, double *ban);

void twins(int *nn, int *jpp, double *x,
	   double *dys, double *dys2, int *jdyss, double *valmd,
	   int *jtmd, int *ndyst, int *jalg, int *method,
	   int *kwan, int *ner, double *ban, double *coef,
	   double *alpha, int *merge, int *trace_lev);
