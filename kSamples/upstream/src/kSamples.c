/* helper functions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R_ext/Rdynload.h>


/* for random number generator in R */
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "myfuns.h"

/* Code based to a large extent on Angie Zhu's original 
   version, last modified 06/11/2015 (FWS)*/

/* dynamically allocates memory for a double array of length n
      and returns the pointer;
   prints out error message if the allocation is unsuccessful */
double *dmalloc(unsigned long n) {
	double *x;
	x = (double*) malloc((size_t) n * sizeof(double));
/*	if (x == NULL) {
		printf("Error: Could not allocate %ld doubles\n", n);
	}
*/
	return(x);
}

/* counts and returns the number of occurrence of a given double 
   number z in a double array dat of length n */
int getCount(double z, double *dat, int n) {
	int i;
	int count = 0;
	
	for (i = 0; i < n; i++) {
		if (dat[i] == z) {
			count++;
		}
	}
	
	return(count);
}

/* computes and returns the sum of elements in a given integer 
   array x of length n */ 
int getSum(int *x, int n) {
	int i; 
	int sum = 0; 
	
	for (i = 0; i < n; i++) { 
		sum += x[i]; 
	} 
	
	return(sum);
}

/* dynamically allocates memory for an integer array with length 
   n and returns its pointer;
   prints out error message if the allocation is unsuccessful */
int *imalloc(unsigned long n) {
	int *x;
	x = (int*) malloc((size_t) n * sizeof(int));
/*	if (x == NULL) {
		printf("Error: Could not allocate %ld ints\n", n);
	}
*/
	return(x);
}

/* produces a copy of double n*m matrix x */
void mcopy(double *x, double *copy, int n, int m) {
	int i;
	n *= m;
	for (i = 0; i < n; i++) { 
		*(copy + i) = *(x + i);
	}
}

/* produces a copy of int n*m matrix X */
void imcopy(int *x, int *copy, int n, int m) {
	int i;
	n *= m;
	for (i = 0; i < n; i++) { 
		*(copy + i) = *(x + i);
	}
}

/* dynamically allocates memory for an array of pointers to 
   double arrays with length n and returns the pointer;
 prints out error message if the allocation is unsuccessful */
 
double **pdmalloc(unsigned long n) {
	double **x;
	x = (double**) malloc((size_t) n * sizeof(double*));
/*	if (x == NULL) {
		printf("Error: Could not allocate %ld pointers 
                   to double arrays\n", n);
	}
*/
	return(x);
}

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * adkTestStat
 * Last modified: 07.21.2015 (FWS)
 ***************************/

/* static variables */
static int k;
static int L;
static int nsum;
static int ncomb;
static int getSmat; 
/* logical value to indicate if smat will be generated */
static int count; 
/* keeps track of how many combinations have been visited */
static int *ns;
static int *ivec;
static int dimst;
static double *xvec;
static double *zstar;
static double *teststat;
static double *pval;
static double *smat; 
/* statistic matrix, recording the statistics*/ 


/* initializes static variables */
void initvals(int kk, double *xxvec, int *nns, double *zzstar, 
              int LL, int *iivec, double *tteststat, 
              double *ppval, int nnsum, int nncomb, 
              int ggetSmat, double *ssmat, int ddimst) {
	k = kk;
	xvec = xxvec;
	ns = nns;
	zstar = zzstar;
	L = LL;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	ncomb = nncomb;
	getSmat = ggetSmat;
	smat = ssmat;
	count = 0;
	dimst = ddimst;
}

/* initializes static variables */
void initvals1(int kk, double *xxvec, int *nns, double *zzstar, 
              int LL){ 
	k = kk;
	xvec = xxvec;
	ns = nns;
	zstar = zzstar;
	L = LL;
}

/* computes the k-sample Anderson-Darling test statistics in
   both original and alternative versions for the nonparametric
   (rank) test described in Scholz F.W. and Stephens M.A. 
   (1987), K-sample Anderson-Darling Tests, Journal of the
   American Statistical Association, Vol 82, No. 399, 
   pp. 918-924
	
	Arguments:
	adk: double array of length 2, stores AkN2 and AakN2
	k: integer, number of samples being compared
	x: double array storing the concatenated samples in the 
        same order as in ns
	ns: integer array storing the k sample sizes, 
         corresponding to x
	zstar: double array storing the l distinct ordered 
             observations in the pooled sample
	L: integer, length of zstar
	
	Outputs:
	when the computation ends, AkN2 and AakN2 are stored in
     the given memory pointed to by adk
*/

void adkTestStat(double *adk, int k, double *x, int *ns){
     	int i;
	int j;

	
	/* fij records the number of observations in the ith 
        sample coinciding with zstar[j], where i = 1, ..., k,
        and j = 1, ..., L                                    */
	/* int fij[k*L]; replaced my next line as per Tomasz Melcer*/
	int *fij = calloc(k*L, sizeof *fij);
	/* lvec is an integer vector with length L, 
		whose jth entry = \sum_{i=1}^{k} f_{ij}, i.e., 
           the multiplicity of zstar[j]                      */
	/* int lvec[L]; replaced my next line as per Tomasz Melcer */
	int *lvec = calloc(L, sizeof *lvec);
	
	/* for computation */
	double mij;
	double maij;
	double innerSum;
	double aInnerSum;
	double bj;
	double baj;
	double tmp;
	
	/* samples is a two-dimensional double array with length
        k; it stores an array of k pointers to double arrays 
        which are the k samples being compared               */
	double **samples;
	
	/* dynamically allocate memory */
	samples = pdmalloc(k);
	nsum = 0;
	for (i = 0; i < k; i++) {
		samples[i] = dmalloc(ns[i]);
		
		for (j = 0; j < ns[i]; j++) {
			samples[i][j] = x[nsum + j];
		}
			
		nsum += ns[i];
	}
	
	/* fij: k*L integer matrix, where L is the length of zstar 
              and k is the number of samples being compared 
	   lvec: integer vector of length L, records the 
               multiplicity of each element of zstar         */	
	for (j = 0; j < L; j++) {
		lvec[j] = 0;
		for (i = 0; i < k; i++) {
			fij[i + j*k] = getCount(zstar[j], samples[i], 
                                                       ns[i]);
			lvec[j] += fij[i + j*k];
		}
	}
	
	adk[0] = adk[1] = 0;
	for (i = 0; i < k; i++) {
		mij = 0;
		maij = 0;
		innerSum = 0;
		aInnerSum = 0;
		
		for (j = 0; j < L; j++) {
			mij += fij[i + j*k];
			maij = mij - (double) fij[i + j*k] / 2.0;
			bj = getSum(lvec, j + 1);
			baj = bj - (double) lvec[j] / 2.0;
			
			if (j < L - 1) {
				tmp = (double) nsum * mij - 
                            (double) ns[i] * bj;
				innerSum = innerSum + 
                                 (double) lvec[j] * tmp * tmp / 
						 (bj * ((double) nsum - bj));
			}
			
			tmp = (double) nsum * maij - 
                                         (double) ns[i] * baj;
			aInnerSum = aInnerSum + 
                            (double) lvec[j] * tmp * tmp / 
					 (baj * (nsum - baj) - 
                             nsum * (double) lvec[j] / 4.0);
		}

		adk[0] = adk[0] + innerSum / ns[i]; /* AkN2*/
		adk[1] = adk[1] + aInnerSum / ns[i]; /* AakN2 */
	}
	
/* k-sample Anderson-Darling test statistics in both original 
   and alternative versions, AkN2 and AakN2, are stored in the 
   given double array adk                                    */
	adk[0] = adk[0] / (double) nsum; /* AkN2*/
	adk[1] = (nsum - 1) * adk[1] / 
               ((double) nsum * (double) nsum); /* AakN2 */
	
	/* free pointers */
	for (i = 0; i < k; i++) {
		free(samples[i]);
	}
	free(samples);
	free(lvec);
	free(fij);
}

 /* *****************************************
  * Project: k-Sample Anderson-Darling Tests
  * functions initvals, exactcomb
  * Last modified: 07.28.2015 (FWS)
  *******************************************/

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */






/* uses recursive backtracking to find all possible
   ways/combinations to divide n elements into k subcollections, 
   where the size of each subcollection is fixed;
   for each combination, the desired test statistics are 
   computed and compared with the observed values.   
   The recursion arises because at each iteration of exactcomb
   a combination of ns[i] is split off from what is left over,
   until nothing more can be split off because only ns[k-1] are
   left over  */
         

     int exactcomb(int now, int *position, int m, int dimst, 
		   void(*testStatFun)(double *teststat, int k, 
                     double *x, int *ns)) {
	/* exactcomb in its initial call starts with now = 0,
   		m = nsum, and position filled with i = 0, 1, 2, ..., nsum-1. */
				
	int i;
	int j;
	
	if (now == k - 1) {
		/* here we have arrived at the end stage when nothing
		   more can be split off */
		
		double xc[nsum];
		double teststatcomb[dimst];
		double *pt;
		/* fills the remaining m=ns[k-1] positions
              		of ivec with now = k-1 (FWS) */
		
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		/* here we equate the pointer pt with that of xc and 
              		by filling the associated array pt we also fill
              		the array xc. (FWS) */
		pt = xc;
		/* here we fill pt=xc first with all the xvec[j]
              		which are supposed to belong to sample 0 according
              		to ivec, then with those belonging to sample 1,
              		and so on. This will give us the full permuted
              		sample sequence (FWS) */
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = xvec[j];
					pt++;
				}
			}
		}
		
					
		/* get test statistics for this combination xc */
		(*testStatFun)(teststatcomb, k, xc, ns);
		
		/* compares test statistics for this combination with
              observed ones */
		for (j = 0; j < dimst; j++) {
			if (getSmat) {
				/* records the test statistics for each 
                         		combination */
				smat[count + ncomb * j] = teststatcomb[j]; 
			}
				/* compares teststatcomb for this combination with 
		   			observed one */
			if (teststatcomb[j] >= teststat[j]) {
				pval[j] = pval[j] + 1;
			}
		}
		
		count++;
		return(2);
		/* this return gets us back to just beyond the point
              		of the last previous call to exactcomb, to find the 
              		next combination at that stage of now (FWS) */
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
			(TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		/* s represents the size of the remainder after the 
              	ns[now] sample values for the current combination 
              	have be chosen. The meaning of the variables a, w, 
              	r is pretty much explained in Knuth, p. 367. 
              	In particular, a[i] = 1 means that the element
              	with index i is designated as part of the chosen
              	combination. (FWS) */
		int r;
		int *a;
		int *w;
		int newposition[s];
		/* this newposition array is meant to replace the
              	position array in the recursive call to exactcomb,
              	to get to the next combination of the k 
              	combinations to be chosen. It is set up below, 
              	right after while (1) {....   (FWS) */

		int *tmp;
		/* this pointer is equated to the pointer of the 
              	newposition array and is used to fill that array. 
              	(FWS) */
		
		/* initializes variables */
		/* this is the start of the C1 step in Knuth's
              	algorithm C (FWS) */
		
		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		/* this is the end of the C1 step in Knuth's 
              		algorithm C (FWS) */
		
		j = r;
		/* the setup of this function assures that j != m at
              		this point since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			/* here we equate the pointers tmp and 
                   		newposition and by filling tmp we fill 
                   		newposition. (FWS) */
			/* visits current combination */
			tmp = newposition;
			/* If indicated by a[i]=1 (relative to the 
                   	current position array and w.r.t. the array a 
                   	in that context), we fill ivec at index
                   	position[i] with the sample index now, that
                   	is under discussion here.
                   	All other position indices are collected
                   	inside the array newposition, by assignment 
                   	via tmp. It amounts to splitting the m 
                   	position elements into two groups of size
                   	ns[now] (the chosen combination for the now 
                   	sample) and s = m-ns[now], the remainder.
                   	(FWS) */
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			/* to get the next combination, as indicated by 
                   	now+1, using the residual position vector 
                   	newposition, but when understanding what 
                   	happens to it, that newposition vector is
                   	referred to as position inside the algorithm
                   	exactcomb. (FWS) */
			exactcomb(now + 1, newposition, s, dimst, testStatFun);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			/* Here we find out whether we have encountered 
                   	the last combination already, and whether we
                   	should step back prior to the last invocation 
                   	of exactcomb, possibly leading to further  
                   	stepping back, until there is no more 
                   	stepping back, i.e., we have traversed all
                   	combination splits. If we do not terminate 
                   	here, we generate the next step in the array
                   	generation, according to Knuth's C2-C7. 
                   	(FWS) */
			
			if (j == m) { 
                     /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		/* This return gets us back to just past the last
              	invocation of exactcomb. We either arrive at 
              	now = k-1 or need to split off further 
              	combinations as needed. (FWS) */ 
		return(0);
	}
}


 /*******************************************
  * Project: k-Sample Anderson-Darling Tests
  * function adkPVal
  * Last modified: 07.21.2015 (FWS)
  *******************************************/



/* computes or estimates (depends on Nsim and the total number
   of all possible combinations) p-values for the observed 
   k-sample Anderson-Darling test statistics in both original 
   and alternative versions for the nonparametric (rank)test 
   described in Scholz F.W. and Stephens M.A. (1987), 
   K-sample Anderson-Darling Tests, Journal of the American 
   Statistical Association, Vol 82, No. 399, pp. 918-924
		
	Arguments:
	pval: double array of length 2, stores estimated p-values 
            for the observed AkN2 and AakN2
	Nsim: integer, number of simulations
	k: integer, number of samples being compared
	x: double array storing the concatenated samples in the 
        same order as in ns
	ns: integer array storing the k sample sizes, 
         corresponding to x
	zstar: double array storing the l distinct ordered 
             observations in the pooled sample
	L: integer, length of zstar
	useExact: integer, 0: not, 1: yes; indicates whether the 
                p-value will be computed via examining all 
                possible combinations	(this occurs when 
                ncomb < Nsim, i.e., the total number of possible 
		     combinations is less than Nsim and the user 
                chooses the exact approach, see R function  
                getAdkPVal for details)
	getA2mat: logical, to indicate if a2mat, a double matrix 
                storing the test statistics of all exact or 
                simulated combinations, will be returned as part 
                of the output	
	ncomb: integer, number of all possible combinations
	a2mat: double matrix, either ncomb * 2 or Nsim * 2, 
             depending on which approach is used, stores the 
             test statistics of all exact or simulated 
             combinations
	
	Outputs:
	when the computation ends, p-values of the observed AkN2 
         and AakN2 are stored in the given memory pointed to by
         pval and the test statistics of all exact or simulated
         combinations are stored in the given memory pointed to
	    by a2mat (1st column: AkN2, 2nd column: AakN2)
*/

void adkPVal(double *pval, int Nsim, int k, double *x, int *ns,
				double *zstar, int L, int useExact,
                      int getA2mat, double ncomb, double *a2mat) {
	int i;
	int j;
	
	int nsum = getSum(ns, k); 
     	/* total sample size = n_1 + ... + n_k */
	double adk[2];
	/* initializes static variables */
	initvals1(k, x, ns, zstar, L);
	
	/* gets observed AkN2 and AakN2 */
	
	(*adkTestStat)(adk, k, x, ns);

	
	pval[0] = pval[1] = 0;
	
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { 
           /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		initvals(k, x, ns, zstar, L, ivec, adk, pval, nsum, 
                    (int) ncomb, getA2mat, a2mat,2);	
		exactcomb(0, position, nsum, 2, adkTestStat);
		
		/* gets exact p-values */
		pval[0] = pval[0] / (double) ncomb;
		pval[1] = pval[1] / (double) ncomb;
		
	
	} else { /* uses Nsim simulations to get p-value */
		double adksim[2];
		double xc[nsum]; /* copy of x */
		
		for (i = 0; i < Nsim; i++) {
			/* gets random permutation xc of x */
			randPerm(nsum, x, xc, ns);

	
		/* gets simulated AkN2 and AakN2 */
		(*adkTestStat)(adksim, k, xc, ns);
		
		/* compares simulated AkN2 and AakN2 with observed 
              ones                                           */
			for (j = 0; j < 2; j++) {
				if (getA2mat) {
		/* records the AD test statistics for each simulated
              combination */
					a2mat[i + Nsim * j] = adksim[j];
				}
				
				if (adksim[j] >= adk[j]) {
					pval[j] = pval[j] + 1;
				}
			}
		}
		
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
		pval[1] = pval[1] / (double) Nsim;
	}
	
	/* finishes using R random number generator */
	PutRNGstate();
	
} 

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Function adkPVal0
 * Last modified: 07.21.2015
 ***************************/

/* wrapper function for function adkPVal to enable R calls 
   to it */


void adkPVal0(double *pval, int *Nsim, int *k, double *x,
              int *ns, double *zstar, int *L, int *useExact,
              int * getA2mat, double *ncomb, double *a2mat){
	adkPVal(pval, *Nsim, *k, x, ns, zstar, *L, *useExact,
                *getA2mat, *ncomb, a2mat);
}

/* wrapper function for function adkTestStat to enable R calls 
   to it */	

void adkTestStat0(double *ans, int *k, double *x, int *ns, double *zstar, int *L){

	initvals1(*k, x, ns, zstar, *L);
	
	adkTestStat(ans, *k, x, ns);
}
/***********************************************
* End Project: k-Sample Anderson-Darling Tests
************************************************/



/***********************************
 * Project:  2*t Contingency Table
 * Function contingency2xtExact
 * Last modified: 07.21.2015
 ***********************************/


/* computes the exact null distribution of the Kruskal-Wallis 
	statistics in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / 
      (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	Avec: integer array of length tnum, storing the column 
            counts with Response "a"
	Bvec: integer array of length tnum, storing the column 
            counts with Response "b"
	tnum: integer, number of columns in the contingency table
	ncomb: integer, number of possible splits of m into the
            sum of tnum nonnegative integers, i.e., 
            choose(m + tnum - 1, tnum - 1)
	results: double vector of length (2 + 2 * ncomb), 
               whose first two entries contain the observed
               value of delta and its p-value, followed by the
               ncomb delta values for each possible split, where
               delta is defined to be \sum (A_i^2 / d_i), 
               followed by the corresponding ncomb probabilities
               (note that some splits have 0 probability since
               a_i <= d_i)
			
	Output:
	the values of delta, \sum (A_i^2 / d_i), observed and for 
      all splits, and their corresponding p-value and
      probabilities are stored in the given memory pointed by 
      results	                                               */

void contingency2xtExact(int *Avec, int *Bvec, int tnum, 
                         int ncomb, int getDist, 
                         double *results) {
     int dvec[tnum]; /* tnum column sums */
	int m = 0; /* row sum for response "a" */
	int n = 0; /* row sum for response "b" */
	int i;
	int j;
	int sum;
	int count;
	int boolean;
	int flag;
	int uvec[tnum - 1]; /* index vector */
	/* xvec: numbers of units with Response "a" in each 
               treatment group                             */
	int xvec[tnum]; 	
	double delta;
	double deltaObserved = 0;
	double prob;
	/* get m, n, dvec, and the observed delta */
	for(i = 0; i < tnum; i++){
        m = m + Avec[i];
		n = n + Bvec[i];
		dvec[i] = Avec[i] + Bvec[i];
		deltaObserved = deltaObserved + Avec[i] * Avec[i] / 
                          (double) dvec[i];
	}
	results[0] = deltaObserved;
	results[1] = 0;

	/* Algorithm using Chase's sequence by Donald Knuth 
		(TAOCP V.4A, 7.2.1.3 Algorithm C);
		goes through all combinations of choosing (tnum-1) 
           from (m + tnum - 1) distinct objects */
	int *a;
	int *w;
	int *pt;
	int mt1 = m + tnum - 1;
	int s = m;
	int r;
	
	/* initializes variables */
	/* this is the start of the C1 step in Knuth's algorithm C 
	   (FWS) */
	a = imalloc(mt1);
	w = imalloc(mt1 + 1);
	
	for (j = 0; j < s; j++) {
		a[j] = 0;
		w[j] = 1;
	}
	
	for (j = s; j < mt1; j++) {
		a[j] = w[j] = 1;
	}
	w[mt1] = 1;
	
	if (s > 0) {
		r = s;
	} else {
		r = tnum - 1;
	}
	/* this is the end of the C1 step in Knuth's algorithm C 
	   (FWS) */
	
	j = r;
	count = 2;
	boolean = 1;
	
	while (boolean) {
		/* visits current combination */
		/* sets up the (tnum - 1) indices, where
		 	1 <= uvec[0] < ... < uvec[tnum - 2] <= mt1 */
		pt = uvec;
		for (i = 0; i < mt1; i++) {
			if (a[i]) {
				*pt = i + 1;
				pt++;
			}
		}
		
		/* computes x_i's , the number of units with response 
              "a" in each of the treatment groups, which are 
               stored in xvec                                  */
		sum = 0;	
		delta = 0;
		for (i = 0; i < tnum - 1; i++) {
			xvec[i] = uvec[i] - i - 1 - sum;
			sum = sum + xvec[i]; 
			delta = delta + xvec[i] * xvec[i] / 
                                        (double) dvec[i];
		}
		
		xvec[tnum-1] = m - sum;	
		delta += xvec[tnum-1] * xvec[tnum-1] / 
                                   (double) dvec[tnum-1];
        /* store delta for this split*/
		if(getDist){
			results[count] = delta;
		}
		
		/* computes the probability associated with current 
               combination */
		prob = 1; /* initializes probability */
		flag = 1;
		i = 0;
		while (flag && i < tnum) {
			if (xvec[i] > dvec[i]) {
				prob = 0;
				flag = 0; 
                       /* gets out of while loop early */
			} else {
				prob *= choose(dvec[i], xvec[i]); 
                      i++;
			}
		}
		
		if (flag) {
			prob = prob / choose(m + n, m);
		}
		/* updating p-value */
		if(delta >= deltaObserved) results[1] += prob;
		/* store probability for this split */
		if(getDist){
			results[count + ncomb] = prob;
		}
		count++;
		/* end of visiting current combination */
		
		/* finds j and branches */
		j = r;
		while(w[j] == 0) {
			w[j] = 1;
			j++;
		}

		if (j == mt1) { 
                 /* terminate point of this algorithm */
			boolean = 0; /* gets out of while loop */
		} else { /* continue */
			w[j] = 0;
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;

					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}

				} else {
					a[j-1] = 1;
					a[j] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}

			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;

					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}

				} else {
					a[j] = 1;
					a[j-1] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
		}
			
	} 	/* end of while(boolean) */
	
	/* frees the pointers */
	free(a);
	free(w);

} 


/* wrapper function for function table2xtExact to enable R 
   calls to it */

void contingency2xtExact0(int *Avec, int *Bvec, int *tnum, 
		int *ncomb, int *getDist, double *ans){
	       contingency2xtExact(Avec, Bvec, *tnum, 
			*ncomb, *getDist, ans);	
}

/* simulates the null distribution of the Kruskal-Wallis 
   statistics in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / 
                                                       (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	dvec: integer array of length tnum, storing the column 
            totals of the 2 x t contingency table
	m: integer, the total number of units with Response "a"
	n: integer, the total number of units with Response "b"
	ncomb: integer, number of possible splits of m into the 
             sum of tnum nonnegative integers, i.e., 
             choose(m + tnum - 1, tnum - 1)
	nsim: integer, number of simulations
	results: double array, storing the simulated values of
              delta, where delta = \sum (A_i^2 / d_i)
			
	Output:
	the simulated values of delta, \sum (A_i^2 / d_i), are
     stored in the given memory pointed to by results	
*/

void contingency2xtSim(int *Avec, int *Bvec, int tnum, int nsim, 
		            int getDist, double *results) {
     int dvec[tnum]; /* tnum column sums */
	int m = 0; /* row sum for response "a" */
	int n = 0; /* row sum for response "b" */
	int i;
	int j;
	int a;
	int nb;
	int k;
	int sum;
	double delta;
	double deltaObserved = 0;
	int pval = 0;
	for(i = 0; i < tnum; i++){
        m = m + Avec[i];
		n = n + Bvec[i];
		dvec[i] = Avec[i] + Bvec[i];
		deltaObserved = deltaObserved + Avec[i] * Avec[i] / 
                                               (double) dvec[i];
	}
	results[0] = deltaObserved;
	results[1] = 0;

	
	/* uses R random number generator */
	GetRNGstate();
	
	for (i = 0; i < nsim; i++) {
		/* initializes variables */
		nb = m + n;
		k = m;
		delta = 0;
		sum = 0;
		
		for (j = 0; j < tnum - 1; j++) {
			nb = nb - dvec[j];
			/* function rhyper in Rmath.h:
			 	random generation for the hypergeometric
                      distribution                           */
			a = (int) rhyper(dvec[j], nb, k); 
			delta = delta + a * a / (double) dvec[j];
			sum = sum + a;
			k = k - a;
		}
		a = m - sum;
		delta = delta + a * a / (double) dvec[tnum - 1];
		if(delta >= deltaObserved) pval +=  1;
		if(getDist){
			results[i+2] = delta;
		}
	}
	results[1] = (double) (pval) / (double) nsim;
	
	/* finishes using R random number generator */
	PutRNGstate();
	
}

/* wrapper function for function contingency2xtSim to enable R calls it */
	


void contingency2xtSim0(int *Avec, int *Bvec, int *tnum, 
		int *nsim, int *getDist, double *ans){
	       contingency2xtSim(Avec, Bvec, *tnum, 
			*nsim, *getDist, ans);	
}
/***********************************
 * End of  2*t Contingency Table 
 * Last modified: 07.21.2015
 ***********************************/




/***************************
 * Project: k-Sample Kruskal-Wallis Test
 * Functions QNTestStat, QNinitvals, QNexact
 * based on KWexact by Angie Zhu 
 * modified 06/11/2015 Fritz Scholz
 ***************************/

/* static int getQNdist; 
            logical value to indicate whether generated QN
               values are recorded in QNvec               
static double *QNvec; */


/* initializes static variables */
void QNinitvals(int kk, double *xxvec, int *nns, 
			int *iivec, double *tteststat, double *ppval, 
                int nnsum, int ggetQNdist, double *QQNvec, int ddimst) {
	k = kk;
	xvec = xxvec;
	ns = nns;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	getSmat = ggetQNdist;
	smat = QQNvec;
	count = 0;
	dimst = ddimst;
}

void QNinitvals1(int kk, double *xxvec, int *nns){
	k = kk;
	xvec = xxvec;
	ns = nns;
}


/* computes the non-normalized k-sample rank score test 
   statistics 
	
	Arguments:
	QN: double array with length 1, stores the non-normalized
         QN test statistic
	k: integer, number of samples being compared
	rx: double array storing the concatenated average rank
         scores of the k samples in the same order as in ns
	ns: integer array storing the k sample sizes, 
          corresponding to rx
   nsum: integer, number of all observations ns[0]+...+ns[k-1]

	Outputs:
	when the computation ends, the non-normalized QN statistic 
      is stored in the given memory pointed to by QN
*/

void QNTestStat(double *QN, int k, double *rx, int *ns) {
    int i;
    int j;
    double Ri;
    int istart, iend;
    QN[0] = 0.0;
    istart = 0;
    for(i = 0; i < k; i++){
		Ri = 0.0;
		iend = istart+ns[i];
		for(j = istart; j < iend;j++){
			Ri += rx[j];
		}
		QN[0] += Ri * Ri / ns[i];
		istart = iend;
	}
    	QN[0] = round(1e8*QN[0])/1e8; 
        /* this avoids computational quirks due to machine
           representation of numbers*/
}


/***************************
 * Project: k-Sample QN Test
 * Function QNpvalue
 * 06/11/2015 Fritz Scholz
 ***************************/

/* for random number generator in R */


/* estimates p-values for the observed k-sample Kruskal-Wallis 
   test statistics.
	
	Arguments:
	pval: double array with length 1, storing the estimated 
            p-value for the observed QN value
	Nsim: integer, number of simulations
	k: integer, number of samples being compared
	rx: double array storing the average rank scores of the
          concatenated samples in the same order as in ns
	ns: integer array, storing the k sample sizes, 
          corresponding to rx
	useExact: integer, 0: not, 1: yes; indicates if the 
                p-value will be computed via examining all
                possible combinations (this occurs when ncomb 
                < Nsim, i.e., the total number of possible
                combinations is less than Nsim and the user
                chooses the exact approach, see R function
                getQNPVal for details)
	getQNdist: logical, to indicate whether the exact or 
                 simulated QNvec will be returned as part of the
                 output
 	ncomb: double, number of all possible combinations
 
		
		
	Outputs:
	when the computation ends, p-values of the observed, 
      non-normalized QN is stored in the memory pointed at by
      pval and the distribution of the non-normalized QN values
      of all exact or simulated combinations is stored in array
      QNvec. The observed non-normalized QN is stored in memory
      pointed at by QNobs.
 */

void QNpvalue(double *pval, int Nsim, int k, double *rx, 
              int *ns, int useExact, int getQNdist, 
              double ncomb, double *QNobs, double *QNvec) {
    
	int i;
     int isim;

	int nsum = getSum(ns, k); 
      	/* total sample size = n_1 + ... + n_k */
	QNinitvals1(k,rx,ns);

	/* get observed test statistic for the average rank score 
        vector rx in non-normalized form */
 	QNTestStat(QNobs, k, rx, ns);
	pval[0] = 0.0;
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { 
             /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
           QNinitvals(k, rx, ns, ivec, QNobs, pval, nsum, 
			   getQNdist, QNvec,1);
		exactcomb(0,position,nsum,1,QNTestStat);

		
		/* gets exact p-values */
		pval[0] = pval[0] / ncomb;
	} else { /* uses Nsim simulations to get p-value */
		double QNsim[1];
		double rc[nsum]; /* copy of rx */
		
		for (isim = 0; isim < Nsim; isim++) {
			randPerm(nsum, rx, rc, ns);

			

			/* gets simulated QN */
 			QNTestStat(QNsim, k, rc, ns);
			/* compares simulated QN with observed one */
			if (QNsim[0] >= QNobs[0]) {
					pval[0] = pval[0] + 1.0;
			}
			
			if (getQNdist) {
				QNvec[isim] = QNsim[0];
			}
		}
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 

/***************************
 * Project: k-Sample QN Test
 * Function QNtest
 * Fritz Scholz, last modified: 06.11.2015
 ***************************/

/* wrapper function for function QNpvalue to enable R calls to it */


void QNtest(double *pval, int *Nsim, int *k, double *rx, 
            int *ns, int *useExact, int *getQNdist, 
            double *ncomb, double *QNobs, double *QNvec){
	   QNpvalue(pval,*Nsim,*k,rx,ns,*useExact,*getQNdist,
		             *ncomb,QNobs,QNvec);

}


/****************************  
 * End of k-Sample QN Test
 ***************************/


/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Functions SteelTestStat, Steelinitvals, 
 * Steelinitvals1, Steelexact, Steelpvalue
 * based on KWexact by Angie Zhu 
 * last modified: 06.11.2015
 ***************************/

/* static variables */
static int k; /* number of samples */
static double *xvec; 
  /* vector of averaged rank scores for all observations */
static int *ns; /* sample sizes for the k samples */
static int *ivec; /* this vector will hold the indices 0,
                     1, ..., k-1, indicating the sample
                     associated with the respective positions. 
                     It indicates the final set of combinations
                     (FWS) */
static double *teststat; /* observed Steel statistic
                           (Steel.obs), not normalized */
static double *pval; /* number of cases where Steel.comb >= 
                        (or <=) Steel.obs */
static int nsum; /* total number of observations on all k
                    samples  */
static int alt; /* indicating the type & direction of the test
                   statistic */
static double *mu; /* array of length k-1 for Wilcoxon means for
                      standardization */
static double *tau; /* array of length k-1 for Wilcoxon std. 
                       devs. for standardization */
/* static int getSteeldist; logical value to indicate whether 
                            generated Steel values are recorded
                            in Steelvec*/
static int count; /* keeps track of how many combinations have 
                     been visited */

/* initializes static variables */
void Steelinitvals(int kk, double *xxvec, int *nns, 
		        int *iivec, double *tteststat, double *ppval, 
                   int nnsum, int aalt, double *mmu, 
                   double *ttau, int ggetSteeldist, 
                   double *SSteelvec, int ddimst) {
	k = kk;
	xvec = xxvec;
	ns = nns;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	alt = aalt;
        mu = mmu;
        tau = ttau;
	getSmat = ggetSteeldist;
	smat = SSteelvec;
	dimst = ddimst;
	count = 0;
}






void Steelinitvals1(int kk, double *xxvec, int *nns, 
		         int aalt, double *mmu, double *ttau){
	k = kk;
	xvec = xxvec;
	ns = nns;
	alt = aalt;
        mu = mmu;
        tau = ttau;
}


/* computes the maximum/minimum/max-absolute pairwise 
   standardized Wilcoxon test statistics when comparing k-1
   samples against the same control sample.
	
	Arguments:
	Steel: double array with length 1, stores the standardized
            version of the Steel test statistic
	k: integer, number of samples being compared, including
        the control sample
	rx: double array storing the concatenated double midrank 
         scores of the k samples in the same order as ns, where 
         the first ns[0] are the controls
	ns: integer array storing the k sample sizes, 
         corresponding to rx
    	nsum: integer, number of all observations ns[0]+...
           +ns[k-1]
     alt: integer with values -1, 0 or 1, indicating which type
          of test statistic to compute. For alt = 1 it computes 
          the maximum standardized Wilcoxon statistic when 
          comparing each treatment sample with the control 
          sample. Standardization is done by mu[i-1] and
          tau[i-1] as mean and standard deviation when dealing
          with treatment i.
          For alt = -1 it computes the minimum standardized
          Wilcoxon statistic when comparing each treatment
          sample with the control sample. 
          For alt = 0 it computes the maximum absolute 
          standardized Wilcoxon statistic when comparing each
          treatment sample with the control sample. 
    	mu:  double array of length k-1, holding the means of the 
          Wilcoxon (Mann-Whitney form) test statistics, i.e., 
          mu[i] = ns[0]*ns[i+1]/2, i=0,...,k-2.
	tau: double array of length k-1, holding the standard 
          deviations of the Wilcoxon (Mann-Whitney form) test
          statistics, in the randomization distribution
  	     conditioned on the tie pattern in all the data, as
          expressed in rx.

	Outputs:
	when the computation ends, the Steel statistic based on
     standardized Wilcoxon statistics is stored in the given
     memory pointed to by Steel
*/

void SteelTestStat(double *Steel, int k, double *rx, int *ns){
             
   	int i;
	int j;
        int m;
    	double Ri;
    	int istart, iend;
	istart = ns[0];
	Ri = 0.0;
	for(i = 1; i < k; i++){
		iend = istart+ns[i];
		for(j = istart; j < iend;j++){
			for( m = 0; m < ns[0]; m++){
				if(rx[m] <= rx[j]){
					if(rx[m] == rx[j]){ Ri += .5; }else{ 
								Ri += 1; }
				}				
			}
		}
		Ri =  (Ri-mu[i-1])/tau[i-1];
		if(alt == 1){
			if(i ==1){ 
				Steel[0] = Ri;}else{
					if(Ri > Steel[0]) Steel[0] = Ri;
			}
		}
		if(alt == -1){
			if(i ==1){ 
				Steel[0]= Ri;}else{
					if(Ri < Steel[0]) Steel[0] = Ri;
			}
		}
		if(alt == 0){
			Ri = fabs(Ri);
			if(i ==1){ 
				Steel[0] = Ri;}else{
					if(Ri > Steel[0]) Steel[0] = Ri;
			}
		}
		istart = iend;
		Ri = 0.0;
	}
}






/* estimates p-values for the observed Steel test statistic.
	
	Arguments:
	pval: double array with length 1, storing the estimated 
            p-value for the observed Steel value
	Nsim: integer, number of simulations
	k: integer, number of samples being compared and control, 
        (k-1) comparisons
	rx: double array storing the midranks of the concatenated
         samples in the same order as ns, the control 
         corresponds to ns[0]
	ns: integer array, storing the k sample sizes,
         corresponding to rx
	useExact: integer, 0: not, 1: yes; indicates if the 
                p-value will be computed via examining all
                possible combinations	(this occurs when ncomb 
                <= Nsim, i.e., the total number of possible 
	           combinations is <= Nsim and the user chooses
                the exact approach)
	getSteeldist: logical, to indicate whether the exact or
                simulated Steelvec will be returned as part of
                the output
 	ncomb: double, number of all possible combinations
	alt: integer -1, 1, or 0 indicating which one-sided or
           two-sided statistic to use. 
		alt = 1, use maximum standardized Wilcoxon statistics
		alt = -1,  use minimum standardized Wilcoxon 
                      statistics
		alt = 0, use maximum absolute standardized Wilcoxon
                    statistics
		
		
	Outputs:
	when the computation ends, p-values of the observed
    	Steel statistic is stored in the memory pointed at by pval
     and the distribution of the standardized Steel statistic 
     values of all exact or simulated combination splits is 
     stored in array Steelvec.
    	The observed standardized Steel statistic is stored in 
     memory pointed at by Steelobs.

*/

void Steelpvalue(double *pval, int Nsim, int k, double *rx, 
                 int *ns, int useExact, int getSteeldist, 
		      double ncomb, int alt, double *mu, 
			 double *tau, double *Steelobs, 
                 double *Steelvec) {
    
	int i;
    	int isim;

	int nsum = getSum(ns, k); 
          /* total sample size = n_1 + ... + n_k */
	Steelinitvals1(k, rx, ns, alt, mu, tau);
	/* get observed test statistic for the average rank score
        vector rx in standardized form */
 	SteelTestStat(Steelobs, k, rx, ns);
		
	pval[0] = 0.0;
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { 
           /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		

		/* initializes static variables */
        	Steelinitvals(k, rx, ns, ivec, Steelobs, pval, nsum,
				alt, mu, tau, getSteeldist, Steelvec,1);
		

		exactcomb(0, position, nsum, 1, SteelTestStat );
		/* gets exact p-values */
		pval[0] = pval[0] / ncomb;
	} else { /* uses Nsim simulations to get p-value */
		double Steelsim[1];
		double rc[nsum]; /* copy of rx */
		
		
		for (isim = 0; isim < Nsim; isim++) {
			/* gets random permutation rc of rx */
			randPerm(nsum, rx, rc, ns);

			/* gets simulated Steel */

 			SteelTestStat(Steelsim, k, rc, ns);
				
			/* compares simulated Steel with observed one */
                        if(alt != -1){ 
				if (Steelsim[0] >= Steelobs[0]) {
					pval[0] = pval[0] + 1.0;
				}
			}else{
				if (Steelsim[0] <= Steelobs[0]) {
					pval[0] = pval[0] + 1.0;
				}
			}
			
			if (getSteeldist) {
				Steelvec[isim] = Steelsim[0];
			}
		}
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 


/***************************************************************
 * Project: Steel Multiple Wilcoxon Test
 * Function Steeltest.c
 * last modified: 06.11.2015
 **************************************************************/

/* wrapper function for function Steelpvalue to enable R calls to it */
	

void Steeltest(double *pval, int *Nsim, int *k, double *rx, 
               int *ns, int *useExact, int *getSteeldist, 
		    double *ncomb, int *alt, double *mu, double *tau,
		    double *Steelobs, double *Steelvec){
	Steelpvalue(pval,*Nsim,*k,rx,ns,*useExact,*getSteeldist,
			*ncomb,*alt,mu,tau,Steelobs,Steelvec);

}

/***************************************
* End of Steel Multiple Wilcoxon Test
*****************************************/

/*************************************************************
 * Project: Steel Multiple Wilcoxon Test Confidence Intervals
 * Function SteelTestStatVec, SteelinitvalsVec, SteelexactVec, 
 * SteelVec
 * 06/11/2015
 *************************************************************/
/* computes the pairwise Wilcoxon test statistics when comparing
   k-1 treatment samples against the same control sample.
	
	Arguments:
	SteelVec: double array of length k-1, stores the 
               Mann-Whitney statistics of the k-1 treatment 
               samples when comparing each with the same
               control sample.
	k: integer, number of samples being compared, including
        the control sample
	rx: double array storing the concatenated double midrank 
         scores of the k samples in the same order as ns, where
         the first ns[0] are the controls
	ns: integer array storing the k sample sizes, 
          corresponding to rx
    	nsum: integer, number of all observations ns[0]+...
           +ns[k-1]

	Outputs:
	when the computation ends, the k-1 Mann-Whitney statistics
     are stored in the given memory pointed to by SteelVec
*/

void SteelTestStatVec(double *SteelVec, int k, double *rx, int *ns){
                 
   	int i;
	int j;
        int m;
    	double Ri;
    	int istart, iend;
	istart = ns[0];
	Ri = 0.0;
	for(i = 1; i < k; i++){
		iend = istart+ns[i];
		for(j = istart; j < iend; j++){
			for( m = 0; m < ns[0]; m++){
				if(rx[m] <= rx[j]){
					if(rx[m] == rx[j]){ Ri += .5; }else{ 
								    Ri += 1; }
				}				
			}
		}
		SteelVec[i-1] = Ri;
		istart = iend;
		Ri = 0.0;
	}
}

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */



/* static variables */
static int k; /* number of samples */
static double *xvec; /* vector of averaged rank scores for all 
                        observations */
static int *ns; /* sample sizes for the k samples */
static int *ivec; /* this vector will hold the indices 0,
                     1, ..., k-1, indicating the sample
                     associated with the respective 
                     positions. It indicates the final set of
                     combinations (FWS) */
static int nsum; /* total number of observations on all k
                    samples */
static double *MannWhitneyStats; 
          /* holds the Mann-Whitney statistics for all visited 
             combinations, stacked on top of each other in
             groups of k-1 */
static int count; /* keeps track of how many combinations have 
                     been visited */

/* initializes static variables */
void SteelinitvalsVec(int kk, double *xxvec, int *nns, 
                      int *iivec, int nnsum, 
				 double *MMannWhitneyStats) {
	k = kk;
	xvec = xxvec;
	ns = nns;
	ivec = iivec;
	nsum = nnsum;
	MannWhitneyStats = MMannWhitneyStats;
	count = 0;
}

/* uses recursive backtracking to find all possible 
   ways/combinations to divide nsum elements into k
   subcollections, where the size of each subcollection 
   is fixed; for each combination, the desired test statistics
   are computed */

int SteelexactVec(int now, int *position, int m) {
			
	int i;
	int ix;
	int ccount;
	int j;
 	int k1 = k - 1;

	if (now == k - 1) {
		
		double rc[nsum];
		double teststatscomb[k1];
		double *pt;
		/* fills the remaining m=ns[k-1] positions of ivec
              with now = k-1 (FWS) */
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		/* here we equate the pointer pt with that of rc and
              by filling the associated array pt we also fill 
              the array rc. (FWS) */
		pt = rc;
		/* here we fill pt=rc first with all the xvec[j] 
              which are supposed to belong to sample 0 according
              to ivec, then with those belonging to sample 1,
              and so on. This will give us the full permuted 
              sample sequence (FWS) */
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = xvec[j];
					pt++;
				}
			}
		}
		
					
		/* get test statistic for this combination rc */
		SteelTestStatVec(teststatscomb, k, rc, ns);
		
		/* records the test statistics for each combination 
              */
		ccount = count * k1;
		for( ix = 0; ix < k1; ix++){
			MannWhitneyStats[ccount+ix] = teststatscomb[ix];
		}
		count++;
		return(2);
		/* this return gets us back to just beyond the point
              of the last previous call to SteelexactVec, to 
              find the next combination at that stage of now
              (FWS) */
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
		   (TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		/* s represents the size of the remainder after the
              ns[now] sample values for the current combination
              have been chosen. The meaning of the variables a, w,
              r are pretty much explained in Knuth, p. 367. 
              In particular, a[i] = 1 means that the element
              with index i is designated as part of the chosen 
              combination. (FWS) */
		int r;
		int *a;
		int *w;
		int newposition[s];	
		/* this newposition array is meant to replace the
              position array in the recursive call 
              SteelexactVec, to get to the next combination of 
              the k combinations to be chosen. It is set up
		    below, right after while (1) {....  (FWS) */

		int *tmp; 
		/* this pointer is equated to the pointer of the 
              newposition array and is used to fill that array. 
              (FWS) */
		
		/* initializes variables */
		/* this is the start of the C1 step in Knuth's 
              algorithm C (FWS) */

		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		/* this is the end of the C1 step in Knuth's
              algorithm C (FWS) */
		j = r;
		/* the setup of this function assures that j != m at
              this point since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			/* here we equate the pointers tmp and 
                   newposition and by filling tmp we fill
                   newposition. (FWS) */
			tmp = newposition;
			/* If indicated by a[i]=1 (relative to the 
                   current position array and w.r.t. the array a 
                   in that context), we fill ivec at index 
                   position[i] with the sample index now, that
                   is under discussion here. All other position 
                   indices are collected inside the array 
                   newposition, by assignment via tmp. 
			   It amounts to splitting the m position
                   elements into two groups of size ns[now] (the
                   chosen combination for the now sample) and s
                   = m-ns[now], the remainder. (FWS) */
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			/* to get the next combination, as indicated by
                   now+1, using the residual position vector 
                   newposition, but when understanding what 
                   happens to it, that newposition vector is 
                   referred to as position inside the algorithm
			   SteelexactVec. (FWS) */
			SteelexactVec(now + 1, newposition, s);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			/* Here we find out whether we have encountered
                   the last combination already, and whether we
                   should step back prior to the last invocation
                   of SteelexactVec, possibly leading to further
                   stepping back, until there is no more 
                   stepping back, i.e., we have traversed all
                   combination splits. If we do not terminate
                   here, we generate the next step in the array 
                   generation, according to Knuth's C2-C7. 
                   (FWS) */
			if (j == m) { 
                      /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		/* This return gets us back to just past the last 
              invocation of SteelexactVec. We either arrive at
              now = k-1 or need to split off further
              combinations as needed.
		   (FWS) */ 
		return(0);
	}
}

/* estimates p-values for the observed Steel test statistic.
	
	Arguments:
	Nsim: integer, number of simulations
	k: integer, number of samples being compared and control, 
        (k-1) comparisons
	rx: double array storing the midranks of the concatenated
         samples in the same order as ns, the control
         corresponds to ns[0]
	ns: integer array, storing the k sample sizes,
         corresponding to rx
	useExact: integer, 0: not, 1: yes; 1 indicates that full
         enumeration of all ncomb combination splits is used in
         getting the Mann-Whitney statistics. This should occur
         only when ncomb <= Nsim. Otherwise Nsim random 
         combination splits are used.
		
	Outputs:
	when the computation ends, the double array 
     MannWhitneyStats will contain in stacked form the 
     enumerated or simulated Mann-Whitney statistics in groups
     of k-1. It will have length ncomb*(k-1) or Nsim*(k-1).
*/

void SteelVec(int Nsim, int k, double *rx, int *ns,
		int useExact, double *MannWhitneyStats) {
    
	int i;
	int ix;
	int ccount;
    	int isim;
	int k1=k-1;

	int nsum = getSum(ns, k); 
        /* total sample size = n_1 + ... + n_k */


	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { 
          /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
        	SteelinitvalsVec(k, rx, ns, ivec, nsum, 
                            MannWhitneyStats);
		SteelexactVec(0, position, nsum);
	} else { /* uses Nsim simulations to get p-value */
		double Steelsim[k1];
		double rc[nsum]; /* copy of rx */
		
		
		for (isim = 0; isim < Nsim; isim++) {
			randPerm(nsum, rx, rc, ns);

			

			/* gets simulated Steel */

 			SteelTestStatVec(Steelsim, k, rc, ns);
			ccount = k1*isim;
			for(ix=0; ix<k1; ix++){
				MannWhitneyStats[ccount+ix] = 
                                          Steelsim[ix];
			}
		}
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 


/*************************************************************
 * Project: Steel Multiple Wilcoxon Test Confidence Intervals
 * Function SteelConf, SteelVec
 * 
 * last modified: 06.11.2015
 *************************************************************/

/* wrapper function for function SteelVec to enable R calls to it */
	
void SteelConf(int *Nsim, int *k, double *rx, int *ns,
		    int *useExact, double *MannWhitneyStats){
	SteelVec(*Nsim, *k, rx, ns, *useExact, MannWhitneyStats);

}

/************************************
 * End of Steel Confidence Intervals
 ************************************/




/* convolution function conv */


/* inserts the value xnew at index place k into a table 
   represented by values 'values' and frequencies 'freqs'. */

void insertxp(double xnew, double pnew, int k, int  *Lt, 
              double *values, double *probs){
/* This function takes a set of strictly increasing values 
   values[0] < ... < values[*Lt-1] with accompanying
   probabilities probs[0], ..., probs[*Lt-1] and inserts a given
   new value xnew at the proper index place k, while shifting
   the higher values and probs up by one index position in the 
   respective arrays.
   The array sizes increase from *Lt to *Lt+1. It is assumed that
   there is sufficient space left. The probability associated 
   with the new values[k] = xnew is set to probs[k]=pnew. Here k 
   can be any of the values 0, 1, 2, ..., *Lt.
*/
	int i;
	if(k <= *Lt-1){
		for(i = *Lt-1; i >= k; i--){
			values[i+1] = values[i];
			probs[i+1] = probs[i];
		}
	}
    values[k] = xnew;
	probs[k] = pnew;
	*Lt = *Lt+1;
}

void convaddtotable(double xnew, double pnew, int *Lt, int M, 
			    double *values, double *probs){
/* this function adds a value xnew to a table of ordered    
   values[0] < ... < values[*Lt-1], increasing *Lt by 1 if a new
   value is inserted at values[k] with probs[k] set to pnew, 
   shifting all other values and probs for indices >= k by one.
   If the value is not new, then only the corresponding probs
   value is incremented by pnew. It requires that *Lt < M.
*/
    	int k1, k2, kk, k;
		if(*Lt > 2){
			k1 = 0;
			k2 = *Lt-1;
			if(xnew < values[k1]){
				k = 0;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew > values[k2]){
				k = *Lt;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[k2]){
              	probs[k2] += pnew;
			}else{ while(k1+1 < k2){
						kk = (int) floor((double)(k2-
                                                   k1)/2) +k1;
						if(xnew < values[kk]){
							k2 = kk;
						} else {
							k1 = kk;
						}
					}
				if(xnew == values[k1]){
					probs[k1] += pnew;
				}else{
					k = k2;
					insertxp(xnew,pnew,k,Lt,values,probs);
				}
			}
		}else if(*Lt <= 0){
			values[0] = xnew;
   			probs[0] = pnew;
			*Lt = 1;
		}else if (*Lt ==1){
			if(xnew < values[0]){
				k  = 0;
				insertxp(xnew,pnew, k,Lt,values,probs);
			}else if(xnew > values[0]){
				k  = 1;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[0]){
				probs[0] += pnew;
			}
		}else if(*Lt == 2){
			if(xnew < values[0]){
				k  = 0;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew > values[1]){
				k = *Lt;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[0]){
				probs[0] += pnew;
			}else if(xnew == values[1]){
				probs[1] += pnew;
			}else{
				k = 1;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}
		}
}

/* convolutes two distributions given by supports x1 and x2 and 
   corresponding probabilities p1 and p2, returning x and p as 
   the resulting distribution */

void conv(double *x1, double *p1, int *n1, double *x2, 
          double *p2, int *n2, double *x, double *p, int *n) {
	int i,j, M;
    double xij, pij;
    M = n1[0]*n2[0];
    n[0] = 0;
    for(i=0; i < n1[0]; i++){
		for(j=0; j < n2[0]; j++){
			xij = x1[i]+x2[j];
			pij = p1[i]*p2[j];
			xij = round(1e8*xij)/1e8;
			convaddtotable(xij,pij, n, M, x, p);
		}
	}
}
/* convolutes two vectors x1 and x2 of respective lengths n1 and 
   and n2 and produces the vector of length n = n1*n2 of all
   possible sums x1[i]+x2[j] */

void convvec(double *x1, int *n1, double *x2, int *n2, 
             double *x, int *n) {
    int i,j;
    /* double xij, pij;  */
    n[0] = 0;
    for(i=0; i < n1[0]; i++){
		for(j=0; j < n2[0]; j++){
			x[j+i*n2[0]] = x1[i]+x2[j];
		}
	}
}


void Harding(int k, int L1, int *nn, int *nvec, double *freq){

	int L, M, i, ii, m, n, P, Q, t, u, s;
	double mnm ;
	L = L1-1;
        M = L/2;

	freq[0] = 1;
	for( i = 1; i < L1; i++ ){ freq[i] = 0; }
	for( i = 1; i <= (k-1); i++ ){
		m = nvec[i-1] - nvec[i];
		n = nvec[i];
		
		if( n+1 <= M){
			P = m + n;
			if( M < P ) P = M;
			for( t = n+1; t <= P; t++){
				for( u = M; u >= t; u--){
					freq[u] = freq[u] - freq[u-t];
				}
			}
			
		}
		Q = M;
		if(m < M) Q = m;
		for(s = 1; s <= Q; s++){
			for(u = s; u <= M; u++){
   				freq[u] = freq[u] + freq[u-s];
			}
		}

		mnm = (double) choose(m+n,m);
		for(ii=0; ii < L1; ii++){
		freq[ii] = freq[ii]/mnm;
	}
		
	}
	

 
	
	if( L % 2 == 0 ){
		for( i = 1; i <= M; i++) freq[M + i] = freq[M - i];

	} else {
		for( i = 1; i <= (M+1); i++) freq[M+i] = freq[M+1-i];
	}
}


/* wrapper function for function Harding to enable R calls to it */
	

void Harding0(int *k, int *L1, int *nn, int *nvec, double *freq){
        		Harding(*k, *L1, nn, nvec, freq);
}

/***************************
 * Project: k-Sample Jonckhere-Terpstra Test
 * Functions JTTestStat, JTinitvals, JTexact
 * based on KWexact by Angie Zhu 
 * modified 08/25/2015 Fritz Scholz
 ***************************/

/* static int getJTdist; 
            logical value to indicate whether generated JT
               values are recorded in JTvec               
static double *JTvec; */


/* initializes static variables */
void JTinitvals(int kk, double *xxvec, int *nns, 
			int *iivec, double *tteststat, double *ppval, 
                int nnsum, int ggetJTdist, double *JJTvec, int ddimst) {
	k = kk;
	xvec = xxvec;
	ns = nns;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	getSmat = ggetJTdist;
	smat = JJTvec;
	count = 0;
	dimst = ddimst;
}

void JTinitvals1(int kk, double *xxvec, int *nns){
	k = kk;
	xvec = xxvec;
	ns = nns;
}


/* computes the non-normalized k-sample rank score test 
   statistics 
	
	Arguments:
	JT: double array with length 1, stores the JT test statistic
	k: integer, number of samples being compared
	rx: double array storing the concatenated scores of the k 
		samples in the same order as in ns
	ns: integer array storing the k sample sizes, 
          	corresponding to rx
   nsum: integer, number of all observations ns[0]+...+ns[k-1]

	Outputs:
	when the computation ends, the JT statistic 
      is stored in the given memory pointed to by JT
*/

void JTTestStat(double *JT, int k, double *rx, int *ns) {
   	int i, j;
        int m, n;
    	int mstart, mend, nstart, nend;
	mstart = 0.0;
	JT[0] = 0.0;


	for(i = 0; i < k-1; i++){

		mend = mstart + ns[i];
		nstart = mend;
		for(j = i+1; j < k; j++){
			nend = nstart+ns[j];
			for(n = nstart; n < nend; n++){
				for( m = mstart; m < mend; m++){
					if(rx[m] <= rx[n]){
						if(rx[m] == rx[n]){ 
							JT[0] += .5; }else{ 
							JT[0] += 1; }

					}				
				}
			}

			nstart = nend;

		}
		mstart = mend;

	}
}

/* wrapper function for function JTTestStat to enable R calls to it */
	

void JTTestStat0(double *JT, int *k, double *rx, int *ns){
        		JTTestStat(JT, *k, rx, ns);

}


/***************************
 * Project: k-Sample JT Test
 * Function JTpvalue
 * 08/25/2015 Fritz Scholz
 ***************************/

/* for random number generator in R */


/* estimates p-values for the observed k-sample 
	Jonckheere-Terpstra test statistics.
	
	Arguments:
	pval: double array with length 1, storing the estimated 
            p-value for the observed JT value
	Nsim: integer, number of simulations
	k: integer, number of samples being compared
	rx: double array storing the scores of the
          concatenated samples in the same order as in ns
	ns: integer array, storing the k sample sizes, 
          corresponding to rx
	useExact: integer, 0: not, 1: yes; indicates if the 
                p-value will be computed via examining all
                possible combinations (this occurs when ncomb 
                < Nsim, i.e., the total number of possible
                combinations is less than Nsim and the user
                chooses the exact approach, see R function
                getJTPVal for details)
	getJTdist: logical, to indicate whether the exact or 
                 simulated JTvec will be returned as part of the
                 output
 	ncomb: double, number of all possible combinations
 
		
		
	Outputs:
	when the computation ends, the p-value of the observed, 
      	JT is stored in the memory pointed at by pval and the 
	distribution of the JT values of all exact or simulated 
	combinations is stored in array JTvec, provided getJTdist
	is not 0. 
	The observed JT is stored in memory pointed at by JTobs.
 */

void JTpvalue(double *pval, int Nsim, int k, double *rx, 
              int *ns, int useExact, int getJTdist, 
              double ncomb, double *JTobs, double *JTvec) {
    
	int i;
     	int isim;

	int nsum = getSum(ns, k); 
      	/* total sample size = n_1 + ... + n_k */
	JTinitvals1(k,rx,ns);

	/* get observed test statistic for the score vector rx */
 	JTTestStat(JTobs, k, rx, ns);
	pval[0] = 0.0;
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { 
             /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
           JTinitvals(k, rx, ns, ivec, JTobs, pval, nsum, 
			   getJTdist, JTvec,1);
		exactcomb(0,position,nsum,1,JTTestStat);

		
		/* gets exact p-values */
		pval[0] = pval[0] / ncomb;
	} else { /* uses Nsim simulations to get p-value */
		double JTsim[1];
		double rc[nsum]; /* copy of rx */
		
		for (isim = 0; isim < Nsim; isim++) {
			randPerm(nsum, rx, rc, ns);

			

			/* gets simulated JT */
 			JTTestStat(JTsim, k, rc, ns);
			/* compares simulated JT with observed one */
			if (JTsim[0] >= JTobs[0]) {
					pval[0] = pval[0] + 1.0;
			}
			
			if (getJTdist) {
				JTvec[isim] = JTsim[0];
			}
		}
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 

/***************************
 * Project: k-Sample JT Test
 * Function JTtest
 * Fritz Scholz, last modified: 08.26.2015
 ***************************/

/* wrapper function for function JTpvalue to enable R calls to it */


void JTtest(double *pval, int *Nsim, int *k, double *rx, 
            int *ns, int *useExact, int *getJTdist, 
            double *ncomb, double *JTobs, double *JTvec){
	   JTpvalue(pval,*Nsim,*k,rx,ns,*useExact,*getJTdist,
		             *ncomb,JTobs,JTvec);

}


/****************************  
 * End of k-Sample JT Test
 ***************************/




void randPerm(int nsum, double *rx, double *rc, int *ns){
/* Prior to using function randPerm, need to execute GetRNGstate().
   After being done with all usages of randPerm during a simulation
   run, need to execute PutRNGstate()  */

	double randy;
	double temp;
	int j, index;

	mcopy(rx, rc, nsum, 1); /* gets a copy of rx */

	/* generate a random permutation of rc by randomly interchanging 
	   values on positions nsum - 1, nsum - 2, ..., ns[0] 
	   (C uses 0-based indexing; elements 0, ... , ns[0] - 1 all belong 
           to the first sample;  
           for details of this algorithm, see "Simulation" by Sheldon M. Ross, 
           e.g. 4b, p.51-52.)  */
	for (j = nsum; j > ns[0]; j--) {
		randy = runif(0, 1);
		while(1 <= randy ) { 
                /* to eliminate rare event randy = 1 */
			randy = runif(0, 1);
		}		
		/* index is an random integer between 0 
                        and j-1 (discrete uniform) */
		index = (int) floor(randy * (double) (j));

		/* interchanges the values at positions j-1 and index */
		temp = rc[j-1];
		rc[j-1] = rc[index];
		rc[index] = temp;
	}


}

static const
R_CMethodDef cMethods[] = {
	{"convvec", (DL_FUNC) &convvec, 6},
	{"adkTestStat0", (DL_FUNC) &adkTestStat0, 6},
	{"adkPVal0", (DL_FUNC) &adkPVal0, 11},
	{"contingency2xtExact0", (DL_FUNC) &contingency2xtExact0, 6},
	{"contingency2xtSim0", (DL_FUNC) &contingency2xtSim0, 6},
	{"convC", (DL_FUNC) &conv, 9},
	{"Harding0", (DL_FUNC) &Harding0, 5},
	{"JTtest", (DL_FUNC) &JTtest, 10},
	{"QNtest", (DL_FUNC) &QNtest, 10},
	{"Steeltest", (DL_FUNC) &Steeltest, 13},
	{"SteelConf", (DL_FUNC) &SteelConf, 6},
	{NULL, NULL, 0}
};

void R_init_kSamples(DllInfo *info)
{
    	R_registerRoutines(info, 
			cMethods, 
			NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);

}



