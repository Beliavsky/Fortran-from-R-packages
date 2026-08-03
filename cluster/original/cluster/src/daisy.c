/* Produced by f2c and Martin's  f2c-clean,v 1.11 2012/05/04 --
   and some pretty editing
*/
#include <math.h>

#include <R.h>
#include <Rinternals.h>
/* -> Rconfig.h, but also Boolean.h RS.h */

/* called _only_ once from ../R/daisy.q
   currently as  .C(cl_daisy, ........)
*/
void cldaisy(int *nn, int *jpp, double *x,
	     double *valmd, /* valmisdat[1:p] */
	     double *weights,/* weights[1:p */
	     int *jtmd, int *jdat,
	     int *vtype,
	     int *ndyst,
	     int *mdata,
	     double *disv)
{
/* Calculating dissimilarities between objects or variables

            nn  = number of objects
            jpp = number of variables used for the calculations

   The following vectors and matrices must be dimensioned in the main program :

       vtype was character originally
       vtype(j) is the type of variable j:
              = 1 (A) for an Asymmetric binary variable
              = 2 (S) for a  Symmetric  binary variable
              = 3 (N) for a  Nominal  variable
              = 4 (O) for an Ordinal  variable
              = 5 (I) for an Interval variable (additive)
              = 6 (T) for a  raTio    variable (log transformed)
     ndyst is the "kind of dissimilarity/distance"  aka 'diss_type'
              = 0  "mixed" / gower == here treated _identically_ to "manhattan" (L1) !
              = 1  "euclidean" (= L2 )
              = 2  "manhattan" (= L1 )
       jtmd[ ] is only read if there are missing values : if(mdata)
       jtmd(j) =  0 if variable j is binary
               = -1 if variable j is not binary and has missing values
               = +1 if variable j is not binary and has no missing values
*/

    int x_dim1 = *nn;
    Rboolean hasna = *mdata != 0;

    /*         calculation of the dissimilarities */
    int nlk = 0;
    if (*jdat == 1) {

/* Case I: `mixed' type variables  <==> ndyst == 0 */

	for (int l = 1; l < *nn; ++l) {
	    for (int k = 0; k < l; ++k) {
		/*  disv(nlk) :=  d(l,k) := Dissimilarity between obs. l & k := dlk/ppa */
		double ppa = 0., dlk = 0.;
		for (int j = 0; j < *jpp; ++j) {
		    if (vtype[j] >= 3) {
			if (hasna) {
			    if (jtmd[j] < 0) {
				if ((x[l + j * x_dim1] == valmd[j]) ||
				    (x[k + j * x_dim1] == valmd[j])) {
				    continue; // goto L420;
				}
			    }
			}
			ppa += weights[j];
			if (vtype[j] == 3) { /* type = "N"ominal */
			    if (x[l + j * x_dim1] != x[k + j * x_dim1]) {
				dlk += weights[j];
			    }
			} else { /* type in {"O", "I", "T"} */
			    dlk += weights[j] * fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
			}
		    } else { /* vtype(j) in {1, 2} <---> binary variable x(*,j) */
			if((x[l + j * x_dim1] != 0.f && x[l + j * x_dim1] != 1.f) ||
			   (x[k + j * x_dim1] != 0.f && x[k + j * x_dim1] != 1.f)) {
			    continue; // goto L420;
			}
			if (vtype[j] == 2 || x[l + j * x_dim1] != 0. || x[k + j * x_dim1] != 0.) {
			    ppa += weights[j];
			}
			if (x[l + j * x_dim1] != x[k + j * x_dim1]) {
			    dlk += weights[j];
			}
		    }
		    /* L420: -- continue -- */
		} // for(j ..)
		if (ppa == 0.f) {
		    disv[nlk] = -1.;
		} else {
		    disv[nlk] = dlk / ppa;
		}
		++nlk;
	    } // for( k ..)
	} // for( l ..)

    } else {

/* Case II : jdat != 1:  all variables are interval scaled
   -------   ~~~~~~~~~ { basically === dysta() in ./dysta.f
                         FIXME: common code! } */
	double pp = (double) (*jpp);
	for (int l = 1; l < *nn; ++l) {
	    for (int k = 0; k < l; ++k) {
		double clk = 0.;
		int npres = 0;
		/*  disv(nlk) :=  d(l,k) := Dissimilarity between obs.  l & k */
		for (int j = 0; j < *jpp; ++j) {
		    if (hasna) {
			if (jtmd[j] < 0) {
			    if ((x[l + j * x_dim1] == valmd[j]) ||
				(x[k + j * x_dim1] == valmd[j])) {
				continue; // goto L530;
			    }
			}
		    }
		    ++npres;
		    if (*ndyst == 1) { /* euclidean */
			clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) *
			       (x[l + j * x_dim1] - x[k + j * x_dim1]);
		    } else { /* manhattan */
			clk += fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
		    }
		    /* L530: ; */
		}
		double rpres = (double) npres;
		if (npres == 0) { /* all missing */
		    disv[nlk] = -1.f;
		} else if (*ndyst == 1) {
		    disv[nlk] = sqrt(clk * (pp / rpres));
		} else {
		    disv[nlk] = clk * (pp / rpres);
		}
		++nlk;
	    } // for( k ..)
	} // for( l ..)
    }
} /* cldaisy */

