/*
**  Person-years calculations, leading to expected survival for a cohort.
**    The output table depends only on factors, not on continuous.
**    This version converted to .Call syntax for memory savings
**
**  Input:
**
**    expected table, a multi-way array
**      efac[edim]  1=is a factor, 0=continuous (time based)
**      edims[edim] the dimension vector of the table; edim is its length
**      ecut[sum(edims)]  the starting point (label) for each dimension.
**                          if it is a factor dim, will be 1:edims[i]
**      expect      the actual table of expected rates
**
**    subject data
**
**      x[edim, n]  where each subject indexes into the expected table
**                       at time 0, n= number of subjects
**      y[n]         the time at risk for each subject
**		status[n]    the status for each subject
**
**    control over output
**
**      times[ntime]    the list of output times
**
**    Output
**
**
*/
#include <math.h>
#include "survprotomoj.h"

/* my habit is to name a S object "charlie2" and the pointer
**  to the contents of the object "charlie"; the latter is
**  used in the computations
*/
SEXP expc(SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2) {
    int i,k;
    int     n,	    edim;
    double  **x;
    double  *data2;
    double  **ecut, *etemp;
    double  hazard;						   /*cum hazard over an interval */
    double     	    etime,	    et2;
    int   indx,
	    indx2;
    double  wt;

    int	    *efac, *edims;
    double  *expect, *y ;
    SEXP      rlist, rlistnames;

	/*my declarations*/

    SEXP    si2;
    double  *si;


    /*
    ** copies of input arguments
    */

    efac  = INTEGER(efac2);
    edims = INTEGER(edims2);
    edim  = LENGTH(edims2);
    expect= REAL(expect2);


    n     = LENGTH(y2);									/*number of individuals */
    x     = dmatrix(REAL(x2), n, edim);
    y     = REAL(y2);									/*follow-up times*/

    /* scratch space */
    data2 = (double *)ALLOC(edim+1, sizeof(double));

    /*si2 = (double *)ALLOC(n, sizeof(double));			/*Si for each individual - a je to prav???*/

	/*
    ** Set up ecut index as a ragged array
    */
    ecut = (double **)ALLOC(edim, sizeof(double *));
    etemp = REAL(ecut2);
    for (i=0; i<edim; i++) {
	ecut[i] = etemp;
	if (efac[i]==0)     etemp += edims[i];
	else if(efac[i] >1) etemp += 1 + (efac[i]-1)*edims[i];
	}

    /*
    ** Create output arrays
    */

    PROTECT(si2 = allocVector(REALSXP, n));					/* Si for each individual*/
    si = REAL(si2);


 	/*initialize Si values*/
    for (i=0; i<n; i++) {
   	si[i] =1;
   	}

    /* compute  for each individual*/
    for (i=0; i<n; i++) {
	/*
	** initialize
	*/
	for (k=0; k<edim; k++)	data2[k] = x[k][i];						/* the individual's values of demographic variables at time 0 */

		/*
	** add up hazard
	*/


	    /* expected calc
	    **  The wt parameter only comes into play for older style US rate
	    **   tables, where pystep does interpolation.
	    ** Each call to pystep moves up to the next 'boundary' in the
	    **  expected table, data2 contains our current position therein
	    */
	    etime = y[i];
	    hazard =0;
	    while (etime >0) {
		et2 = pystep(edim, &indx, &indx2, &wt, data2, efac,
			     edims, ecut, etime, 1);
		if (wt <1) hazard+= et2*(wt*expect[indx] +(1-wt)*expect[indx2]);
		else       hazard+= et2* expect[indx];
		for (k=0; k<edim; k++)
		    if (efac[k] !=1) data2[k] += et2;
		etime -= et2;

		}
		si[i] =  exp(-hazard);
	}

    /*
    ** package the output
    */
    PROTECT(rlist = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(rlist,0, si2);


    PROTECT(rlistnames= allocVector(STRSXP, 1));
    SET_STRING_ELT(rlistnames, 0, mkChar("surv"));


    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(3);					/*kolk mora bit tu stevilka??  kolikor jih je +2 (rlist, rlistnames)*/
    return(rlist);
    }
