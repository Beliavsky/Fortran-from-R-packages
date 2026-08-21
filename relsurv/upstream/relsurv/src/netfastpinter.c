/*
**  calculation of various quantities needed for the rs.surv function (for PP method and Ederer II method) - sums over individuals at each time
**
**    This version converted to .Call syntax for memory savings
**
**  Input:
**
**
**      efac[edim]  1=is a factor, 0=continuous (time based)  (edim is the number of variables in population mortality tables, usually 3 (age,sex,year), efac tells if they change in time, usually 1,0,1 (age and year change, sex does not))
**      edims[edim] the dimension vector of the population mortality table; edim is its length   (for example 111, 2, 40  : 111 ages, 2 sexes, 40 years)
**      ecut[sum(edims)]  the starting point (label) for each dimension, if factor variable, then NULL.
**							for example, for age:   0.00,   365.24,   730.48,  1095.72,  1460.96 ...
**      expect      the actual population mortality table (values - hazards per day)
**
**    subject data
**
**      x[edim, n]  where each subject indexes into the population mortality table  at time 0, n= number of subjects: a matrix - one row per individual - his value of age, sex and year at time of diagnosis
**      y[n]         the time at risk (follow-up time) for each subject
**		status[n]    the status for each subject: 0 (censored) or 1 (death)
**
**    Output
**
**	  dnisi:    sum(dNi/Spi) at each follow-up time
**	  yisi:    sum(Yi/Spi) at each follow-up time
**	  yidlisi:    sum(YidLambdapi/Spi) at each follow-up time
**	  dnisisq:    sum(dNi/Spi^2) at each follow-up time - needed for the variance
**	  yi:     sum(Yi) at each follow-up time - number at risk at that time
**	  dni:    sum(dNi) at each follow-up time - number of events at that time
**    yidli:    sum(YidLambdapi/Spi) at each follow-up time
**
*/
#include <math.h>
#include "survprotomoj.h"

/* using thernau's habit: name a S object "charlie2" and the pointer
**  to the contents of the object "charlie"; the latter is
**  used in the computations
*/
SEXP netfastpinter(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2,SEXP ys2, SEXP status2,     SEXP times2) {
    int i,j,k;
    int     n,
	    edim,
	    ntime;
    double  **x;
    double  *data2, *si, *sitt;
    double  **ecut, *etemp;
    double  hazard, hazspi;						   /*cum hazard over an interval, also weigthed hazard */
    double     thiscell,
	    etime,
	    time,
	    et2;
    int   indx,
	    indx2;
    double  wt;

    int	    *efac, *edims, *status;
    double  *expect, *y,*ys, *times;
    SEXP      rlist, rlistnames;

	/*my declarations*/

    SEXP    yidli2, dnisi2,yisi2,yidlisi2,yi2,dni2,dnisisq2,yisitt2,yidlisitt2,yidlisiw2;
    double  *yidli, *dnisi,*yisi,*yidlisi,*yi,*dni,*dnisisq, *yisitt,*yidlisitt,*yidlisiw;


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
    ys	  = REAL(ys2);
    status = INTEGER(status2);								/* status */
    times = REAL(times2);
    ntime = LENGTH(times2);								/*length of times for reportint */

    /* scratch space */
    data2 = (double *)ALLOC(edim+1, sizeof(double));

    si = (double *)ALLOC(n, sizeof(double));			/*Si for each individual - this is a pointer, the values are called using s[i]*/
	sitt = (double *)ALLOC(n, sizeof(double));			/*Si at the beg. of the interval for each individual */
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

	PROTECT(yidli2 = allocVector(REALSXP, ntime));		/*sum Yi dLambdai for each time* - length=length(times2)*/
    yidli = REAL(yidli2);
	PROTECT(dnisi2 = allocVector(REALSXP, ntime));		/*sum dNi/Si for each time* - length=length(times2)*/
    dnisi = REAL(dnisi2);
    PROTECT(yisi2 = allocVector(REALSXP, ntime));		/*sum Yi/Si for each time* - length=length(times2)*/
    yisi = REAL(yisi2);
    PROTECT(yisitt2 = allocVector(REALSXP, ntime));		/*add tt*/
    yisitt = REAL(yisitt2);
    PROTECT(yidlisi2 = allocVector(REALSXP, ntime));		/*sum yi/Si dLambdai for each time* - length=length(times2)*/
    yidlisi = REAL(yidlisi2);
    PROTECT(yidlisitt2 = allocVector(REALSXP, ntime));		/*add tt*/
    yidlisitt = REAL(yidlisitt2);
    PROTECT(yidlisiw2 = allocVector(REALSXP, ntime));		/*add w*/
    yidlisiw = REAL(yidlisiw2);
	PROTECT(yi2 = allocVector(REALSXP, ntime));					/* sum yi at each time*/
    yi = REAL(yi2);
    PROTECT(dni2 = allocVector(REALSXP, ntime));		/*sum Yi dLambdai for each time* - length=length(times2)*/
    dni = REAL(dni2);
    PROTECT(dnisisq2 = allocVector(REALSXP, ntime));		/*sum yi/Si dLambdai for each time* - length=length(times2)*/
    dnisisq = REAL(dnisisq2);


 	/*initialize Si values*/
    for (i=0; i<n; i++) {
   	si[i] =1;
   	sitt[i] =1;
   	}


	/*initialize output values*/
    for (j=0; j<ntime; j++) {
	yidli[j] =0;
	dnisi[j] =0;
	yisi[j]=0;
	yisitt[j]=0;
	yidlisi[j]=0;
	yidlisitt[j]=0;
	yidlisiw[j]=0;
	yi[j]=0;
	dni[j]=0;
	dnisisq[j]=0;
	}

	time =0;
	for (j=0; j<ntime ; j++) {			/* loop in time */


    thiscell = times[j] - time;

    /* compute  for each individual*/
    for (i=0; i<n; i++) {
	if(y[i]>= times[j]){				// if still at risk
	/*
	** initialize
	*/
	for (k=0; k<edim; k++){
		data2[k] = x[k][i];						/* the individual's values of demographic variables at time 0 */
		if (efac[k] !=1) data2[k] += time;  /* add time to time changing variables */
	}

	/*
	** add up hazard
	*/


	    /* expected calc
	    **  The wt parameter only comes into play for older style US rate
	    **   tables, where pystep does interpolation.
	    ** Each call to pystep moves up to the next 'boundary' in the
	    **  expected table, data2 contains our current position therein
	    */
	    etime = thiscell;
	    hazard =0;
	    hazspi=0;					//integration of haz/si
	    while (etime >0) {
		et2 = pystep(edim, &indx, &indx2, &wt, data2, efac,
			     edims, ecut, etime, 1);
		hazspi+= et2* expect[indx]/(si[i]*exp(-hazard));		//add the integrated part
		if (wt <1) hazard+= et2*(wt*expect[indx] +(1-wt)*expect[indx2]);
		else       hazard+= et2* expect[indx];
		for (k=0; k<edim; k++)
		    if (efac[k] !=1) data2[k] += et2;
		etime -= et2;
		}
		sitt[i] = si[i];				// si at the beginning of the interval
		si[i] = si[i]* exp(-hazard);


		if(ys[i]<times[j]){		// if start of observation before this time  (change in 224 from <=)
			yisi[j]+=1/si[i];
			yisitt[j]+=1/sitt[i];
			yidlisi[j]+=hazard/si[i];
			yidlisitt[j]+=hazard/sitt[i];
			yidlisiw[j]+=hazspi;
			yidli[j]+=hazard;
			yi[j]+=1;
			if(y[i]==times[j]){
				dnisi[j]+=status[i]/si[i];
				dni[j]+=status[i];
				dnisisq[j]+=status[i]/(si[i]*si[i]);
			}	// if this person died at this time
		  } // if start of observation before this time
	   	  } // if still at risk
	    }// loop through individuals
	    time  += thiscell;
	}// loop through times

    /*
    ** package the output
    */
    PROTECT(rlist = allocVector(VECSXP, 10));					//number of variables
  	SET_VECTOR_ELT(rlist,0, dnisi2);
 	 SET_VECTOR_ELT(rlist,1, yisi2);
 	 SET_VECTOR_ELT(rlist,2, yidlisi2);
 	 SET_VECTOR_ELT(rlist,3, dnisisq2);
 	 SET_VECTOR_ELT(rlist,4, yi2);
 	 SET_VECTOR_ELT(rlist,5, dni2);
 	 SET_VECTOR_ELT(rlist,6, yidli2);
 	 SET_VECTOR_ELT(rlist,7, yisitt2);							/*added tt*/
	 SET_VECTOR_ELT(rlist,8, yidlisitt2);						/*added tt*/
	 SET_VECTOR_ELT(rlist,9, yidlisiw2);						/*added w*/




    PROTECT(rlistnames= allocVector(STRSXP, 10));					//number of variables
    SET_STRING_ELT(rlistnames, 0, mkChar("dnisi"));
    SET_STRING_ELT(rlistnames, 1, mkChar("yisi"));
    SET_STRING_ELT(rlistnames, 2, mkChar("yidlisi"));
 	SET_STRING_ELT(rlistnames, 3, mkChar("dnisisq"));
    SET_STRING_ELT(rlistnames, 4, mkChar("yi"));
    SET_STRING_ELT(rlistnames, 5, mkChar("dni"));
 	SET_STRING_ELT(rlistnames, 6, mkChar("yidli"));
    SET_STRING_ELT(rlistnames, 7, mkChar("yisitt"));				/*added tt*/
    SET_STRING_ELT(rlistnames, 8, mkChar("yidlisitt"));				/*added tt*/
    SET_STRING_ELT(rlistnames, 9, mkChar("yidlisiw"));				/*added w*/



    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(12);					/*number of variables + 2*/
    return(rlist);
    }
