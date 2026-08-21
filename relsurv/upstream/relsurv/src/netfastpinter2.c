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
 **	  dci:    sum(dCi) at each follow-up time - number of censorings at that time
 **    yidli:    sum(YidLambdapi/Spi) at each follow-up time
 **
 */
#include <math.h>
#include "survprotomoj.h"

/* using thernau's habit: name a S object "charlie2" and the pointer
 **  to the contents of the object "charlie"; the latter is
 **  used in the computations
 */
SEXP netfastpinter2(   SEXP   efac2,   SEXP edims2,
                       SEXP   ecut2,     SEXP   expect2,
                       SEXP   x2, 	SEXP   y2,SEXP ys2, SEXP status2,     SEXP times2, SEXP myprec2) {
  int i,j,k,jfine;
  int     n,
  edim,
  ntime,
  nprec;
  double  **x;
  double  *data2, *si, *sitt;
  double  **ecut, *etemp;
  double  hazard;						   /*cum hazard over an interval */
double     thiscell,
time,
et2,
fyisi,										/* fyisi and fyidlisi are the values in the finer division of the interval, ftime is the tiny time in those intervals */
fyidlisi,
fyidlisi2,
fyisi2,
ftime,
fthiscell,
fint,
sisum,
sisumtt,
lambdapi,
lambdapi2,
timestart;
int   indx,
indx2;
double  wt;

int	    *efac, *edims, *status;
double  *expect, *y,*ys, *times, *myprec;
SEXP      rlist, rlistnames;

/*my declarations*/

SEXP    yidli2, dnisi2,yisi2,yidlisi2,yi2,dni2,dci2,dnisisq2,yisitt2,yidlisitt2,yidlisiw2;
double  *yidli, *dnisi,*yisi,*yidlisi,*yi,*dni,*dci,*dnisisq, *yisitt,*yidlisitt,*yidlisiw;


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
myprec = REAL(myprec2);
//nprec = LENGTH(myprec);



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
PROTECT(dci2 = allocVector(REALSXP, ntime));
dci = REAL(dci2);
PROTECT(dnisisq2 = allocVector(REALSXP, ntime));		/*sum yi/Si dLambdai for each time* - length=length(times2)*/
dnisisq = REAL(dnisisq2);


/*initialize Si values*/
for (i=0; i<n; i++) {
  si[i] =1;
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
  dci[j]=0;
  dnisisq[j]=0;
}

time =0;
timestart=0;

//for (j=0; j<nprec ; j++) {			/* loop in time */
//fthiscell = myprec[j];
//}

for (j=0; j<ntime ; j++) {			/* loop in time */
//for (j=0; j<2 ; j++) {

thiscell = times[j] - time;

  /* add an additional, tinier division for integral calculation. Keep values only at the end of the less fine division (j). For now, precision is fixed to 0.1*/

  ftime=0;
  fyisi=0;
  fyidlisi=0;
  fyisi2=0;
  fyidlisi2=0;
  hazard=0;
  jfine=0;
  fint=0;
  /*initialize Sitt values*/
  for (i=0; i<n; i++) {
    sitt[i] =si[i]; 				// si at the beginning of the crude interval
  }

  while(ftime<thiscell) {		// start the finer division
    jfine+=1;

    /*initialize output values for those that happen only at event times -  I need them at the end of the crude interval - hence, I set them back to zero everytime I start a new fine interval*/

    dnisi[j] =0;
    dni[j]=0;
    dci[j]=0;
    dnisisq[j]=0;

    timestart = time + ftime;	// time elapsed from the start of the study to the beginning of this interval - this is the time at which the population tables are evaluated


    /*temporary - precision is set to 1!*/
    fthiscell=myprec[0];				//the length of this fine interval is min(precision, time to the end of crude interval)
    //fthiscell=0.1;
    if((thiscell-ftime)<fthiscell) {
      fthiscell=thiscell-ftime;
    }

    sisum=0;
    sisumtt=0;
    /* compute  for each individual within the finer division*/
    for (i=0; i<n; i++) {

      if(y[i]>= times[j]) {				// if still at risk - this is the same throughout the time intervals -  the crude fine intervals are at event and censoring times. Spi must be calculated also for those entering later (period...)
        /*
         ** initialize
         */
        for (k=0; k<edim; k++) {
          data2[k] = x[k][i];						/* the individual's values of demographic variables at time 0 */
        if (efac[k] !=1) data2[k] += timestart;  /* add time to time changing variables */
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

        /* while (etime >0) {*/			//this loop is needed if changes can happen between the interval points.
        et2 = pystep2(edim, &indx, &indx2, &wt, data2, efac, edims, ecut, fthiscell, 1);
        lambdapi = expect[indx];
        lambdapi2 = expect[indx2];
        if(ys[i]<times[j]) {			//he has entered before the crude interval - this guy is at risk for the whole interval - contributes to the values on this interval
          fyidlisi+= lambdapi/si[i];
          fyidlisi2+= lambdapi/(si[i]*exp(-fthiscell* lambdapi));
          fyisi+=1/si[i];
          fyisi2+=1/(si[i]*exp(-fthiscell* lambdapi));
          if (wt <1) hazard+= fthiscell*(wt*lambdapi +(1-wt)*lambdapi2);
          else       hazard+= fthiscell* lambdapi;				//length of the time interval * hazard on this interval
        } // if start of observation before this time

        /*for (k=0; k<edim; k++)
         if (efac[k] !=1) data2[k] += et2;*/
        /*etime -= et2;
        }*/

        si[i] = si[i]* exp(-fthiscell* lambdapi);		//the value of SPi at the end of this fine interval - calculated for all not censored yet, even those not yet at risk

        if(ys[i]<=times[j]) {			//he has entered before the crude interval - this guy is at risk for the whole interval - contributes to the values on this interval
          sisum+=1/si[i];
          sisumtt+=1/sitt[i];

          if(jfine==1) {					//count the number at risk only on the first fine interval
            yi[j]+=1;
          }
        }

        // if this person died at this time
        if(y[i]==times[j]) {
          dnisi[j]+=status[i]/si[i];
          dni[j]+=status[i];
          dci[j]+=(1-status[i]);
          dnisisq[j]+=status[i]/(si[i]*si[i]);
        }


      } // if still at risk
    }// loop through individuals
    fint+= (fyidlisi/fyisi/2 + fyidlisi2/fyisi2/2)*fthiscell;			//the value under the integral at the end of the fine time interval: the product of the value at the beginning * the length of the time interval
    ftime+= fthiscell;
  }// loop through fine times

  yisi[j]=sisum;						//sum of 1/si at the end of the crude interval
  yisitt[j]=sisumtt;					//sum of 1/si at the beginning of the crude interval
  yidlisi[j]=hazard/sisum;				//the total hazard divided by the si at the end of the crude interval
  yidlisitt[j]=hazard/sisumtt;					// the total hazard divided by the si at the beginning of the crude interval
  yidlisiw[j]=fint;						//this is now my best shot at the integrated value on the interval
  yidli[j]=hazard;						//total hazard (yidlambda) on this interval



  time  += thiscell;
}// loop through crude times					AAAA

/*
 ** package the output
 */
PROTECT(rlist = allocVector(VECSXP, 11));					//number of variables
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
SET_VECTOR_ELT(rlist,10, dci2);




PROTECT(rlistnames= allocVector(STRSXP, 11));					//number of variables
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
SET_STRING_ELT(rlistnames, 10, mkChar("dci"));



setAttrib(rlist, R_NamesSymbol, rlistnames);

unprotect(13);					/*number of variables + 2*/
return(rlist);
}
