/*
 **  calculation of various quantities needed for the rs.surv function - sums over individuals at each time
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

/* using thernau's habit: name a S object "charlie2" and the pointer
 **  to the contents of the object "charlie"; the latter is
 **  used in the computations
 */
SEXP netweiDM(   SEXP   efac2,   SEXP edims2,
                 SEXP   ecut2,     SEXP   expect2,
                 SEXP   x2, 	SEXP   y2, 	SEXP ys2, SEXP status2,     SEXP times2) {
  int i,j,k;
  int     n,
  edim,
  ntime;
  double  **x;
  double  *data2, *si, *si2;
  double  **ecut, *etemp;
  double  hazard;						   /*cum hazard over an interval */
double     thiscell,
etime,
time,
et2;
int   indx,
indx2;
double  wt;

int	    *efac, *edims, *status;
double  *expect, *y, *ys, *times;
SEXP      rlist, rlistnames;

/*my declarations*/

SEXP    yidli2, dnisi2,yisi2,yidlisi2,yi2,dni2,sidli2,sidliD2,dnisisq2,yisisq2,sis2,sisD2,yisidli2,yisis2,yidsi2,sit2;
double  *yidli, *dnisi,*yisi,*yidlisi,*yi,*dni,*sidli,*sidliD,*dnisisq,*yisisq,*sis,*sisD,*yisidli,*yisis,*yidsi,*sit;


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

si = (double *)ALLOC(n, sizeof(double));			/*Si for each individual - to je zdaj pointer, vrednosti klicem s s[i]*/
si2 = (double *)ALLOC(n, sizeof(double));			/*Si for each individual - to je zdaj pointer, vrednosti klicem s s[i]*/

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
PROTECT(yidlisi2 = allocVector(REALSXP, ntime));		/*sum yi/Si dLambdai for each time* - length=length(times2)*/
yidlisi = REAL(yidlisi2);
PROTECT(yi2 = allocVector(REALSXP, ntime));					/* sum yi at each time*/
yi = REAL(yi2);
PROTECT(dni2 = allocVector(REALSXP, ntime));		/*sum Yi dLambdai for each time* - length=length(times2)*/
dni = REAL(dni2);
PROTECT(sidli2 = allocVector(REALSXP, ntime));		/*sum dNi/Si for each time* - length=length(times2)*/
sidli = REAL(sidli2);
PROTECT(sidliD2 = allocVector(REALSXP, ntime));		/*sum dNi/Si for each time* - length=length(times2)*/
sidliD = REAL(sidliD2);
PROTECT(yisisq2 = allocVector(REALSXP, ntime));		/*sum Yi/Si for each time* - length=length(times2)*/
yisisq = REAL(yisisq2);
PROTECT(dnisisq2 = allocVector(REALSXP, ntime));		/*sum yi/Si dLambdai for each time* - length=length(times2)*/
dnisisq = REAL(dnisisq2);
PROTECT(sis2 = allocVector(REALSXP, ntime));					/* sum of Si at each time*/
sis = REAL(sis2);
PROTECT(sisD2 = allocVector(REALSXP, ntime));					/* sum of Si at each time*/
sisD = REAL(sisD2);
PROTECT(yisidli2 = allocVector(REALSXP, ntime));					/* sum of Si*dLambdai*Yi at each time*/
yisidli = REAL(yisidli2);
PROTECT(yisis2 = allocVector(REALSXP, ntime));					/* sum of Si*Yi at each time*/
yisis = REAL(yisis2);
PROTECT(sit2 = allocVector(REALSXP, n));					/* Si for each individual*/
sit = REAL(sit2);
PROTECT(yidsi2 = allocVector(REALSXP, ntime));					/* sum of dSi*Yi at each time*/
yidsi = REAL(yidsi2);



/*initialize Si values*/
for (i=0; i<n; i++) {
  si[i] =1;
  si2[i] =1;
  sit[i]=0;
}


/*initialize output values*/
for (j=0; j<ntime; j++) {
  yidli[j] =0;
  dnisi[j] =0;
  yisi[j]=0;
  yidlisi[j]=0;
  yi[j]=0;
  dni[j]=0;
  sidli[j]=0;
  sidliD[j]=0;
  dnisisq[j]=0;
  yisisq[j]=0;
  sis[j]=0;
  sisD[j]=0;
  yisidli[j]=0;
  yisis[j]=0;
  yidsi[j]=0;
}

time =0;
for (j=0; j<ntime ; j++) {			/* loop in time */


thiscell = times[j] - time;

  /* compute  for each individual*/
  for (i=0; i<n; i++) {
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
    while (etime >0) {
      et2 = pystep(edim, &indx, &indx2, &wt, data2, efac,
                   edims, ecut, etime, 1);
      //sit[i]+=1/expect[indx]*(si[i]* exp(-hazard)- si[i]* exp(-hazard + et2*expect[indx]));
      if(expect[indx]==0) expect[indx]=0.000000001;
      if (wt <1) hazard+= et2*(wt*expect[indx] +(1-wt)*expect[indx2]);
      else       hazard+= et2* expect[indx];

      for (k=0; k<edim; k++)
        if (efac[k] !=1) data2[k] += et2;
        etime -= et2;

    }
    sit[i]+=si[i]*(1-exp(-hazard))/(hazard/thiscell);
    si[i] = si[i]* exp(-hazard);
    sis[j]+=si[i];
    sidli[j]+=hazard*si[i];

    si2[i] = si2[i]* exp(-hazard);
    if(y[i]>= times[j]){
      if(ys[i]==times[j]){
        si2[i]=1;
      }
      if(ys[i]<times[j]){
        sisD[j]+=si2[i];
        sidliD[j]+=hazard*si2[i];
      }
    }
    if(y[i]>= times[j]){
      yidsi[j]+=exp(-hazard);
      yidli[j]+=hazard;
      yisidli[j]+=hazard*si[i];
      yi[j]+=1;
      yisi[j]+=1/si[i];
      yisisq[j]+=1/(si[i]*si[i]);
      yisis[j]+=si[i];
      yidlisi[j]+=hazard/si[i];
      if(y[i]==times[j]){
        dnisi[j]+=status[i]/si[i];
        dni[j]+=status[i];
        dnisisq[j]+=status[i]/(si[i]*si[i]);
      }
    }
  }
  time  += thiscell;
}

/*
 ** package the output
 */
PROTECT(rlist = allocVector(VECSXP, 16));
SET_VECTOR_ELT(rlist,0, yidli2);
SET_VECTOR_ELT(rlist,1, yidsi2);
SET_VECTOR_ELT(rlist,2, dnisi2);
SET_VECTOR_ELT(rlist,3, yisi2);
SET_VECTOR_ELT(rlist,4, yidlisi2);
SET_VECTOR_ELT(rlist,5, sidli2);
SET_VECTOR_ELT(rlist,6, yi2);
SET_VECTOR_ELT(rlist,7, dnisisq2);
SET_VECTOR_ELT(rlist,8, yisisq2);
SET_VECTOR_ELT(rlist,9, dni2);
SET_VECTOR_ELT(rlist,10, sis2);
SET_VECTOR_ELT(rlist,11, yisidli2);
SET_VECTOR_ELT(rlist,12, yisis2);
SET_VECTOR_ELT(rlist,13, sit2);
SET_VECTOR_ELT(rlist,14, sidliD2);
SET_VECTOR_ELT(rlist,15, sisD2);

PROTECT(rlistnames= allocVector(STRSXP, 16));
SET_STRING_ELT(rlistnames, 0, mkChar("yidli"));
SET_STRING_ELT(rlistnames, 1, mkChar("yidsi"));
SET_STRING_ELT(rlistnames, 2, mkChar("dnisi"));
SET_STRING_ELT(rlistnames, 3, mkChar("yisi"));
SET_STRING_ELT(rlistnames, 4, mkChar("yidlisi"));
SET_STRING_ELT(rlistnames, 5, mkChar("sidli"));
SET_STRING_ELT(rlistnames, 6, mkChar("yi"));
SET_STRING_ELT(rlistnames, 7, mkChar("dnisisq"));
SET_STRING_ELT(rlistnames, 8, mkChar("yisisq"));
SET_STRING_ELT(rlistnames, 9, mkChar("dni"));
SET_STRING_ELT(rlistnames, 10, mkChar("sis"));
SET_STRING_ELT(rlistnames, 11, mkChar("yisidli"));
SET_STRING_ELT(rlistnames, 12, mkChar("yisis"));
SET_STRING_ELT(rlistnames, 13, mkChar("sit"));
SET_STRING_ELT(rlistnames, 14, mkChar("sidliD"));
SET_STRING_ELT(rlistnames, 15, mkChar("sisD"));


setAttrib(rlist, R_NamesSymbol, rlistnames);

unprotect(18);					/*kolk mora bit tu stevilka??  kolikor jih je +2??*/
return(rlist);
}
