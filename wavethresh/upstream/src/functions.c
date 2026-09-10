#include <R.h> 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>


/* For boundary condition handling */
#define PERIODIC        1
#define SYMMETRIC       2

/* For the type of wavelet decomposition */
#define WAVELET     1   /* The standard decomposition */
#define STATION     2   /* The stationary decomposition */

/* Threshold types */
#define HARD    1
#define SOFT    2

/*
 * ACCESSC handles negative accesses, as well as those that exceed the number
 * of elements
 */

#define ACCESS(image, size, i, j)       *(image + (i)*(size) + (j))
#define ACCESSC(c, firstC, lengthC, ix, bc) *(c+reflect(((ix)-(firstC)),(lengthC),(bc)))
#define ACCESSD(l, i)   *(Data + (*LengthData*(l)) + (i))
#define POINTD(l,i) (Data + (*LengthData*(l)) + (i))
#define POINTC(l,i) (Carray +(*LengthData*(l)) + (i))

/*
 * The next three are exclusively for the stationary wavelet packet algorithm
 * WPST
 */
#define NPKTS(level, nlev)  (1 << (2*(nlev-level)))
#define PKTLENGTH(level)    (1 << level)

#define ACCWPST(a, level, avixstart, pkix, i) *((a) + *(avixstart+(level))+(pkix)*PKTLENGTH(level)+i)

/* Optimiser parameters */

#define R   0.61803399  /* The golden ratio for bisection searches */
#define Cons    (1.0-R)     /* For bisection searches          */

/* These next 3 are for the ipndacw code */
#define ACCESSW(w,j,k)  *(*(w+j)+k)
#define max(a,b)        ((a) > (b) ? (a) : (b))
#define min(a,b)        ((a) > (b) ? (b) : (a))

/*
 * The next 5 are for the swt2d code
 */
 

#define ACCESS3D(ar, d1, d12, ix1, ix2, ix3)    *(ar + (ix3)*(d12)+ (ix2)*(d1)+(ix1))

#define TYPES   0
#define TYPEH   1
#define TYPEV   2
#define TYPED   3

/*
 * End of the swt2d  macro code
 */


/* The code starts here !! */

/*
 * Do wavelet cross-validation in C
 */

void CWaveletCV(double *noisy, int *nnoisy, double *UniversalThresh,
    double *C, double *D, int *LengthD, double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *tol, int *maxits, double *xvthresh,
    int *interptype, int *error)
/* double *noisy;		 The noisy data				*/
/* int *nnoisy;			 Length of noisy data			*/
/* double *UniversalThresh;	 The universal threshold		*/
/* double *C;			 Workspace for scaling coefficients	*/
/* double *D;			 Workspace for wavelet coefficients	*/
/* int *LengthD;		 Length of D workspace			*/
/* double *H;  			 The wavelets to use                 	*/
/* int *LengthH;   		 The length of the filter       	*/
/* int *levels;			 Number of levels			*/
/* int *firstC, *lastC, *offsetC; array indexing info for C vector	*/
/* int *firstD, *lastD, *offsetD; array indexing info for D vector	*/
/* int *ntt;   			 The threshold type			*/ 
/* int *ll;    			 lowest level to threshold; all above too */
/* int *bc;			 The boundary conditions           	*/
/* double *tol; 		 Tol that causes algorithm termination	*/ 
/* int *maxits;			 Max no. of its permitted in optimization */
/* double *xvthresh;		 Returned cross-validatory threshold	*/
/* int *interptype;		 1=noise interpolate, 2=std interpolate */
/* int *error;			 There was an error!                 	*/
{
register int verbose=0;
register int iterations=0;
double ax, bx,cx;
double x0, x1, x2, x3;
/* NOT NEEDED
double fa,fb,fc
 * END */
double f1,f2;
double ssq, tmp;

void Call_Crsswav(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *interptype, int *error);

ax = 0.0;
bx = *UniversalThresh/2.0;
cx = *UniversalThresh;
x0 = ax;
x3 = cx;

if (*error != 0)    {
    verbose=1;
    *error = 0;
    }
else
    verbose=0;

if (verbose)    {
    Rprintf("Entered WaveletCV\n");
    }
if (fabs(cx - bx) > fabs(bx - ax))  {
    x1 = bx;
    x2 = bx + Cons*(cx-bx);
    }
else    {
    x2 = bx;
    x1 = bx - Cons*(bx-ax);
    }

/* NOT NEEDED 

Call_Crsswav(noisy, nnoisy, &ax, C, D, LengthD, H,
    LengthH, levels, firstC, lastC, offsetC, firstD,
    lastD, offsetD, ntt, ll, bc, &ssq, interptype, error);
 
fa = ssq;

if (*error != 0)    {
    *error += 1000;
    return;
    }

Call_Crsswav(noisy, nnoisy, &bx, C, D, LengthD, H, LengthH, levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    ntt, ll, bc, &ssq, interptype, error);
fb = ssq;

if (*error != 0)    {
    *error += 1100;
    return;
    }

Call_Crsswav(noisy, nnoisy, &cx, C, D, LengthD, H, LengthH, levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    ntt, ll, bc, &ssq, interptype, error);
fc = ssq;

if (*error != 0)    {
    *error += 1200;
    return;
    }
 * END OF NOT NEEDED */

Call_Crsswav(noisy, nnoisy, &x1, C, D, LengthD, H, LengthH, levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    ntt, ll, bc, &ssq, interptype, error);
f1 = ssq;

if (*error != 0)    {
    *error += 1300;
    return;
    }

Call_Crsswav(noisy, nnoisy, &x2, C, D, LengthD, H, LengthH, levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    ntt, ll, bc, &ssq, interptype, error);
f2 = ssq;

if (*error != 0)    {
    *error += 1400;
    return;
    }

/*
 *  Next is the MAIN iterative loop.
 *  As well as checking to see if the solution converges, we need to keep
 *  an eye on the maximum number of iterations.
 */

while((fabs(x3-x0) > *tol*(fabs(x1) + fabs(x2))) && iterations++ < *maxits) {
    if (verbose)    {
    Rprintf("x0=%lf, x1=%lf, x2=%lf, x3=%lf\n", x0,x1,x2,x3);
    Rprintf("f1=%lf, f2=%lf\n", f1,f2);
    /*
    fflush(stdout);
    */
        }
    if (f2 < f1)    {
        x0 = x1;
        x1 = x2;
        x2 = R*x1 + Cons*x3;
        f1 = f2;
        Call_Crsswav(noisy, nnoisy, &x2, C, D, LengthD, H,
            LengthH, levels, firstC, lastC, offsetC, firstD, lastD,
            offsetD, ntt, ll, bc, &ssq, interptype, error);
        f2 = ssq;

        if (*error != 0)    {
            *error += 1500;
            return;
            }
    }

    else    {
        x3 = x2;
        x2 = x1;
        x1 = R*x2 + Cons*x0;
        f2 = f1;
        Call_Crsswav(noisy, nnoisy, &x1, C, D, LengthD, H,
            LengthH, levels, firstC, lastC, offsetC, firstD, lastD,
            offsetD, ntt, ll, bc, &ssq, interptype, error);
        f1 = ssq;

        if (*error != 0)    {
            *error += 1600;
            return;
            }
    }
  }

/*
 * Check to see if we've exceeded maximum iterations and return error, if so
 *
 * Also, return a value in tol that indicates how close to the tolerance
 * we are. 
 */


if (iterations >= *maxits)	{
	*error = 1700;
	*tol = fabs(x3-x0)/(fabs(x1)+fabs(x2));
	return;
	}

if (f1 < f2)
    tmp = x1;
else
    tmp = x2;

x1 = tmp/sqrt(1 - log(2.0)/log((double)*nnoisy));

*xvthresh = x1;
return;
}
/*
 * Wrapper to call Crsswav or Crsswav2 depending on the value of interptype
 * This allows one to easily change which type of interpolation one does. 
 *
 */

void Call_Crsswav(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *interptype, int *error)
/*---------------------
 * Argument description
double *noisy::   The noisy data - power of 2              
int *nnoisy::     The number of noisy data elements, must be power of 2
double *value::   The threshold value at which to estimate CV Score    
double *C::
double *D::
int *LengthD::
double *H::   The wavelets to use                  
int *LengthH::    The length of the filter             
int *levels::
int *firstC, *lastC, *offsetC::
int *firstD, *lastD, *offsetD::
int *ntt::    The threshold type                   
int *ll::     The lowest level to threshold; all levels above too  
int *bc::     The boundary conditions              
double *ssq::     The answer!                      
int *interptype::
int *error::  There was an error!                  
---------------------*/
{
void Crsswav(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *error);

void Crsswav2(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *error);

switch(*interptype) {

    case 1:
        Crsswav(noisy, nnoisy, value, C, D, LengthD,
            H, LengthH, levels, firstC, lastC, offsetC, firstD,
            lastD, offsetD, ntt, ll, bc, ssq, error);
        break;
    case 2:
        Crsswav2(noisy, nnoisy, value, C, D, LengthD,
            H, LengthH, levels, firstC, lastC, offsetC, firstD,
            lastD, offsetD, ntt, ll, bc, ssq, error);
        break;

    default:
        *error += 3000;
        break;
    }
return;
}
/*
 * Do rsswav in C
 *
 * This version interpolates the noisy data and compares it to the
 * function values at the index points. (Crsswav2.c does it the other way
 * round - i.e. interpolates the reconstructed function and compares these
 * to the noisy values). Maybe this function is not as good as Crsswav2.c
 * because we smooth the noise before comparison.
 */

void Crsswav(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *error)
/*---------------------
 * Argument description
 *---------------------
double *noisy::   The noisy data - power of 2              
int *nnoisy::     The number of noisy data elements, must be power of 2
double *value::   The threshold value at which to estimate CV Score    
double *C::
double *D::
int *LengthD::
double *H::   The wavelets to use                  
int *LengthH::    The length of the filter             
int *levels::
int *firstC, *lastC, *offsetC::
int *firstD, *lastD, *offsetD::
int *ntt::    The threshold type                   
int *ll::     The lowest level to threshold:: all levels above too  
int *bc::     The boundary conditions              
double *ssq::     The answer!                      
int *error::  There was an error!                  
---------------------*/
{
register int nodd,i;
int type;
int Dlevels;
int *levs_to_do;
int qlevels;
int local_levels;
double *interps;
double ssq1=0.0;
double tmp;

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
void Cthreshold(double *D, int *LengthD, int *firstD, int *lastD, int *offsetD,
	int *Dlevels, int *ntt, double *value, int *levels,
	int *qlevels, int *bc, int *error);

/*
Rprintf("Crsswav\n");
Rprintf("LengthH is %ld\n", *LengthH);
Rprintf("levels is %ld\n", *levels);
Rprintf("ll is %ld\n", *ll); 
fflush(stdout);
*/

/* Get memory for levels to do array */

local_levels = *levels-1; 

qlevels = local_levels - *ll;

if ((levs_to_do = (int *)malloc((unsigned)qlevels*sizeof(int)))==NULL){
    *error = 1;
    return;
    }
else    for(i=0; i < qlevels; ++i)
        *(levs_to_do+i) = *ll+i;

/* Get memory for interps */

nodd = *nnoisy/2; 

if ((interps = (double *)malloc((unsigned)nodd*sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }

type = (int)WAVELET;    /* Only for wavelet transform */

/* Do the ODD analysis and reconstruction */
/*
 * Now copy odd elements to C array 
 */

for(i=0; i<nodd; ++i)   {
    *(C+i) = *(noisy+(2*i));
    /* Rprintf("Copied %lf to C[%d]\n", *(C+i), i); */
    }

/*
 * Do wavelet decomposition of odd, then threshold, then reconstruction.
 */

wavedecomp(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC,
    firstD, lastD, offsetD,
    &type, bc, error); 

/*Rprintf("Done wavelet decomposition\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n"); */

if (*error != 0)    {
    *error += 10;
    return;
    }
/* Threshold */

Dlevels = local_levels - 1;
Cthreshold(D, LengthD, firstD, lastD, offsetD, &Dlevels, ntt, value,
    levs_to_do, &qlevels, bc, error);

if (*error != 0)    {
    *error += 20;
    return;
    }
/*
Rprintf("Done thresholding\n");
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");

*/
/* Reconstruct */

waverecons(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    &type, bc, error);

if (*error != 0)    {
    *error += 30;
    return;
    }
/*
Rprintf("Done reconstruction\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
Rprintf("-----\n");
*/

/* Now work out interpolated evens from noisy */

for(i=0; i<nodd-1; ++i) {
    *(interps+i+1) = (*(noisy+2*i+1) + *(noisy+2*i+3))/2;
    }
*(interps) = *(noisy+1);

/*
for(i=0; i<nodd; ++i)
    Rprintf("interps[%d] is %lf\n", i, *(interps+i));
*/

/* Now form sum of squares */

for(i=0; i<nodd; ++i)   {
    tmp = *(interps+i) - *(C+i);
    ssq1 += tmp*tmp;
    /*Rprintf("ssq1[%d] is %lf\n", i, ssq1); */
    }

*ssq = ssq1;

/* Do the EVEN analysis and reconstruction */
/*
 * Now copy even elements to C array 
 */

for(i=0; i<nodd; ++i)   {
    *(C+i) = *(noisy+(2*i+1));
    /* Rprintf("Copied %lf to C[%d]\n", *(C+i), i); */
    }

/*
 * Do wavelet decomposition of even, then threshold, then reconstruction.
 */

wavedecomp(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC,
    firstD, lastD, offsetD,
    &type, bc, error); 

/*
Rprintf("Done wavelet decomposition\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");
*/
if (*error != 0)    {
    *error += 40;
    return;
    }
/* Threshold */

Dlevels = local_levels - 1;
Cthreshold(D, LengthD, firstD, lastD, offsetD, &Dlevels, ntt, value,
    levs_to_do, &qlevels, bc, error);

if (*error != 0)    {
    *error += 50;
    return;
    }
/*
Rprintf("Done thresholding\n");
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");
*/

/* Reconstruct */

waverecons(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    &type, bc, error);

if (*error != 0)    {
    *error += 60;
    return;
    }

/* Now work out interpolated odds from noisy */

for(i=0; i<nodd-1; ++i) {
    *(interps+i+1) = (*(noisy+2*i) + *(noisy+2*i+2))/2;
    }
*(interps) = *(noisy);

/* Now form sum of squares */

ssq1 = 0.0;

for(i=0; i<nodd; ++i)   {
    tmp = *(interps+i) - *(C+i);
    ssq1 += tmp*tmp;
/*  Rprintf("ssq1[%d] is %lf\n", i, ssq1); */
    }

*ssq = (*ssq + ssq1)/2.0;
free((char *)levs_to_do);
free((char *)interps);
}
/*
 * Do rsswav in C
 *
 * This version interpolates the reconstructed function and compares it to the
 * noisy data at the noisy data index points. (Crsswav.c does it the other way
 * round - i.e. interpolates the noise and compares this to the reconstruced
 * function on the index points of the reconstructed function)
 * Maybe this function is better than Crsswav.c
 * because we don't smooth the noise before comparison.
 */

void Crsswav2(double *noisy, int *nnoisy, double *value,
    double *C, double *D, int *LengthD,
    double *H, int *LengthH,
    int *levels, int *firstC, int *lastC, int *offsetC,
    int *firstD, int *lastD, int *offsetD,
    int *ntt, int *ll, int *bc, double *ssq, int *error)
/*---------------------
 * Argument description
 *---------------------
double *noisy::   The noisy data - power of 2              
int *nnoisy::     The number of noisy data elements, must be power of 2
double *value::   The threshold value at which to estimate CV Score    
double *C::
double *D::
int *LengthD::
double *H::   The wavelets to use                  
int *LengthH::    The length of the filter             
int *levels::
int *firstC, *lastC, *offsetC::
int *firstD, *lastD, *offsetD::
int *ntt::    The threshold type                   
int *ll::     The lowest level to threshold:: all levels above too  
int *bc::     The boundary conditions              
double *ssq::     The answer!                      
int *error::  There was an error!                  
 ---------------------*/
{
register int nodd,i;
int type;
int Dlevels;
int *levs_to_do;
int qlevels;
int local_levels;
double *interps;
double ssq1=0.0;
double tmp;

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
void Cthreshold(double *D, int *LengthD, int *firstD, int *lastD, int *offsetD,
	int *Dlevels, int *ntt, double *value, int *levels,
	int *qlevels, int *bc, int *error);

/*
Rprintf("Crsswav\n");
Rprintf("LengthH is %ld\n", *LengthH);
Rprintf("levels is %ld\n", *levels);
Rprintf("ll is %ld\n", *ll); 
fflush(stdout);
*/

/* Get memory for levels to do array */

local_levels = *levels-1; 

qlevels = local_levels - *ll;

if ((levs_to_do = (int *)malloc((unsigned)qlevels*sizeof(int)))==NULL){
    *error = 1;
    return;
    }
else    for(i=0; i < qlevels; ++i)
        *(levs_to_do+i) = *ll+i;

/* Get memory for interps */

nodd = *nnoisy/2; 

if ((interps = (double *)malloc((unsigned)nodd*sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }

type = (int)WAVELET;    /* Only for wavelet transform */

/* Do the ODD analysis and reconstruction */
/*
 * Now copy odd elements to C array 
 */

for(i=0; i<nodd; ++i)   {
    *(C+i) = *(noisy+(2*i));
    /* Rprintf("Copied %lf to C[%d]\n", *(C+i), i); */
    }

/*
 * Do wavelet decomposition of odd, then threshold, then reconstruction.
 */

wavedecomp(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC,
    firstD, lastD, offsetD,
    &type, bc, error); 

/*Rprintf("Done wavelet decomposition\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n"); */

if (*error != 0)    {
    *error += 10;
    return;
    }
/* Threshold */

Dlevels = local_levels - 1;
Cthreshold(D, LengthD, firstD, lastD, offsetD, &Dlevels, ntt, value,
    levs_to_do, &qlevels, bc, error);

if (*error != 0)    {
    *error += 20;
    return;
    }
/*
Rprintf("Done thresholding\n");
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");

*/
/* Reconstruct */

waverecons(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    &type, bc, error);

if (*error != 0)    {
    *error += 30;
    return;
    }
/*
Rprintf("Done reconstruction\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
Rprintf("-----\n");
*/

/* This bit is different to Crsswav.
 * We interpolate the reconstructed function returned in C
 */

for(i=0; i<nodd-1; ++i) {
    *(interps+i) = (*(C+i+1) + *(C+i))/2.0;
    }
*(interps+nodd-1) = (*C + *(C+nodd-1))/2.0;

/*
for(i=0; i<nodd; ++i)
    Rprintf("interps[%d] is %lf\n", i, *(interps+i));
*/

/* Now form sum of squares */

for(i=0; i<nodd; ++i)   {
    tmp = *(interps+i) - *(noisy+2*i+1);
    ssq1 += tmp*tmp;
    /*Rprintf("ssq1[%d] is %lf\n", i, ssq1); */
    }

*ssq = ssq1;

/* Do the EVEN analysis and reconstruction */
/*
 * Now copy even elements to C array 
 */

for(i=0; i<nodd; ++i)   {
    *(C+i) = *(noisy+(2*i+1));
    /* Rprintf("Copied %lf to C[%d]\n", *(C+i), i); */
    }

/*
 * Do wavelet decomposition of even, then threshold, then reconstruction.
 */

wavedecomp(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC,
    firstD, lastD, offsetD,
    &type, bc, error); 

/*
Rprintf("Done wavelet decomposition\n");
for(i=0; i< *LengthC; ++i)
    Rprintf("C[%d] is %lf\n", i, *(C+i));
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");
*/
if (*error != 0)    {
    *error += 40;
    return;
    }
/* Threshold */

Dlevels = local_levels - 1;
Cthreshold(D, LengthD, firstD, lastD, offsetD, &Dlevels, ntt, value,
    levs_to_do, &qlevels, bc, error);

if (*error != 0)    {
    *error += 50;
    return;
    }
/*
Rprintf("Done thresholding\n");
for(i=0; i< *LengthD; ++i)
    Rprintf("D[%d] is %lf\n", i, *(D+i));
Rprintf("-------\n");
*/

/* Reconstruct */

waverecons(C, D, H, LengthH, &local_levels,
    firstC, lastC, offsetC, firstD, lastD, offsetD,
    &type, bc, error);

if (*error != 0)    {
    *error += 60;
    return;
    }

/* Now work out interpolated odds from noisy */

for(i=0; i<nodd-1; ++i) {
    *(interps+i) = (*(C+i+1) + *(C+i))/2.0;
    }
*(interps+nodd-1) = (*C + *(C+nodd-1))/2.0;

/* Now form sum of squares */

ssq1 = 0.0;

for(i=0; i<nodd; ++i)   {
    tmp = *(interps+i) - *(noisy+2*i);
    ssq1 += tmp*tmp;
/*  Rprintf("ssq1[%d] is %lf\n", i, ssq1); */
    }

*ssq = (*ssq + ssq1)/2.0;
free((char *)levs_to_do);
free((char *)interps);
}
/* 
 * Do manual thresholding in C
 */

#define HardThreshold(coef, threshold) (fabs(coef) > (threshold) ? (coef):(0.0))

void Cthreshold(double *D, int *LengthD, int *firstD, int *lastD, int *offsetD,
	int *Dlevels, int *ntt, double *value, int *levels,
	int *qlevels, int *bc, int *error)
{
register int i,j, local_level;
double cough;
double *din;
int reflect(int n, int lengthC, int bc);
double SoftThreshold(double cough, double threshold);

/*
Rprintf("Cthreshold\n");
Rprintf("LengthD is %ld\n", *LengthD);
Rprintf("ntt is %ld\n", *ntt);
Rprintf("value is %lf\n", *value);
Rprintf("qlevels is %ld\n", *qlevels);
*/

*error = 0;

/*
 * Check that threshold value is positive or zero
 */

if (*value < 0.0)   {
    *error = 3;
    return;
    }

/*
 * Check to see that the levels we are asked to look at are legal
 */

for(i=0; i<*qlevels; ++i)   {
    if (*(levels+i) > *Dlevels) {
        *error = 1;
        return;
        }
    }
/*
 * Now do the thresholding
 */

if (*ntt == HARD)   {
    for(i=0; i<*qlevels; ++i)   {
        local_level =  *(levels+i);

        /*
         * Make din point to correct place in D array
         */

        din = D+*(offsetD+local_level);

        /*
         * Now go through this array doing the thresholding
         */

        for(j= *(firstD+local_level); j<= *(lastD+local_level); ++j){
            cough = ACCESSC(din, (int)*firstD, (int)*LengthD, j,
                (int)*bc);
            cough = HardThreshold(cough, *value);
            ACCESSC(din, (int)*firstD, (int)*LengthD, j, (int)*bc)
                = cough;
            }
        }
    }
else if (*ntt == SOFT)  {
    for(i=0; i<*qlevels; ++i)   {
        local_level =  *(levels+i);

        /*
         * Make din point to correct place in D array
         */

        din = D+*(offsetD+local_level);

        /*
         * Now go through this array doing the thresholding
         */

        for(j= *(firstD+local_level); j<= *(lastD+local_level); ++j){
            cough = ACCESSC(din, (int)*firstD, (int)*LengthD, j,
                (int)*bc);
            cough = SoftThreshold(cough, *value);
            ACCESSC(din, (int)*firstD, (int)*LengthD, j, (int)*bc)
                = cough;
            }
        }
    }
else    {
    *error = 2;
    return;
    }
}

double SoftThreshold(double cough, double threshold)
{
register double s=1.0;

if (cough < 0.0)
    s = -1.0;

if (fabs(cough) > threshold)    {
    return(s*(fabs(cough) - threshold));
    }
else
    return(0.0);
}
/*
 * Function that estimates function with removed observation
 */

void EstWitRem(double *ynoise, int *Lynoise, int *removed, double *thresh,
	double *H, int *LengthH, int *ntt, int *ll, double *answer, int *error)
/*---------------------
 * Argument description
 *---------------------
double *ynoise::      The data                 
int *Lynoise::        The length of the data           
int *removed::        The index to remove from the data        
double *thresh::
double *H::       The wavelets to use                                  
int *LengthH::   The length of the filter                             
int *ntt::       The threshold type                                   
int *ll::        The lowest level to threshold; all levels above too  
double *answer::
int *error::      Possible errors              
 ---------------------*/
{
register int i; /* Register int?                */
int nleft, nright;  /* The number of data points to the left & right */
int nleftExtend;    /* The length of the leftEx vector      */
int nrightExtend;   /* The length of the rightEx vector     */
double *leftEx;     /* Array that contains left and it's extension  */
double *rightEx;    /* Array that contains right and it's extension */
int Dlevels;
int *levs_to_do;
int qlevels;
int local_levels;
int bc;
int type;
double *C, *D;
int LengthC, LengthD;
int levels;
int *firstC, *lastC, *offsetC;
int *firstD, *lastD, *offsetD;

void simpleWT(double *TheData, int *ndata, double *H, int *LengthH,
    double **C, int *LengthC, double **D, int *LengthD, int *levels,
    int **firstC, int **lastC, int **offsetC,
    int **firstD, int **lastD, int **offsetD,
    int *type, int *bc, int *error);

void Cthreshold(double *D, int *LengthD, int *firstD, int *lastD, int *offsetD,
	int *Dlevels, int *ntt, double *value, int *levels,
	int *qlevels, int *bc, int *error);

void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);

int LargerPowerOfTwo(int n);    /* A function that returns next larger  */
                		/* power of two than it's argument  */

/* No errors yet */
*error = 0;

/*
 * Compute number of elements in left and right sections
 */

--(*removed);

nleft = *removed;   /* To the left there will be "removed" elements */

nright = *Lynoise - *removed - 1;   /* To the right there is this   */


/*
 * Now to do the wavelet transform we have to make two vectors that are
 * a power of two in length and that are just inter than twice nleft and nright
 */

nleftExtend = LargerPowerOfTwo(2*nleft);
nrightExtend= LargerPowerOfTwo(2*nright);

/* Now check that we can do the thresholding using the ll number of levels */

/*
 * Get memory for these extensions
 */

if ((leftEx = (double *)malloc((size_t)nleftExtend*sizeof(double)))==NULL){
    *error = 2003;
    return;
    }

if ((rightEx = (double *)malloc((size_t)nrightExtend*sizeof(double)))==NULL){
    *error = 2004;
    return;
    }

/*
 * Now fill these extensions up
 */

for(i=0; i<nleft; ++i)  {
    *(leftEx + (nleftExtend-nleft) + i) = *(ynoise+i);
    *(leftEx + (nleftExtend-nleft)-1- i) = *(ynoise+i);
    }

for(i=0; i< (nleftExtend-2*nleft); ++i) /* Fill up left of left with first el*/
    *(leftEx + i) = *(ynoise+nleft-1);

for(i=0; i<nright; ++i) {
    *(rightEx + i) = *(ynoise+*removed+i+1);
    *(rightEx +2*nright-1-i) = *(ynoise+*removed+i);
    }

for(i=0; i< (nrightExtend-2*nright); ++i)
    *(rightEx+i+2*nright) = *(ynoise+*removed+1);

/*
 * Now we can work out the wavelet transforms of leftExtend and rightExtend
 */

simpleWT(leftEx, &nleftExtend, H, LengthH,
    &C, &LengthC, &D, &LengthD, &levels,
    &firstC, &lastC, &offsetC,
    &firstD, &lastD, &offsetD,
    &type, &bc, error);

if (*error != 0)    {
    return;
    }

local_levels = levels;

qlevels = local_levels - *ll;


if (qlevels>0)  {
    if ((levs_to_do = (int *)malloc((size_t)qlevels*sizeof(int)))==NULL){
        *error = 2005;
        return;
        }
    else    for(i=0; i < qlevels; ++i)
            *(levs_to_do+i) = *ll+i;

    /* Threshold */
    Dlevels = local_levels - 1;
    Cthreshold(D, &LengthD, firstD, lastD, offsetD, &Dlevels, ntt, thresh,
        levs_to_do, &qlevels, &bc, error);

    if (*error != 0)    {
        return;
        }

    free((void *)levs_to_do);
    }

else if (qlevels <0)    {
    *error = 2002;  /* ll, the depth of thresholding exceeded the number
             * of levels that were available for this smaller
             * wavelet transform 
             */
    return;
    }

waverecons(C, D, H, LengthH, &local_levels,
        firstC, lastC, offsetC, firstD, lastD, offsetD,
        &type, &bc, error);

if (*error != 0)        {
        return;
        }

/* Now transfer them back to leftEx */

for(i=0; i<nleftExtend; ++i)
    *(leftEx + i) = *(C+i);

/* Now free the various arrays */

free((void *)C);    free((void *)D);
free((void *)firstC);   free((void *)lastC);    free((void *)offsetC);
free((void *)firstD);   free((void *)lastD);    free((void *)offsetD);

/* Now repeat everything for the right one */

simpleWT(rightEx, &nrightExtend, H, LengthH,
    &C, &LengthC, &D, &LengthD, &levels,
    &firstC, &lastC, &offsetC,
    &firstD, &lastD, &offsetD,
    &type, &bc, error);

if (*error != 0)    {
    return;
    }

local_levels = levels;

qlevels = local_levels - *ll;

if (qlevels>0)  {
    if ((levs_to_do = (int *)malloc((size_t)qlevels*sizeof(int)))==NULL){
        *error = 2005;
        return;
        }
    else    for(i=0; i < qlevels; ++i)
            *(levs_to_do+i) = *ll+i;

    /* Threshold */
    Dlevels = local_levels - 1;
    Cthreshold(D, &LengthD, firstD, lastD, offsetD, &Dlevels, ntt, thresh,
        levs_to_do, &qlevels, &bc, error);

    if (*error != 0)    {
        return;
        }

    free((void *)levs_to_do);
    }
else if (qlevels<0) {
    *error = 2001;
    return; /* ll was too large for this smaller transform */
    }

waverecons(C, D, H, LengthH, &local_levels,
        firstC, lastC, offsetC, firstD, lastD, offsetD,
        &type, &bc, error);

if (*error != 0)        {
        return;
        }

/* Now transfer them back to rightEx */

for(i=0; i<nrightExtend; ++i)
    *(rightEx + i) = *(C+i);

/* Now free the various arrays */

free((void *)C);    free((void *)D);
free((void *)firstC);   free((void *)lastC);    free((void *)offsetC);
free((void *)firstD);   free((void *)lastD);    free((void *)offsetD);

/* Now work out what the removed value should be        */

*answer = *(leftEx + nleftExtend - 1) + *rightEx;
*answer /= 2.0;

/* Now free leftEx and rightEx */
free((void *)leftEx);   free((void *)rightEx);
*error = 0;
return;
}

/*
 * LargerPowerOfTwo:    Returns smallest power of two larger than n
 */

int LargerPowerOfTwo(int n)
{
register int cnt=0;
n--;

while(n>>=1)
    ++cnt;

n = 1;
++cnt;
while(cnt--)
    n<<=1;
return(n);
}
/*
 * Do wavelet cross-validation in C
 */

void FullWaveletCV(double *noisy, int *nnoisy, double *UniversalThresh,
	double *H, int *LengthH, int *ntt, int *ll, double *tol,
	double *xvthresh, int *error)
/*---------------------
 * Argument description
 *---------------------
double *noisy::
int *nnoisy::
double *UniversalThresh::
double *H::   The wavelets to use                  
int *LengthH::    The length of the filter             
int *ntt::    The threshold type                   
int *ll::     The lowest level to threshold; all levels above too  
double *tol::
double *xvthresh::
int *error::  There was an error!                  
 ---------------------*/
{
int verbose=0;
double ax, bx,cx;
double x0, x1, x2, x3;
/* NOT NEEDED
double fa,fb,fc; */
double f1,f2;
double ssq;

int mRi;    /* This is required as an argument to GetRSS, but we don't
         * make use of it here */

void GetRSS(double *ynoise, int *Lynoise, double *thresh, double *H,
	int *LengthH, int *ntt, int *ll, double *rss, int *smallestRSSindex,
	int *verbose, int *error);

ax = 0.0;
bx = *UniversalThresh/2.0;
cx = *UniversalThresh;
x0 = ax;
x3 = cx;

if (*error != 0)    {
    verbose=1;
    *error = 0;
    }
else
    verbose=0;

if (verbose)    {
    Rprintf("Entered FullWaveletCV\n");
    }
if (fabs(cx - bx) > fabs(bx - ax))  {
    x1 = bx;
    x2 = bx + Cons*(cx-bx);
    }
else    {
    x2 = bx;
    x1 = bx - Cons*(bx-ax);
    }


if (verbose)    {
    Rprintf("About to enter GetRSS for the first time\n");
    }

/* NOT NEEDED 

GetRSS(noisy, nnoisy, &ax, H, LengthH, ntt, ll, &ssq, &mRi, &verbose, error);
if (verbose)    {
    Rprintf("Just left GetRSS for the first time\n");
    }
fa = ssq;

if (*error != 0)    {
    *error += 1000;
    return;
    }

GetRSS(noisy, nnoisy, &bx, H, LengthH, ntt, ll, &ssq, &mRi, &verbose, error);
fb = ssq;

if (*error != 0)    {
    *error += 1100;
    return;
    }

GetRSS(noisy, nnoisy, &cx, H, LengthH, ntt, ll, &ssq, &mRi, &verbose, error);
fc = ssq;

if (*error != 0)    {
    *error += 1200;
    return;
    }

 * END OF NOT NEEDED */

GetRSS(noisy, nnoisy, &x1, H, LengthH, ntt, ll, &ssq, &mRi, &verbose, error);
f1 = ssq;

if (*error != 0)    {
    *error += 1300;
    return;
    }

GetRSS(noisy, nnoisy, &x2, H, LengthH, ntt, ll, &ssq, &mRi, &verbose, error);
f2 = ssq;

if (*error != 0)    {
    *error += 1400;
    return;
    }

while(fabs(x3-x0) > *tol*(fabs(x1) + fabs(x2))) {
    if (verbose)    {
    Rprintf("x0=%lf, x1=%lf, x2=%lf, x3=%lf\n", x0,x1,x2,x3);
    Rprintf("f1=%lf, f2=%lf\n", f1,f2);
    /*
    fflush(stdout);
    */
        }
    if (f2 < f1)    {
        x0 = x1;
        x1 = x2;
        x2 = R*x1 + Cons*x3;
        f1 = f2;
        GetRSS(noisy, nnoisy, &x2, H, LengthH, ntt, ll, &ssq, &mRi,
            &verbose, error);
        f2 = ssq;

        if (*error != 0)    {
            *error += 1500;
            return;
            }
    }

    else    {
        x3 = x2;
        x2 = x1;
        x1 = R*x2 + Cons*x0;
        f2 = f1;
        GetRSS(noisy, nnoisy, &x1, H, LengthH, ntt, ll, &ssq, &mRi,
            &verbose, error);
        f1 = ssq;

        if (*error != 0)    {
            *error += 1600;
            return;
            }
    }
  }

if (f1 < f2)
    *xvthresh = x1;
else
    *xvthresh = x2;

return;
}
#define GRSTART 4       /* The first index to remove for GetRSS */

void GetRSS(double *ynoise, int *Lynoise, double *thresh, double *H,
	int *LengthH, int *ntt, int *ll, double *rss, int *smallestRSSindex,
	int *verbose, int *error)
{
int removed, local_removed;
int minRSSix;
double TheMinRSS;
double answer;
void EstWitRem(double *ynoise, int *Lynoise, int *removed, double *thresh,
	double *H, int *LengthH, int *ntt, int *ll, double *answer, int *error);

/* No error yet!! */
*error = 0;

*rss = 0.0;
minRSSix = 0;
TheMinRSS = 0.0;

for(removed=GRSTART; removed<= *Lynoise-GRSTART+1; ++removed)   {


    local_removed = removed;
    EstWitRem(ynoise, Lynoise, &local_removed, thresh, H, LengthH,
        ntt, ll, &answer, error);

    if (*error != 0)    {
        return;
        }
    answer -= *(ynoise+removed-1);

    *rss += answer*answer;

    if (removed==GRSTART)   {
        minRSSix = removed;
        TheMinRSS = answer* answer;
        }
    else if (TheMinRSS > answer*answer) {
        minRSSix = removed; 
        TheMinRSS = answer*answer;
        }

    if (*verbose>1)
      Rprintf("GetRSS: Removed is %d, ynoise[%d] is %lf RSS is %lf\n",     /*MAN: changed %ld to %d since declared as int (L1329) */
        removed, removed, *(ynoise+removed-1), *rss);

    }
*rss /= (*Lynoise - 4);
*smallestRSSindex = minRSSix;
return;
}
/*
 * ImageDecomposeStep   -   Take an image and do a one level decomp
 *
 * Error Codes
 *
 *  0   -   Ok.
 *
 *  1   -   Memory error for (afterC) temporary image
 *
 *  2   -   Memory error for (afterD) temporary image
 *
 *  3   -   Memory error for (ccopy) temporary row store
 *
 *  4   -   Memory error for (ccopy_out) temporary row store
 *
 *  5   -   Memory error for (dcopy_out) temporary row store
 *
 *  6-9 -   Memory errors for (afterCC,afterCD,afterDC,afterDD)
 *          store for the answers
 */



void ImageDecomposeStep(double *C, int Csize, int firstCin, double *H,
	int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
	double **cc_out, double **cd_out, double **dc_out, double **dd_out,
	int bc, int type, int *error)
/*---------------------
 * Argument description
 *---------------------
double *C::   Input data image                 
int Csize::  Size of image (side length)              
int firstCin::   Index number of first element in input "C" image 
double *H::   Filter coefficients                  
int LengthH::    Length of filter                 
 Details about output image 
int LengthCout:: Length of C part of output image         
int firstCout::  Index number of first element in output "C" image    
int lastCout::   Index number of last element             
int LengthDout:: Length of D part of output image         
int firstDout::  Index number of first element in output "D" image    
int lastDout::   Index number of last element             
double **cc_out:: Smoothed output image                
double **cd_out:: Horizontal detail                    
double **dc_out:: Vertical detail                  
double **dd_out:: Diagonal detail                  
int bc::     Method of boundary correction            
int type::   Type of transform, wavelet or stationary     
int *error::     Error code                       
 *---------------------*/
{
register int j,row,col;
double *ccopy;  /* Used to copy input data to convolution routines  */
double *ccopy_out;/* Used to copy output data to afterC after conv. */
double *dcopy_out;/* Used to copy output data to afterD after conv. */
double *afterC; /* Temporary store for image data after C convolution   */
double *afterD; /* Temporary store for image data after D convolution   */
double *afterCC,*afterCD,*afterDC,*afterDD; /* Results      */
int step_factor;    /* This should always be 1 for the WAVELET trans*/

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

*error = 0;

step_factor = 1;

/* Get memory for afterC */

if ((afterC = (double *)malloc((unsigned)(Csize*LengthCout*sizeof(double))))==NULL){
    *error = 1;
    return;
    }

/* Get memory for afterD */

if ((afterD = (double *)malloc((unsigned)(Csize*LengthDout*sizeof(double))))==NULL){
        *error = 2;
        return;
        }

/* Get memory for row of image to pass to convolution routines */

if ((ccopy = (double *)malloc((unsigned)(Csize*sizeof(double)))) == NULL) {
    *error = 3;
    return;
    }

/* Get memory for output row after C convolution */

if ((ccopy_out = (double *)malloc((unsigned)(LengthCout*sizeof(double))))==NULL) {
    *error = 4;
    return;
    }

/* Get memory for output row after D convolution */

if ((dcopy_out = (double *)malloc((unsigned)(LengthDout*sizeof(double))))==NULL) {
        *error = 5;
        return;
        }



/* Do convolutions on rows of C */

for(row=0; row < (int)Csize; ++row) {

    /* Copy row of C into ccopy */

    for(j=0; j<Csize; ++j)
        *(ccopy+j) = ACCESS(C, Csize, row, j);

    /* Now convolve this row with C filter */

    convolveC(ccopy, (int)Csize, (int)firstCin,H, (int)LengthH, ccopy_out,
        (int)firstCout, (int)lastCout,
        (int)type, step_factor, (int)bc); 

    /* Now convolve this row with D filter */

    convolveD(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH, dcopy_out,
        (int)firstDout, (int)lastDout,
        (int)type, step_factor, (int)bc);

    /* Copy answer back to arrays */

    for(j=0; j<(int)LengthCout; ++j)
        ACCESS(afterC, (int)LengthCout, row, j) = *(ccopy_out + j);

    for(j=0; j<(int)LengthDout; ++j)
        ACCESS(afterD, (int)LengthDout, row, j) = *(dcopy_out + j);

    }


/* Now we have to apply both the C and D filters to afterC and afterD.
 * We get four answers. First we get the necessary memory
 */

if ((afterCC = (double *)malloc((unsigned)(LengthCout*LengthCout*sizeof(double))
        ))==NULL)   {
    *error = 6;
    return;
    }

if ((afterCD = (double *)malloc((unsigned)(LengthDout*LengthCout*sizeof(double))
        ))==NULL)   {
    *error = 7;
    return;
    }

if ((afterDC = (double *)malloc((unsigned)(LengthCout*LengthDout*sizeof(double))
        ))==NULL)   {
    *error = 8;
    return;
    }

if ((afterDD = (double *)malloc((unsigned)(LengthDout*LengthDout*sizeof(double))
        ))==NULL)   {
    *error = 9;
    return;
    }


/* Link this memory to the returning pointers */

*cc_out = afterCC;
*cd_out = afterCD;
*dc_out = afterDC;
*dd_out = afterDD;


/* Apply the filters, first to afterC to get afterCC and afterCD */

for(col=0; col < (int)LengthCout; ++col)    {

    /* Copy column to ccopy */

    for(j=0; j<(int)Csize; ++j)
        *(ccopy + j) = ACCESS(afterC, (int)LengthCout, j, col);

    /* Apply C filter */

    convolveC(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH, ccopy_out,
        (int)firstCout, (int)lastCout,
        (int)type, step_factor, (int)bc);

    /* Apply D filter */

    convolveD(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH, dcopy_out,
        (int)firstDout, (int)lastDout,
        (int)type, step_factor, (int)bc);

    /* Copy answer back */

    for(j=0; j<(int)LengthCout; ++j)
        ACCESS(afterCC, (int)LengthCout, j, col) = *(ccopy_out+j);

    for(j=0; j<(int)LengthDout; ++j)
        ACCESS(afterCD, (int)LengthCout, j, col) = *(dcopy_out+j);
    }

/* Apply the filters, now to afterD to get afterDC and afterDD */

for(col=0; col < (int)LengthDout; ++col)    {

    /* Copy column to ccopy */

    for(j=0; j<(int)Csize; ++j)
        *(ccopy + j) = ACCESS(afterD, (int)LengthDout, j, col);

    /* Apply C filter */

    convolveC(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH, ccopy_out,
        (int)firstCout, (int)lastCout,
        (int)type, step_factor, (int)bc);

    /* Apply D filter */

    convolveD(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH, dcopy_out,
        (int)firstDout, (int)lastDout,
        (int)type, step_factor, (int)bc);

    /* Copy answer back */

    for(j=0; j<(int)LengthCout; ++j)
        ACCESS(afterDC, (int)LengthDout, j, col) = *(ccopy_out+j);

    for(j=0; j<(int)LengthDout; ++j)
        ACCESS(afterDD, (int)LengthDout, j, col) = *(dcopy_out+j);
    }

/* That should be it ! */

free((char *) afterD);
free((char *) afterC);
free((char *)dcopy_out);
free((char *)ccopy_out);
free((char *)ccopy);

return;
}


void StoIDS(double *C, int *Csize, int *firstCin, double *H, int *LengthH,
	int *LengthCout, int *firstCout, int *lastCout,
	int *LengthDout, int *firstDout, int *lastDout,
	double *ImCC, double *ImCD, double *ImDC, double *ImDD,
	int *bc, int *type, int *error)
{
register int i,j;
double *cc_out, *cd_out, *dc_out, *dd_out;

void ImageDecomposeStep(double *C, int Csize, int firstCin, double *H,
	int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
	double **cc_out, double **cd_out, double **dc_out, double **dd_out,
	int bc, int type, int *error);

ImageDecomposeStep(C, *Csize, *firstCin, H, *LengthH,
    *LengthCout, *firstCout, *lastCout,
    *LengthDout, *firstDout, *lastDout,
    &cc_out, &cd_out, &dc_out, &dd_out, *bc, *type,
    error);

/* Copy images */

for(i=0; i<(int)*LengthDout; ++i)   {
    for(j=0; j<(int)*LengthDout; ++j)
        ACCESS(ImDD, (int)*LengthDout, i, j) = ACCESS(dd_out,
            *LengthDout, i, j);

    for(j=0; j<(int)*LengthCout; ++j)
        ACCESS(ImDC, (int)*LengthDout, j,i) = ACCESS(dc_out,
            *LengthDout, j,i);
    }

for(i=0; i<(int)*LengthCout; ++i)   {
    for(j=0; j<(int)*LengthDout; ++j)
        ACCESS(ImCD, (int)*LengthCout, j,i) = ACCESS(cd_out,
            *LengthCout, j,i);

    for(j=0; j<(int)*LengthCout; ++j)
        ACCESS(ImCC, (int)*LengthCout,j,i) = ACCESS(cc_out,
            *LengthCout, j,i);
    }

free((void *)cc_out);
free((void *)cd_out);
free((void *)dc_out);
free((void *)dd_out);
}

void StoIRS(double *ImCC, double *ImCD, double *ImDC, double *ImDD,
	int *LengthCin, int *firstCin,
	int *LengthDin, int *firstDin,
	double *H, int *LengthH,
	int *LengthCout, int *firstCout, int *lastCout,
	double *ImOut, int *bc, int *error)
{
void ImageReconstructStep(double *ImCC, double *ImCD, double *ImDC,
	double *ImDD,
	int LengthCin, int firstCin,
	int LengthDin, int firstDin,
	double *H, int LengthH,
	int LengthCout, int firstCout, int lastCout,
	double *ImOut, int *bc, int *error);

*error = 0;

ImageReconstructStep(ImCC, ImCD, ImDC, ImDD,
    *LengthCin, *firstCin, 
    *LengthDin, *firstDin,
    H, *LengthH,
    *LengthCout, *firstCout, *lastCout,
    ImOut, bc,
    error);

}

void ImageReconstructStep(double *ImCC, double *ImCD, double *ImDC,
	double *ImDD,
	int LengthCin, int firstCin,
	int LengthDin, int firstDin,
	double *H, int LengthH,
	int LengthCout, int firstCout, int lastCout,
	double *ImOut, int *bc, int *error)
{
register int i,j;
double *c_in;
double *d_in;
double *c_out;
double *toC;
double *toD;
int type=WAVELET;       /* The type of the WAVELET transform    */

void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);

/* Get memory for c_in and d_in */

if ((c_in = (double *)malloc((unsigned)LengthCin*sizeof(double)))==NULL) {
        *error = 1;
        return;
        }

if ((d_in = (double *)malloc((unsigned)LengthDin*sizeof(double)))==NULL) {
        *error = 2;
        return;
        }

if ((c_out = (double *)malloc((unsigned)LengthCout*sizeof(double)))==NULL) {
        *error = 3;
        return;
        }

if ((toC = (double *)malloc((unsigned)LengthCin*LengthCout*sizeof(double)))
    == NULL)    {
        *error = 4;
        return;
        }
        
/* Now apply C and D filters to CC and CD to obtain toC */

for(i=0; i<LengthCin; ++i)  {

    for(j=0; j < LengthDin; ++j)
        *(d_in + j) = ACCESS(ImCD, (int)LengthCin, j, i);

    for(j=0; j < LengthCin; ++j)
        *(c_in + j) = ACCESS(ImCC, (int)LengthCin, j, i);

    conbar(c_in, (int)LengthCin, (int)firstCin,
           d_in, (int)LengthDin, (int)firstDin,
           H, (int)LengthH,
           c_out, (int)LengthCout, (int)firstCout, (int)lastCout,
        type, (int)*bc);

    for(j=0; j < LengthCout; ++j)
        ACCESS(toC, (int)LengthCout, i, j) = *(c_out+j); 
    }

/* Now magically we can use c_in and d_in again, but we now need a toD
 * to store the answer in */

if ((toD = (double *)malloc((unsigned)LengthDin*LengthCout*sizeof(double)))== NULL)    {
    *error = 5;
    return;
    }

/* Now apply C and D filters to DC and DD to obtain toD */

for(i=0; i<LengthDin; ++i)      {

        for(j=0; j < LengthDin; ++j)
                *(d_in + j) = ACCESS(ImDD, (int)LengthDin, j, i);

        for(j=0; j < LengthCin; ++j)
                *(c_in + j) = ACCESS(ImDC, (int)LengthDin, j, i);

        conbar(c_in, (int)LengthCin, (int)firstCin,
               d_in, (int)LengthDin, (int)firstDin,
               H, (int)LengthH,
               c_out, (int)LengthCout, (int)firstCout, (int)lastCout,
        type, (int)*bc);

        for(j=0; j < LengthCout; ++j)
                ACCESS(toD, (int)LengthCout, i, j) = *(c_out+j);
        }


/* Now apply C and D filters to toC and toD to get ImOut */

for(i=0; i<LengthCout; ++i) {

    for(j=0; j< LengthDin; ++j)
        *(d_in + j) = ACCESS(toD, (int)LengthCout, j, i);

    for(j=0; j<LengthCin; ++j)
        *(c_in + j) = ACCESS(toC, (int)LengthCout, j, i);

    conbar(c_in, (int)LengthCin, (int)firstCin,
               d_in, (int)LengthDin, (int)firstDin,
               H, (int)LengthH,
               c_out, (int)LengthCout, (int)firstCout, (int)lastCout,
        type, (int)*bc);

    for(j=0; j<LengthCout; ++j)
        ACCESS(ImOut, (int)LengthCout, i, j)  = *(c_out+j);
    }
}
/* AV_BASIS Do the basis averaging */

/*
 * Error codes
 *
 * 1    -   Memory error in creating cl
 * 2    -   Memory error in creating cr
 * 3    -   Memory error in creating packet (getpacket)
 */


double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wst::	The stationary wavelet decomposition         
double *wstC::  The stationary wavelet decomposition         
int nlevels::  	The original length of the data         
int level::     The level to reconstruct             
int ix1::       The "left" packet index              
int ix2::       The "right" packet index             
double *H::     The filter                       
int LengthH::   The length of the filter             
int *error::    Error code                       
 *---------------------*/
{
register int i;
double *cl;
double *cr;
double *genericC;
double *genericD;
void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);
double *getpacket(double *wst, int nlevels, int level, int index, int *error);
double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error);
int LengthC;
int LengthCin;

void rotateback(double *book, int length);

*error = 0;

/*
 * Now we must create cl and cr. These will contain the reconstructions
 * from the left and right packets respectively. The length of these
 * vectors depends upon the level we're at.
 */

LengthC = 1 << (level+1);
LengthCin = 1 << level;

/*
 * Create cl and cr
 */

if ((cl = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 1;
    return(NULL);
    }

if ((cr = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 2;
    return(NULL);
    }
    
/*
 * What we do next depends on the level.
 *
 * If level is zero then we've recursed all the way down to the bottom of
 * the tree. And we can reconstruct the 2-vectors one-up-the-tree by using
 * good old conbar().
 *
 * If the level is not zero then we construct at that stage using conbar()
 * but to obtain the Cs we recurse. 
 */

if (level != 0) {

    /* Get C's at this level by asking the next level down. */

    genericC = av_basis(wst, wstC, nlevels, level-1, 2*ix1, 2*ix1+1,
            H, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genericD = getpacket(wst, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0, 
           H, LengthH,
           cl, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genericC);
    free((void *)genericD);

    /* Now do the RHS */
    
    genericC = av_basis(wst, wstC, nlevels, level-1, 2*ix2, 2*ix2+1,
        H, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genericD = getpacket(wst, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0,
           H, LengthH,
           cr, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(cr, LengthC);

    /* Can get rid of generics now */

    free((void *)genericC);
    free((void *)genericD);
    }

else    {
    /* Have to really do it! */

    genericC = getpacket(wstC, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    genericD = getpacket(wst, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0, 
           genericD, LengthCin, 0, 
           H, LengthH,
           cl, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genericC);
    free((void *)genericD);

    genericC = getpacket(wstC, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    genericD = getpacket(wst, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    conbar(genericC, LengthCin, 0,
           genericD, LengthCin, 0,
           H, LengthH,
           cr, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(cr, LengthC);

    free((void *)genericC);
    free((void *)genericD);
    }

for(i=0; i<LengthC; ++i)
    *(cl+i) = ((double)0.5)*( *(cl+i) + *(cr+i) );

/*
 *  Return the answer in cl (which has to be freed later)
 *  Destroy pointer to cr, as it is not needed now
 *
 */

free((void *)cr);
return(cl);
}


/* Note this is different to what is in wavepackst.c */
#define AVBPOINTD(w, l,i) (w + (nlevels*(i)) + (l))

/*
 * A C version of getpacket
 *
 * Warning. The argument list is the same as for the WaveThresh Splus version
 * except that nlevels here should be one more!
 */

double *getpacket(double *wst, int nlevels, int level, int index, int *error)
/* int nlevels::    This looks like it should be nlevels+1 for some reason */
{
register int i;
double *packet;
int PacketLength;

PacketLength = 1 << level;

if ((packet = (double *)malloc((unsigned)PacketLength*sizeof(double)))==NULL){
    *error = 3;
    return(NULL);
    }


for(i=0; i< PacketLength; ++i)
    *(packet+i) = *AVBPOINTD(wst, level, (index*PacketLength+i));

return(packet);
}

/* Wrapper for av_basis */


void av_basisWRAP(double *wst, double *wstC, int *LengthData, int *level,
	double *H, int *LengthH, double *answer, int *error)
{
register int i;
int nlevels;
double *acopy;
double *av_basis(double *wst, double *wstC, int nlevels, int level,
	int ix1, int ix2, double *H, int LengthH, int *error);

nlevels = 2 + (int)*level;

acopy =  av_basis(wst, wstC, nlevels, (int)*level, 0, 1, H,
        (int)*LengthH, error);

for(i=0; i< (int)*LengthData; ++i)
    *(answer+i) = *(acopy+i);

free((void *)acopy);
}
/*
 * CONBAR: Does the reconstruction convolution
 */

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc)
{
register int n,k;
register int cfactor;
double sumC, sumD;

int reflect(int n, int lengthC, int bc);

switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumC = 0.0;

    while( cfactor*k <= n ) {

        sumC += *(H + n - cfactor*k)*ACCESSC(c_in, firstCin, LengthCin,
                    k, bc);

        ++k;
        }

    /* Now do D part */

    k = CEIL(n-1);

    sumD = 0.0;

    while( cfactor*k <= (LengthH +n -2) )   {

        sumD += *(H+1+cfactor*k-n) * ACCESSC(d_in, firstDin, LengthDin,
                    k, bc);

        ++k;

        }

    if (n & 1)      /* n odd */
        sumC -= sumD;
    else
        sumC += sumD;

    ACCESSC(c_out, firstCout, LengthCout, n, bc) = sumC;
    }

}
/*
 * CONBARL: Wrapper called by SPlus conbar() to call C conbar.
 */

void conbarL(double *c_in, int *LengthCin, int *firstCin,
	double *d_in, int *LengthDin, int *firstDin,
	double *H, int *LengthH,
	double *c_out, int *LengthCout, int *firstCout, int *lastCout,
	int *type, int *bc)
{
int LLengthCin;
int LfirstCin;
int LLengthDin;
int LfirstDin;
int LLengthH;
int LLengthCout;
int LfirstCout;
int LlastCout;
int Ltype;
int Lbc;
void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);

LLengthCin = (int)*LengthCin;
LfirstCin = (int)*firstCin;
LLengthDin = (int)*LengthDin;
LfirstDin = (int)*firstDin;
LLengthH = (int)*LengthH;
LLengthCout = (int)*LengthCout;
LfirstCout = (int)*firstCout;
LlastCout = (int)*lastCout;
Ltype = (int)*type;
Lbc = (int)*bc;


conbar(c_in, LLengthCin, LfirstCin,
       d_in, LLengthDin, LfirstDin,
       H, LLengthH,
       c_out, LLengthCout, LfirstCout, LlastCout, Ltype, Lbc);
}
/*
 * CONVOLVE -   Do filter H filter convolution with boundary
 */

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_in::    	Input data                       
int LengthCin::   	Length of this array                 
int firstCin::    	The first C value                    
double *H::   		Filter                       
int LengthH::     	Length of filter                 
double *c_out::   	Output data                      
int firstCout::   	First index of C array               
int lastCout::    	Last index of C array                
int type::		Type of wavelet decomposition            
int step_factor::	For stationary wavelets only             
int bc::      		Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sum;
register int k;
register int count_out;
register int m;
register int cfactor;   /* This determines what sort of dilation we do  */
            /* and depends on the type argument     */
int reflect(int n, int lengthC, int bc);

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstCout; k<=lastCout; ++k)  {
    sum = 0.0;


    for(m=0; m<LengthH; ++m)    {


        sum += *(H+m) * ACCESSC(c_in, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        }

    *(c_out + count_out) = sum;
    ++count_out;
    }
}
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_in::    	Input data                       
int LengthCin::   	Length of this array                 
int firstCin::    	The first C value                    
double *H::   		Filter                       
int LengthH::     	Length of filter                 
double *d_out::   	Output data                      
int firstDout::   	First index of C array               
int lastDout::    	Last index of C array                
int type::		Type of wavelet decomposition            
int step_factor::	For stationary wavelets only             
int bc::      		Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sum;
double tmp;
register int k;
register int count_out;
register int m;
register int cfactor;

int reflect(int n, int lengthC, int bc);

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstDout; k<=lastDout; ++k)  {
    sum = 0.0;


    for(m=0; m<LengthH; ++m)    {

        tmp = ACCESSC(c_in, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);
        
        if (m&1)    /* odd */
            sum += *(H+m) *  tmp;
        else
            sum -= *(H+m) *  tmp;
        
        }

    *(d_out + count_out) = sum;
    ++count_out;
    }
}


/* Works out reflection, as REFLECT, but reports access errors */
int reflect(int n, int lengthC, int bc)
{

if ((n >= 0) && (n < lengthC))
    return(n);
else if (n<0)   {
    if (bc==PERIODIC)   {
        /*
        n = lengthC+n;
        */
        n = n%lengthC + lengthC*((n%lengthC)!=0);
        if (n < 0)      {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: left info from right\n");
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }

    else if (bc==SYMMETRIC) {
        n = -1-n;
        if (n >= lengthC)       {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }

    else    {
        REprintf("reflect: Unknown boundary correction");
        REprintf("value of %d\n", bc);
        error("This should not happen. Stopping.\n");
        }

    }
else    {
    if (bc==PERIODIC)   {
        /*
        Rprintf("periodic extension, was %d (%d) now ",n,lengthC);
        n = n - lengthC; 
        */
        n %= lengthC;
        /*
        Rprintf("%d\n", n);
        */
        if (n >= lengthC)   {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
            REprintf("reflect: right info from left\n");
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }
    else if (bc==SYMMETRIC) {
        n = 2*lengthC - n - 1;
        if (n<0)        {
            REprintf("reflect: access error (%d,%d)\n",
                n,lengthC);
	    error("This should not happen. Stopping.\n");
            }
        else
            return(n);
        }
    else    {
        REprintf("reflect: Unknown boundary correction\n");
	error("This should not happen. Stopping.\n");
        }


    }
/* Safety */
REprintf("reflect: SHOULD NOT HAVE REACHED THIS POINT\n");
error("This should not happen. Stopping.\n");
return(0); /* for lint only */
}

/* Rotate a vector */

/* Vector: a_1, a_2, a_3, ..., a_{n-1}, a_n

    becomes

       a_2, a_3, a_4, ..., a_n, a_1

   rotateback() does the opposite

*/

void rotater(double *book, int length)
{
register int i;
double tmp;

tmp = *book;

for(i=0; i<length-1; ++i)
        *(book+i) = *(book+i+1);

*(book+length-1) = tmp;
}

void rotateback(double *book, int length)
{
register int i;
double tmp;

tmp = *(book+length-1);

for(i= length-1; i>0; --i)
    *(book+i) = *(book+i-1);

*book = tmp;
}

/*
 * Does a simple wavelet transform
 *
 * This is just like the ordinary periodic wavelet transform
 *
 * The purpose of this function is simplicity. All you need supply is the
 * data and some pointers for the arguments
 *
 * This function uses Calloc to create the arrays:
 *
 *   firstC,lastC,offsetC,firstD,lastD,offsetD,C,D
 *
 * When you have used their contents it is a good idea for you to destroy
 * the memory associated with these arrays. To do this call
 *
 * free((char *)C);     - This frees the memory associated with the pointer
 *
 * [We recommend you declare these arrays in the calling program like
 *
 *  double *C, *D;
 *  int *firstC, *lastC, *offsetC, *firstD, *lastD, *offsetD;
 *
 * Then PASS the ADDRESS of these to this function, e.g. 
 * &C, &D etc.]
 */

void simpleWT(double *TheData, int *ndata, double *H, int *LengthH,
    double **C, int *LengthC, double **D, int *LengthD, int *levels,
    int **firstC, int **lastC, int **offsetC,
    int **firstD, int **lastD, int **offsetD,
    int *type, int *bc, int *error) 
/*---------------------
 * Argument description
 *---------------------
double *TheData::	The data to transform; must be a power of two els     
int *ndata::  		The length of the data               
double *H::   		The wavelet filter that you want to use      
int *LengthH::    	The length of the wavelet filter         
 The following arguments are the answer/output
double **C::  		A pointer to the array of C answers is returned  
int *LengthC::    	The length of the C array is returned        
double **D::  		A pointer to the array of D answers is returned  
int *LengthD::    	The length of the D array is returned        
int *levels::     	The number of levels of the transform is returned    
int **firstC,**lastC,**offsetC:: These are computed and returned   
int **firstD,**lastD,**offsetD:: These are computed and returned   
int *type::   		This is filled in with type WAVELET          
int *bc::     		This is filled in with PERIODIC          
int *error::  		Returns any error condition              
 *---------------------*/
{
int *lfC,*llC,*loC; /* Local versions of firstC,lastC,offsetC   */
int *lfD,*llD,*loD; /* Local versions of firstD,lastD,offsetD   */
double *lC, *lD;    /* Local versions of C and D            */
int cnt,i;

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
int IsPowerOfTwo(int n);

/* No errors yet */

*error = 0;

/* Fill in type of transform and type of boundary handling conditions   */

*type = (int)WAVELET;
*bc = (int)PERIODIC;

/* Now work out the size of the arrays needed for the transform     */

*levels = (int)IsPowerOfTwo(*ndata);

/* Now create memory for first/last and offset */

/* Now create memory for first/last and offset */

if ((lfC = (int *)malloc((size_t)(*levels+1)*sizeof(int)))==NULL) {
    *error = 3001;
    return;
    }

if ((llC = (int *)malloc((size_t)(*levels+1)*sizeof(int)))==NULL) {
    *error = 3002;
    return;
    }

if ((loC = (int *)malloc((size_t)(*levels+1)*sizeof(int)))==NULL) {
    *error = 3003;
    return;
    }

if ((lfD = (int *)malloc((size_t)(*levels)*sizeof(int)))==NULL)   {
    *error = 3004;
    return;
    }

if ((llD = (int *)malloc((size_t)(*levels)*sizeof(int)))==NULL)   {
    *error = 3005;
    return;
    }

if ((loD = (int *)malloc((size_t)(*levels)*sizeof(int)))==NULL)   {
    *error = 3006;
    return;
    }





/* Now fill up these arrays */

*lfC = *llC = 0;
cnt = 1;

for(i=1; i<*levels+1; ++i)  {
    *(lfC+i) = 0;
    *(llC+i) = *(llC+i-1) + cnt;
    cnt<<=1;
    } 

*(loC+*levels+1-1) = 0;
for(i=*levels+1-2; i>=0; --i)   {
    *(loC+i) = *(loC+i+1) + *(llC+i+1)+1;
    }

*lfD = *llD = 0;
cnt = 1;

for(i=1; i<*levels; ++i)    {
    *(lfD+i) = 0;
    *(llD+i) = *(llD+i-1) + cnt;
    cnt<<=1;
    }

*(loD+*levels-1) = 0;
for(i=*levels-2; i>=0; --i) {
    *(loD+i) = *(loD+i+1) + *(llD+i+1)+1;
    }

/* Now we have to create the C and D arrays */

*LengthC = *loC + 1;
*LengthD = *loD + 1;

if ((lC = (double *)calloc((size_t)*LengthC,(size_t)sizeof(double)))==NULL) {
    *error = 3007;
    return;
    }

if ((lD = (double *)calloc((size_t)*LengthD,(size_t)sizeof(double)))==NULL) {
    *error = 3008;
    return;
    }
/* Calloc should already zero these arrays */

for(i=0; i<*ndata; ++i)
    *(lC+i) = *(TheData+i);

/* Sorted, now do the wavelet transform */

wavedecomp(lC, lD, H, LengthH, levels,
    lfC, llC, loC, lfD, llD, loD, type, bc, error);

if (*error != 0)    {
    *error = 3009;
    return;
    }

/* Now we can return all the answers. To do this we have to link the information
 * in the l* arrays to the real ones
 */

*C = lC;
*D = lD;
*firstC = lfC;
*lastC = llC;
*offsetC = loC;
*firstD = lfD;
*lastD = llD;
*offsetD = loD;

/* That's it, time to go home */

return;
}

void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
 double *C::       Input data, and the subsequent smoothed data 
 double *D::       The wavelet coefficients                     
 double *H::       The smoothing filter H                       
 int *LengthH::    Length of smoothing filter                   
 int *levels::     The number of levels in this decomposition   
 int *firstC::     The first possible C coef at a given level   
 int *lastC::      The last possible C coef at a given level    
 int *offsetC::    Offset from C[0] for certain level's coeffs  
 int *firstD::     The first possible D coef at a given level   
 int *lastD::      The last possible D coef at a given level    
 int *offsetD::    Offset from D[0] for certain level's coeffs  
 int *type::       The type of wavelet decomposition        
 int *bc::         Method of boundary correction        
 int *error::      Error code                                   
 *---------------------*/
{
register int next_level,at_level;
register int step_factor;   /* Controls width of filter for station */
register int verbose;   /* Controls message printing, passed in error var*/

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

if (*error == 1l)   /* Error switches on verbosity */
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
    }
        
if (verbose) Rprintf("Decomposing into level: ");

*error = 0;

step_factor = 1;    /* This variable should *always* be 1 for standard
             * wavelets. It should start at 1 for stationary
             * wavelets and multiply itself by 2 each stage
             */

for(next_level = *levels - 1; next_level >= 0; --next_level)    {

    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level + 1;

/* For stationary wavelets we need to define a step factor.
 * This widens the span of the filter. At the top level (*levels->*levels-1)
 * it is one, as usual. Then for the next step it becomes 2, then 4 etc.
 */

    convolveC( (C+*(offsetC+at_level)),
        (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
        (int)(*(firstC+at_level)),
        H,
        (int)*LengthH,
        (C+*(offsetC+next_level)),
        (int)(*(firstC+next_level)),
        (int)(*(lastC+next_level)) , (int)*type,
        step_factor, (int)*bc);

    convolveD( (C+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                H,
                (int)*LengthH,
        (D+*(offsetD+next_level)),
        (int)(*(firstD+next_level)),
        (int)(*(lastD+next_level)), (int)*type,
        step_factor, (int)*bc );

    if (*type == STATION)
        step_factor *= 2;   /* Any half decent compiler should
                     * know what to do here ! */
    }
if (verbose)
    Rprintf("\n");
return;
}

void accessDwp(double *Data, int *LengthData, int *nlevels, int *level,
	double *answer, int *error)
/*---------------------
 * Argument description
 *---------------------
double *Data::        This is a 2D array. Top level contains data 
int *LengthData::     Length of Data, this is power of 2              
int *nlevels::        The number of levels in this decomposition   
int *level::          Which level you want to extract      
double *answer::      The level of coefficients            
int *error::          Error code                   
 *---------------------*/
{
register int i;

*error = 0;

/*
 * Check variable integrity
 */

if (*level < 0) {
    *error =4000;
    return;
    }

else if (*level > *nlevels) {
    *error = 4001;
    return;
    }

for(i=0; i< *LengthData; ++i)
    *(answer+i) = ACCESSD(*level, i);


}


void wavepackde(double *Data, int *LengthData, int *levels, double *H,
	int *LengthH)
/*---------------------
 * Argument description
 *---------------------
double *Data::        This is a 2D array. Top level contains data 
int *LengthData::     Length of Data, this is power of 2              
int *levels::         The number of levels, 2^(levels+1)=LengthData   
double *H::           The filter to use                   
int *LengthH::        Length of filter                
 *---------------------*/
{
int startin, outstart1, outstart2;
/*
int i,j;
*/
void wvpkr(double *Data, int startin, int lengthin, int outstart1,
	int outstart2, int level, double *H, int LengthH, int *LengthData);

/*
Rprintf("This routine is wavepackde\n");
Rprintf("Length of data is %ld\n", *LengthData);
Rprintf("Number of levels is %ld\n", *levels);
Rprintf("Data array is:\n");
for(i= (int)*levels; i>=0; --i)
    for(j=0; j< *LengthData; ++j)   {
        Rprintf("Level %d, Item %d is %lf\n", i,j, ACCESSD(i,j));
        }
*/

startin = 0;
outstart1 = 0;
outstart2 = ((int)*LengthData)/2;

wvpkr(Data, startin, (int)*LengthData, outstart1, outstart2, (int)*levels, H,
 (int)*LengthH, LengthData);
}

void wvpkr(double *Data, int startin, int lengthin, int outstart1,
	int outstart2, int level, double *H, int LengthH, int *LengthData)
/* int level;  The level where we're at         */
{
int lengthout;
void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

lengthout = lengthin/2;


convolveC( POINTD(level, startin), lengthin, 0, H, LengthH,
    POINTD(level-1, outstart1), 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

convolveD( POINTD(level, startin), lengthin, 0, H, LengthH,
    POINTD(level-1, outstart2), 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

if (lengthout==1)
    return;
else    {
    /*
     * Now apply both filters to the LOW pass filtered data
     */
    wvpkr(Data, outstart1, lengthout, outstart1, outstart1+lengthout/2,
        level-1, H, LengthH, LengthData); 
    /*
     * Now apply both filters to the HIGH pass filtered data
     */
    wvpkr(Data, outstart2, lengthout, outstart2, outstart2+lengthout/2,
        level-1, H, LengthH, LengthData); 
    }
}
/* WAVEPACKRECON    -   inverse swt             */
/*
 * Error codes
 *
 *  1   -   As the reconstruction is built up the vector
            ldata should contain a doubling sequence (apart
            from the first two numbers which should be the
            same). This error is returned if this is not the case.

 *  2   -   memory error on creating c_in
 *  3   -   memory error on creating c_out
 */


void wavepackrecon(double *rdata, int *ldata, int *nrsteps, int *rvector, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *rdata::	The transformed data, packets are packed together    
int *ldata:: 	Array of lengths of packets in rdata         
int *nrsteps::	The number of reconstruction steps (also length of ldata array)
int *rvector::	Integer whose binary decomp reveals rotate/not instruction
double *H::	Filter                       
int *LengthH::	Length of filter                 
int *error::	Error code                       
 *---------------------*/
{
register int i,j;
register int msb;
register int ldctr;
int LengthCin;
int LengthCout;
int LengthDin;
double *c_in;
double *c_out;
void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);
void rotateback(double *book, int length);

/* Set error code to zero as no error has occured yet!  */

*error = 0;

/* We can use conbar to do all our hard work for us         */
/* This is the reconstruction step in the ordinary DWT. The only    */
/* modification that we have to make is to rotate the data at each  */
/* step if we need to. This information is stored in "rvector" (rotate  */
/* vector). This is a single integer whose information is stored in */
/* binary form. Each bit refers to a rotate/non rotate operation and    */
/* should be applied by the following method:               */
/*
 *
 *  a. do conbar
 *  b. check next most sig bit of rvector and rotate if 1
 *
 * And start with the most significant bit.
 */

/*
Rprintf("Rvector is %d\n", (int)*rvector);
*/

/* First let's generate the MSB */

msb = 0x01 << ((int)*nrsteps-1);

/* Get initial C data stored    */

LengthCin = (int)*(ldata+0);

ldctr = LengthCin;  /* ldctr measures how far aint rdata we have gone */

if ((c_in = (double *)malloc((unsigned)LengthCin*sizeof(double)))==NULL) {
        *error = 2;
        return;
        }

for(j=0; j< LengthCin; ++j)
    *(c_in+j) = *(rdata+j);

LengthCout = LengthCin;

c_out=calloc(LengthCout,sizeof(double));    /* MAN: added initialization.  Hopefully shouldn't have any bad
                                                consequences... */
for(i=0; i< (int)*nrsteps; ++i) {
    LengthCout *= 2;

    if (i != 0)
        free((void *)c_out);

    if ((c_out=(double *)malloc((unsigned)LengthCout*sizeof(double)))==NULL)    {
        *error = 3;
        return;
        }

    /* Now store D data at this level   */
    LengthDin = (int)*(ldata+(i+1));

    /* Don't need to store cos we can put rdata+ldctr straight in for d_in
     * for(j=0; j < LengthDin; ++j)
     *  *(d_in+j) = *(rdata+ldctr+j);
     */

    conbar(c_in, LengthCin, 0,
           rdata+ldctr, LengthDin, 0,
           H, (int)*LengthH,
           c_out, LengthCout, 0, LengthCout-1,
           WAVELET, PERIODIC);

    ldctr += LengthDin; /* update cos we've moved aint rdata */ 

    /* O.k. chaps, c_out must now become c_in, and we should check
       that the lengths match */

    /*
    Rprintf("LengthCout is %d\n", LengthCout);
    Rprintf("i is %d\n", i);
    Rprintf("nrsteps is %d\n", (int)*nrsteps);
    */
    /*Rprintf("ldata+i+2 is %d\n", (int)*(ldata+i+2));*/
    if (i+1 != (int)*nrsteps && LengthCout != (int)*(ldata+i+2))    {
        *error = 1;
        return;
        }

    /* Do we rotate back ? */

    if (msb & (int)*rvector)    {
        /*
        Rprintf("Rotating\n");
        */
        rotateback(c_out, LengthCout);
        }

    /*
    Rprintf("msb is: %d\n", msb);
    */

    msb >>= 1;

    /* Now c_in <- c_out */

    free((void *)c_in);

    if ((c_in = (double *)malloc((unsigned)LengthCout*sizeof(double)))==NULL)   {
        *error = 2;
        return;
        }

    for(j=0; j<LengthCout; ++j)
        *(c_in+j) = *(c_out+j);

    LengthCin = LengthCout;

    }
/* Now copy c_out into rdata - they should be the same length at this stage */

for(j=0; j<LengthCout; ++j)
    *(rdata+j) = *(c_out+j);

free((void *)c_out);
free((void *)c_in);

}
/*
 * Routine to perform the stationary wavelet decomposition in a wavelet
 * packet fashion.
 */

void wavepackst(double *Carray, double *Data, int *LengthData, int *levels,
	double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *Carray::      Will contain bottom most Cs             
double *Data::        This is a 2D array. Zeroeth level contains data 
int *LengthData::     Length of Data, this is power of 2              
int *levels::         The number of levels, 2^(*levels)=LengthData    
double *H::	      The filter to use                   
int *LengthH::        Length of filter                
int *error::          Error code, if non-zero then it's a mem error   
 *---------------------*/
{
int startin, outstart1, outstart2;
register int i;
double *book;

void wvpkstr(double *Carray, double *Data, int startin, int lengthin,
	int outstart1, int outstart2, int level, double *H, int LengthH,
	int *LengthData, double *book, int *error);

*error = 0;

/*
Rprintf("This routine is wavepackst\n");
Rprintf("Length of data is %ld\n", *LengthData);
Rprintf("Number of levels is %ld\n", *levels);
Rprintf("Data array is:\n");
for(i= (int)*levels; i>=0; --i)
    for(j=0; j< *LengthData; ++j)   {
        Rprintf("Level %d, Item %d is %lf\n", i,j, ACCESSD(i,j));
        }
*/

/* Create a bookeeping vector. That contains the C,C' level smooths
   thoughout the algorithm */

if ((book = (double *)malloc((unsigned)*LengthData*sizeof(double)))==NULL){
    *error = 1;
    return;
    }
    
/* Copy original data to book keeping vector */

for(i=0; i< *LengthData; ++i)
    *(book+i) = *POINTD(*levels, i);


startin = 0;
outstart1 = 0;
outstart2 = ((int)*LengthData)/2;

wvpkstr(Carray, Data, startin, (int)*LengthData, outstart1, outstart2,
    (int)*levels, H, (int)*LengthH, LengthData, book, error);

if (*error != 0)
    return;
else
    free((void *)book);
}

void wvpkstr(double *Carray, double *Data, int startin, int lengthin,
	int outstart1, int outstart2, int level, double *H, int LengthH,
	int *LengthData, double *book, int *error)
/* int level;  The level where we're at         */
{
register int i;
int lengthout;
double *book1, *book2;

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);
void rotater(double *book, int length);

/*
Rprintf("wvpkstr entry\n");
Rprintf("lengthout is %d\n", lengthout);
*/

lengthout = lengthin/2;

if ((book1 = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL) {
    *error = 1;
    return;
    }

else if ((book2 = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL){
    *error = 1;
    return;
    }


convolveC(book, lengthin, 0, H, LengthH,
    book1, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);


for(i=0; i < lengthout; ++i)
    * POINTC(level-1, (outstart1+i)) = *(book1+i);

/*
Rprintf("book1 coefficients \n");
for(i=0; i<lengthout; ++i)
    Rprintf("%lf ", *(book1+i));;
Rprintf("\n");
*/

convolveD( book, lengthin, 0, H, LengthH,
    POINTD(level-1, outstart1), 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

/* Now cycle book around one, and do the convolutions again */

/* COMMENT OUT
 * tmp = *(book+lengthin-1);
 *
 * for(i=lengthin-1; i>0; --i)
 *  *(book+i) = *(book+i-1);
 * book = tmp;
 */

/* COMMENT OUT (replaced by rotater function) tmp = *book;
 * for(i=0; i<lengthin-1; ++i)
 *  *(book+i) = *(book+i+1);
 *
 * *(book+lengthin-1) = tmp;
 */

rotater(book, lengthin);

convolveC( book, lengthin, 0, H, LengthH,
    book2, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

for(i=0; i < lengthout; ++i)
    *POINTC(level-1, (outstart2+i)) = *(book2+i);

/*
Rprintf("book2 coefficients \n");
for(i=0; i<lengthout; ++i)
    Rprintf("%lf ", *(book2+i));
Rprintf("\n");
*/


convolveD( book, lengthin, 0, H, LengthH,
    POINTD(level-1, outstart2), 0, lengthout-1, 
    WAVELET, 1, PERIODIC);


if (lengthout!=1)   {
    /*
     * Now pass book1 and book2 to the next level 
     */
    wvpkstr(Carray, Data, outstart1, lengthout, outstart1,
        outstart1+lengthout/2, level-1, H, LengthH, LengthData,
        book1, error); 

    if (*error != 0)
        return;

    wvpkstr(Carray, Data, outstart2, lengthout, outstart2,
        outstart2+lengthout/2, level-1, H, LengthH, LengthData,
        book2, error); 

    if (*error != 0)
        return;
    }


free((void *)book1);
free((void *)book2);
}
/*
 * waverecons:  Do 1D wavelet reconstruction
 */
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
 double *C::       Input data, and the subsequent smoothed data 
 double *D::       The wavelet coefficients                     
 double *H::       The smoothing filter H                       
 int *LengthH::    Length of smoothing filter                   
 int *levels::     The number of levels in this decomposition   
 int *firstC::     The first possible C coef at a given level   
 int *lastC::      The last possible C coef at a given level    
 int *offsetC::    Offset from C[0] for certain level's coeffs  
 int *firstD::     The first possible D coef at a given level   
 int *lastD::      The last possible D coef at a given level    
 int *offsetD::    Offset from D[0] for certain level's coeffs  
 int *type::       The type of wavelet decomposition        
 int *bc::         Method of boundary correction        
 int *error::      Error code                                   
 *---------------------*/
{
register int next_level, at_level;
register int verbose;   /* Printing messages, passed in error       */

void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);

if (*error == 1)
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
    }

if (verbose) Rprintf("Building level: ");

*error = 0;

for(next_level = 1; next_level <= *levels; ++next_level)    {

    
    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level - 1; 

    conbar( (C+*(offsetC+at_level)),
        (int)(*(lastC+at_level) - *(firstC+at_level) + 1),
        (int)(*(firstC+at_level)),
        (D+*(offsetD+at_level)),
        (int)(*(lastD+at_level) - *(firstD+at_level) + 1),
        (int)(*(firstD+at_level)),
        H,
        (int)*LengthH,
        (C+*(offsetC+next_level)),
        (int)(*(lastC+next_level) - *(firstC+next_level)+1),
                (int)(*(firstC+next_level)),
                (int)(*(lastC+next_level)),
        (int)(*type),
        (int)(*bc) );
    }
if (verbose)
    Rprintf("\n");

return;
}

/*
 * Functions to do complex arithmetic
 *
 */

/*
 * Addition: a+ib + c+id = a+c +i(b+d) = e + i f
 */

void comadd(double a, double b, double c, double d, double *e, double *f)
{
*e = a+c;
*f = b+d;
}

/*
 * Subtraction: a+ib - c+id = a+c -i(b+d) = e + i f
 */

void comsub(double a, double b, double c, double d, double *e, double *f)
{
*e = a-c;
*f = b-d; 
}

/*
 * Multiplication: (a+ib)(c+id) = ac-bd +i(bc+ad) = e + i f
 */

void commul(double a, double b, double c, double d, double *e, double *f)
{
*e = (a*c - b*d);
*f = (b*c + a*d);
}


/*
 * Division: (a+ib)(c+id) = (ac+bd +i(bc-ad))/(c^2+d^2) = e + i f
 */

void comdiv(double a, double b, double c, double d, double *e, double *f)
{
double tmp;

tmp = c*c + d*d;

*e = (a*c + b*d)/tmp;
*f = (b*c - a*d)/tmp;
}

/*
 * Complex wavelet version
 */

/*
 * COMCBR: Does the reconstruction convolution
 */

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc)
{
register int n,k;
register int cfactor;
double sumCR, sumCI, sumDR, sumDI;
double a,b,c,d,e,f;

switch(type)    {

    case WAVELET:   /* Standard wavelets */
        cfactor = 2;
        break;

    case STATION:   /* Stationary wavelets */
        cfactor = 1;
        break;

    default:    /* This should never happen */
        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
        break;
    }


/* Compute each of the output C */

for(n=firstCout; n<=lastCout; ++n)  {

    /* We want  n+1-LengthH <= 2*k to start off */


    k = CEIL(n+1-LengthH);

    sumCR = 0.0;
    sumCI = 0.0;
    sumDR = 0.0;
    sumDI = 0.0;

    while( cfactor*k <= n ) {

        a = *(HR + n - cfactor*k);
        b = *(HI + n - cfactor*k);

        c = ACCESSC(c_inR, firstCin, LengthCin, k, bc);
        d = ACCESSC(c_inI, firstCin, LengthCin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumCR += e;
        sumCI += f;

        /* Now D part */

        a = *(GR + n - cfactor*k);
        b = *(GI + n - cfactor*k);

        c = ACCESSC(d_inR, firstDin, LengthDin, k, bc);
        d = ACCESSC(d_inI, firstDin, LengthDin, k, bc);

        commul(a,b,c,d, &e, &f);

        sumDR += e;
        sumDI += f;

        ++k;
        }

    sumCR += sumDR;
    sumCI += sumDI;

    ACCESSC(c_outR, firstCout, LengthCout, n, bc) = sumCR;
    ACCESSC(c_outI, firstCout, LengthCout, n, bc) = sumCI;
    }

}

/*
 * This routine is identical to the convolve.c routine except it
 * does it for complex wavelets.

 * COMCONC  -   Do filter H filter convolution with boundary
 */


void comconC(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *HR, double *HI, int LengthH,
	double *c_outR, double *c_outI,
	int LengthCout, int firstCout, int lastCout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_inR::   Input data (real)                    
double *c_inI::   Input data (imaginary)               
int LengthCin::   Length of this array                 
int firstCin::     <-- MAN: added since missing...     
double *HR::      Lowpass Filter                   
double *HI::      Lowpass Filter                   
int LengthH::     Length of filter                 
double *c_outR::  Output data (real)                   
double *c_outI::  Output data (imaginary)              
int LengthCout::  Length of above array                
int firstCout::   First index of C array               
int lastCout::    Last index of C array                
int type::        Type of wavelet decomposition            
int step_factor:: For stationary wavelets only             
int bc::          Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sumR,sumI;
double a,b,c,d,e,f;
register int k;
register int count_out;
register int m;
register int cfactor;   /* This determines what sort of dilation we do  */
            /* and depends on the type argument     */

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstCout; k<=lastCout; ++k)  {
    sumR = 0.0;
    sumI = 0.0;


    for(m=0; m<LengthH; ++m)    {

        a = *(HR + m);  /* real part */
        b = *(HI + m);  /* imaginary part */

        c = ACCESSC(c_inR, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);

        commul(a,b,c,d,&e, &f);

        sumR += e;
        sumI += f;
        }

    *(c_outR + count_out) = sumR;
    *(c_outI + count_out) = sumI;
    ++count_out;
    }
}

void comconD(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *GR, double *GI, int LengthH,
	double *d_outR, double *d_outI,
	int LengthDout, int firstDout, int lastDout,
	int type, int step_factor, int bc)
/*---------------------
 * Argument description
 *---------------------
double *c_inR::   Input data (real)                    
double *c_inI::   Input data (imaginary)               
int LengthCin::   Length of this array                 
int firstCin::     <-- MAN: added since missing...     
double *GR::      Lowpass Filter                   
double *GI::      Lowpass Filter                   
int LengthH::     Length of filter                 
double *d_outR::  Output data (real)                   
double *d_outI::  Output data (imaginary)              
int LengthDout::  Length of above array                
int firstDout::   First index of C array               
int lastDout::    Last index of C array                
int type::        Type of wavelet decomposition            
int step_factor:: For stationary wavelets only             
int bc::          Method of boundary correction PERIODIC, SYMMETRIC    
 *---------------------*/
{
double sumR, sumI;
double a,b,c,d,e,f;
register int k;
register int count_out;
register int m;
register int cfactor;

count_out = 0;

switch(type)    {

    case WAVELET:   /*  Ordinary wavelets   */
            cfactor = 2;    /* Pick every other coefficient */
            break;

    case STATION:   /* Stationary wavelets  */
            cfactor = 1;    /* Pick every coefficient   */
            break;


    default:    /* This is an error, one of the above must have */
            /* been picked */
            /* However, this must be tested in a previous   */
            /* routine.                 */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
        }

for(k=firstDout; k<=lastDout; ++k)  {
    sumR = 0.0;
    sumI = 0.0;


    for(m=0; m<LengthH; ++m)    {

        a = *(GR+m);

        b = *(GI+m);

        c = ACCESSC(c_inR, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
            ((step_factor*m)+(cfactor*k)),bc);
        /*
        c = ACCESSC(c_inR, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);
        d = ACCESSC(c_inI, firstCin, LengthCin,
                (cfactor*k+(step_factor*(1-m))),bc);

        Rprintf("%d: (%lf, %lf)* (%lf, %lf)\n", a,b,c,d);
        */
        commul(a,b,c,d,&e,&f);

        sumR += e;
        sumI += f;
        }

    *(d_outR + count_out) = sumR;
    *(d_outI + count_out) = sumI;
    ++count_out;
    }
}


/*
 * Complex version of wavelet transform
 */

void comwd(double *CR, double *CI, int *LengthC,
	   double *DR, double *DI, int *LengthD,
	double *HR, double *HI, double *GR, double *GI, int *LengthH,
	int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
double *CR::             Input data, and the subsequent smoothed data 
double *CI::             Input data, and the subsequent smoothed data 
int *LengthC::           Length of C array                            
double *DR::             The wavelet coefficients                     
double *DI::             The wavelet coefficients                     
int *LengthD::           Length of D array                            
double *HR::             The smoothing filter H                       
double *HI::             The smoothing filter H                       
double *GR::             The highpass filter H                       
double *GI::             The highpass filter H                       
int *LengthH::           Length of smoothing filter                   
int *levels::            The number of levels in this decomposition   
int *firstC::            The first possible C coef at a given level   
int *lastC::             The last possible C coef at a given level    
int *offsetC::           Offset from C[0] for certain level's coeffs  
int *firstD::            The first possible D coef at a given level   
int *lastD::             The last possible D coef at a given level    
int *offsetD::           Offset from D[0] for certain level's coeffs  
int *type::       	 The type of wavelet decomposition        
int *bc::        	 Method of boundary correction        
int *error::             Error code                                   
 *---------------------*/
{
register int next_level,at_level;
register int step_factor;   /* Controls width of filter for station */
register int verbose;   /* Controls message printing, passed in error var*/

if (*error == 1)   /* Error switches on verbosity */
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
        break;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
        break;
    }
        
if (verbose) Rprintf("Decomposing into level: ");

*error = 0;

step_factor = 1;    /* This variable should *always* be 1 for standard
             * wavelets. It should start at 1 for stationary
             * wavelets and multiply itself by 2 each stage
             */

for(next_level = *levels - 1; next_level >= 0; --next_level)    {

    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level + 1;

/* For stationary wavelets we need to define a step factor.
 * This widens the span of the filter. At the top level (*levels->*levels-1)
 * it is one, as usual. Then for the next step it becomes 2, then 4 etc.
 */

    comconC( (CR+*(offsetC+at_level)),
           (CI+*(offsetC+at_level)),
        (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
        (int)(*(firstC+at_level)),
        HR, HI,
        (int)*LengthH,
        (CR+*(offsetC+next_level)),
        (CI+*(offsetC+next_level)),
        (int)(*(lastC+next_level) - *(firstC+next_level)+1),
        (int)(*(firstC+next_level)),
        (int)(*(lastC+next_level)) , (int)*type,
        step_factor, (int)*bc);

    comconD( (CR+*(offsetC+at_level)),
           (CI+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                GR, GI,
                (int)*LengthH,
        (DR+*(offsetD+next_level)),
        (DI+*(offsetD+next_level)),
        (int)(*(lastD+next_level) - *(lastD+next_level)+1),
        (int)(*(firstD+next_level)),
        (int)(*(lastD+next_level)), (int)*type,
        step_factor, (int)*bc );

    if (*type == STATION)
        step_factor *= 2;   /* Any half decent compiler should
                     * know what to do here ! */
    }
if (verbose)
    Rprintf("\n");
return;
}

/*
 * comwr:  Do 1D complex wavelet reconstruction
 */

void comwr(double *CR, double *CI, int *LengthC,
	   double *DR, double *DI, int *LengthD,
	double *HR, double *HI, double *GR, double *GI, int *LengthH,
	int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/*---------------------
 * Argument description
 *---------------------
double *CR::             Input data, and the subsequent smoothed data 
double *CI::             Input data, and the subsequent smoothed data 
int *LengthC::           Length of C array                            
double *DR::             The wavelet coefficients                     
double *DI::             The wavelet coefficients                     
int *LengthD::           Length of D array                            
double *HR::             The smoothing filter H                       
double *HI::             The smoothing filter H                       
double *GR::             The highpass filter H                       
double *GI::             The highpass filter H                       
int *LengthH::           Length of smoothing filter                   
int *levels::            The number of levels in this decomposition   
int *firstC::            The first possible C coef at a given level   
int *lastC::             The last possible C coef at a given level    
int *offsetC::           Offset from C[0] for certain level's coeffs  
int *firstD::            The first possible D coef at a given level   
int *lastD::             The last possible D coef at a given level    
int *offsetD::           Offset from D[0] for certain level's coeffs  
int *type::       	 The type of wavelet decomposition        
int *bc::        	 Method of boundary correction        
int *error::             Error code                                   
 *---------------------*/

{
register int next_level, at_level;
register int verbose;   /* Printing messages, passed in error       */

if (*error == 1)
    verbose = 1;
else
    verbose = 0;

switch(*bc) {

    case PERIODIC:  /* Periodic boundary conditions */
        if (verbose) Rprintf("Periodic boundary method\n");
        break;

    case SYMMETRIC: /* Symmetric boundary conditions */
        if (verbose) Rprintf("Symmetric boundary method\n");
        break;

    default:    /* The bc must be one of the above */
        Rprintf("Unknown boundary correction method\n");
        *error = 1;
        return;
        break;
    }

switch(*type)   {

    case WAVELET:   /* Standard wavelets */
        if (verbose) Rprintf("Standard wavelet decomposition\n");
        break;

    case STATION:   /* Stationary wavelets */
        if (verbose) Rprintf("Stationary wavelet decomposition\n");
        break;

    default:    /* The type must be of one the above */
        if (verbose) Rprintf("Unknown decomposition type\n");
        *error = 2;
        return;
        break;
    }

if (verbose) Rprintf("Building level: ");

*error = 0;

for(next_level = 1; next_level <= *levels; ++next_level)    {

    
    if (verbose)
        Rprintf("%d ", next_level);

    at_level = next_level - 1; 

    comcbr( (CR+*(offsetC+at_level)),
        (CI+*(offsetC+at_level)),
        (int)(*(lastC+at_level) - *(firstC+at_level) + 1),
        (int)(*(firstC+at_level)),
        (int)(*(lastC+at_level)),
        (DR+*(offsetD+at_level)),
        (DI+*(offsetD+at_level)),
        (int)(*(lastD+at_level) - *(firstD+at_level) + 1),
        (int)(*(firstD+at_level)),
        (int)(*(lastD+at_level)),
        HR, HI, GR, GI,
        (int)*LengthH,
        (CR+*(offsetC+next_level)),
        (CI+*(offsetC+next_level)),
        (int)(*(lastC+next_level) - *(firstC+next_level)+1),
                (int)(*(firstC+next_level)),
                (int)(*(lastC+next_level)),
        (int)(*type),
        (int)(*bc) );
    }
if (verbose)
    Rprintf("\n");

return;

}
/*
 * Emulate the WavDE function in C (but not plotting information)
 * and don't return the wavelet coefficients.
 */

#define HARDTHRESH(w,t) ( fabs((w)) > (t) ? (w) : (0.0))

void CWavDE(double *x, int *n, double *minx, double *maxx, int *Jmax,
	double *threshold, double *xout, double *fout, int *nout,
	double *PrimRes,
	double *SFx, double *SFy, int *lengthSF,
	double *WVx, double *WVy, int *lengthWV,
	int *kmin, int *kmax, int *kminW, int *kmaxW,
	double *xminW, double *xmaxW,
	double *phiLH, double *phiRH, double *psiLH, double *psiRH,
	int *verbose, int *error)
/*---------------------
 * Argument description
 *---------------------
double *x::   		The data                     
int *n::  		The length of the data               
double *minx::		The min of the data                  
double *maxx::		The max of the data                  
int *Jmax::		The number of levels in the expansion        
double *threshold::	Threshold value for thresholding the wv coefs
	vvv Output Variables vvv
double *xout::    	The grid on which the density estimate is defined     
double *fout::    	The density estimate defined on the above grid    
int *nout::		The length of the grid               
	vvv Input variables again vvv
double *PrimRes:: 	The primary resolution               
double *SFx::     	The grid on which the scaling function is defined    
double *SFy::     	The scaling function                 
int *lengthSF::   	The length of the grid               
double *WVx::     	The grid on which the wavelet is defined     
double *WVy::     	The wavelet function                 
int *lengthWV::   	The length of the grid               
int *kmin::   		minimum k for scaling function coefficient comp. 
int *kmax::   		maximum k for scaling function coefficient comp. 
int *kminW::  		as above but for each wavelet level (1:Jmax)     
int *kmaxW::  		as above but for each wavelet level (1:Jmax)     
double *xminW::   	minimum x value for each level for wavelet       
double *xmaxW::   	maximum x value for each level for wavelet       
double *phiLH::   	left hand end of support of Phi          
double *phiRH::   	right hand end of support of Phi         
double *psiLH::   	left hand end of support of psi          
double *psiRH::   	right hand end of support of psi         
int *verbose::    	Print messages or not?               
int *error::  		Error codes                      
 *---------------------*/

/* Error codes

    0   -   O.k.
    1   -   Memory error

 */
{
register int i,k,l,j,twopowjp1;
register int la;
double atmp;
double *a;
double sum;
double divisor;
double widthSF,widthWV;
double evalF(double *Fx, double *Fy, int *lengthF, double widthF, double x);
double xmin, xmax; /* Note these are not the same as maxx and minx */
double SFYscale, WVYscale; /* I forgot to multiply by p^{1/2} etc. */

if (*verbose > 1)
    Rprintf("Entered CWavDE function\n");

*kmin = (int)floor(*minx - *phiRH/ *PrimRes);
*kmax = (int)ceil(*maxx - *phiLH/ *PrimRes);

if (*verbose > 1)
    Rprintf("kmin is %d, kmax is %d\n", *kmin, *kmax);       /*MAN: changed %ld to %d since declared as int (L1329) */

la = (int)(*kmax - *kmin) + 1;

if ((a = (double *)malloc((unsigned)(sizeof(double)*la)))==NULL)    {
    *error = 1;
    return;
    }

/* Now compute the widths of the wavelet/scaling function supports */

widthSF = *(SFx+(int)*lengthSF-1) - *SFx;
widthWV = *(WVx+(int)*lengthWV-1) - *WVx;

/* 
 * Now work out all of the scaling function coefficients
 */

k = (int)*kmin;

/* I forgot to multiply by p^{1/2} ! */
SFYscale = sqrt(*PrimRes);

for (i=0; i<la; ++i)    {
    sum = 0.0;
    for(l=0; l<(int)*n; ++l)    {
        sum += evalF(SFx, SFy, lengthSF, widthSF,
            (*PrimRes* *(x+l))-(double)k);
        }
    *(a+i) = SFYscale*sum/(double)*n;
    ++k;
    }
/*
 * Now compute the wavelet supports
 */

for(j=0; j< (int)*Jmax; ++j)    {

    twopowjp1 = 1 << (j+1);     /* MAN: added parentheses for bit shift */

    divisor = *PrimRes*(double)twopowjp1;

    *(kminW+j) = (int)floor(*minx - *psiRH/divisor);
    *(kmaxW+j) = (int)ceil(*maxx - *psiLH/divisor);
    *(xminW+j) = (double)*(kminW+j) + *psiLH/divisor;
    *(xmaxW+j) = (double)*(kmaxW+j) + *psiRH/divisor;
    }

/* Now figure out range of x values over which the density estimate
 * is compactly supported
 */

xmin = (double)*kmin + *phiLH/ *PrimRes;
xmax = (double)*kmax + *phiRH/ *PrimRes;

for(j=0; j< (int)*Jmax; ++j)    {
    if (*(xminW + j) < xmin)
        xmin = *(xminW + j);

    if (*(xmaxW + j) > xmax)
        xmax = *(xmaxW + j);
    }

divisor = (xmax-xmin)/(double)(*nout-1);

for(i=0; i< (int)*nout; ++i)    {
    *(fout+i) = 0.0;
    *(xout+i) = xmin + (double)i*divisor; 
    }

k = *kmin;


for (i=0; i<la; ++i)    {
    for(l=0; l<(int)*nout; ++l) {
        *(fout+l) += *(a+i) *
            evalF(SFx, SFy, lengthSF, widthSF,
                ((*PrimRes * *(xout+l))-(double)k));
        }
    ++k;
    }

for(l=0; l<(int)*nout; ++l)
    *(fout+l) = SFYscale * *(fout+l);

/* That was the easy part. Now we have to repeat the same operations
 * as above for the wavelets themselves
 *
 * Go round in a loop and get the wavelet coefficients for each level
 * reuse a
 */

free((void *)a);

for(j=0; j<(int)*Jmax; ++j) {
    if (*verbose > 0)
        Rprintf("Wavelet step: level %d\n", j);

    twopowjp1 = 1 << (j+1);           /* MAN: added parentheses for bit shift */

    divisor = *PrimRes*(double)twopowjp1;

    WVYscale = sqrt(divisor);

    la = (int)(*(kmaxW+j) - *(kminW+j)) + 1;

    if ((a = (double *)malloc((unsigned)(sizeof(double)*la)))==NULL) {
        *error = 1;
        return;
        }

    /* Now compute the coefficients for this level j */

    k = *(kminW+j);
    for(i=0; i<la; ++i) {
        sum = 0.0;
        for(l=0; l<(int)*n; ++l)    {
            sum += evalF(WVx, WVy, lengthWV, widthWV,
               ((double)twopowjp1* *PrimRes * *(x+l))
               -(double)k);
            }
        sum *= WVYscale;
        *(a+i) = HARDTHRESH(sum/(double)*n, *threshold);
        atmp = WVYscale * *(a+i);
        for(l=0; l<(int)*nout; ++l)     {
            *(fout+l) += atmp *
                evalF(WVx, WVy, lengthWV, widthWV,
                ((double)twopowjp1* *PrimRes * *(xout+l))-(double)k);
            }
        ++k;
        }
    free((void *)a);
    }



*error = 0;
}

void SCevalF(double *Fx, double *Fy, int *lengthF, double *widthF,
	double *x, int *nx, double *answer)
/*---------------------
 * Argument description
 *---------------------
double *Fx::	  Grid upon which function is defined          
double *Fy::	  Function definition                  
int *lengthF::    Length of above grids                
double *widthF::  Width end Fx - start Fx              
double *x::       Vector x to evaluate function at         
int *nx::         Length of x                      
double *answer::  The answer (again vector of length x)        
 *---------------------*/
{
register int i;
double evalF(double *Fx, double *Fy, int *lengthF, double widthF, double x);

for(i=0; i< (int)*nx; ++i)  {
    *(answer+i) = evalF(Fx, Fy, lengthF, *widthF, *(x+i));
    }
}


/*
 * Evaluate the function Fx,Fy at x. Function is deemed to be zero outside
 * of range of Fx
 *
 * Fx must be a strictly increasing equally spaced design
 */

double evalF(double *Fx, double *Fy, int *lengthF, double widthF, double x)
/*---------------------
 * Argument description
 *---------------------
double *Fx::	  Grid upon which function is defined          
double *Fy::	  Function definition                  
int *lengthF::    Length of above grids                
double widthF::  Width end Fx - start Fx              
double x::       Vector x to evaluate function at         
 *---------------------*/
{
register int il,ir;
double a;
double fp;

/*
 * First get approximate index of point to return
 */

if (x < *Fx || x > *(Fx + (int)*lengthF - 1))
    return(0.0);

/*
 * From VALGRIND check, changed the next line from this
a = (double)((int)*lengthF - 1) * (x - *Fx)/widthF;
 *
 * To this next one immediately after this comment.
 * This is because if x was equal to the RH end then the ratio would be
 * one, and then il would be the last element of the array and ir could be
 * OUTSIDE the array.
 */

a = (double)((int)*lengthF - 2) * (x - *Fx)/widthF;

/* Now a should always be >= 0, since we've already rejected any
   possible negatives, so we don't have to use floor & ceil here
   Just (int) will do. */


il = (int)a;
ir = il+1;

fp = a - (double)il;

return( ((1.0-fp)* *(Fy+il)) + (fp* *(Fy+ir)) );
}

#define MAX(a,b)    ( (a) < (b) ? (b) : (a))
#define MIN(a,b)    ( (a) < (b) ? (a) : (b))

void CScalFn(double *v, double *ans, int *res, double *H, int *lengthH)
{
register int k,n;
double sum;
int b,e;

for(n=0; n< (int)*res; ++n) {
    sum = 0.0;
    b = MAX(0, (int )ceil( ((double)(n+1- *lengthH))/2.0));
        e = MIN(*res, (int )floor(((double) n)/2.0));
    for(k=b; k<= e; ++k)    {
        sum += *(H+n-2*k) * *(v+k);
        }
    *(ans+n) = sum;
    }
}

/* Perform tensor product wavelet transform */

void tpwd(double *image, int *nrow, int *ncol, int *levr, int *levc, 
    int *firstCr, int *lastCr, int *offsetCr,
    int *firstDr, int *lastDr, int *offsetDr,
    int *firstCc, int *lastCc, int *offsetCc,
    int *firstDc, int *lastDc, int *offsetDc,
    int *type, int *bc,
    double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *image::   The image to decompose           
int *nrow::       The number of rows in the image      
int *ncol::       The number of cols in the image      
int *levr::       The number of levels as rows in the image    
int *levc::       The number of levels as cols in the image    
int *firstCr::    The first possible C coef at a given level   
int *lastCr::     The last possible C coef at a given level    
int *offsetCr::   Offset from C[0] for certain level's coeffs  
int *firstDr::    The first possible D coef at a given level   
int *lastDr::     The last possible D coef at a given level    
int *offsetDr::   Offset from D[0] for certain level's coeffs  
int *firstCc::    The first possible C coef at a given level   
int *lastCc::     The last possible C coef at a given level    
int *offsetCc::   Offset from C[0] for certain level's coeffs  
int *firstDc::    The first possible D coef at a given level   
int *lastDc::     The last possible D coef at a given level    
int *offsetDc::   Offset from D[0] for certain level's coeffs  
int *type::       The type of wavelet decomposition            
int *bc::         Method of boundary correction                
double *H::       The wavelet filter               
int *LengthH::    The length of the wavelet filter     
int *error::      0=no error, various errors possible      
 *---------------------*/
{
register int i,j;

double *C;      /* temporary store for input/output data    */
double *D;      /* temporary store for wavelet coefficients */


void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);

*error = 0;

if ((C = (double *)malloc(2*(unsigned)*ncol *sizeof(double)))==NULL)    {
    *error = 1;
    return;
    }
if ((D = (double *)malloc((unsigned)*ncol *sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }

/* First do the wavelet transform across all rows in the image for each row
 */

for(i=0; i< *nrow; ++i) {

    /* Copy the row across - sorry there is probably a more efficient
       way to do this in-place, but what the hell. */


    for(j=0; j< *ncol; ++j) {
        *(D+j) = 0.0;
        *(C+j) = ACCESS(image, *ncol, i, j);
        }

    /* Now do the jolly old wavelet transform */


    wavedecomp(C, D, H, LengthH, levc,
        firstCc, lastCc, offsetCc,
        firstDc, lastDc, offsetDc,
        type, bc, error);


    if (*error != 0)
        return;

    /* And put the answers back in the image array  */

    ACCESS(image, *ncol, i, 0) = *(C+ (*ncol*2)-2);


    for(j=1; j< *ncol; ++j) {
        ACCESS(image, *ncol, i, j) = *(D+j-1);
        } 

    }

free(C);
free(D);

/* Now do it the other way around */

if ((C = (double *)malloc(2*(unsigned)*nrow *sizeof(double)))==NULL)    {
    *error = 1;
    return;
    }
if ((D = (double *)malloc((unsigned)*nrow *sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }
    
/* Second do the wavelet transform across all cols in the image for each col
 */

for(j=0; j< *ncol; ++j) {

    /* Copy the row across - sorry there is probably a more efficient
       way to do this in-place, but what the hell. */


    for(i=0; i< *nrow; ++i) {
        *(D+i) = 0.0;
        *(C+i) = ACCESS(image, *ncol, i, j);
        }

    /* Now do the jolly old wavelet transform */


    wavedecomp(C, D, H, LengthH, levr,
        firstCr, lastCr, offsetCr,
        firstDr, lastDr, offsetDr,
        type, bc, error);


    if (*error != 0)
        return;

    /* And put the answers back in the image array  */

    ACCESS(image, *ncol, 0, j) = *(C+ (*nrow*2)-2);


    for(i=1; i< *nrow; ++i) {
        ACCESS(image, *ncol, i, j) = *(D+i-1);
        } 

    }

free(C);
free(D);
    
}

/* Inverse tensor product wavelet transform             */

void tpwr(double *image, int *nrow, int *ncol, int *levr, int *levc, 
    int *firstCr, int *lastCr, int *offsetCr,
    int *firstDr, int *lastDr, int *offsetDr,
    int *firstCc, int *lastCc, int *offsetCc,
    int *firstDc, int *lastDc, int *offsetDc,
    int *type, int *bc,
    double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *image::   The tpwd coefficients to reconstruct     
int *nrow::       The number of rows in the image      
int *ncol::       The number of cols in the image      
int *levr::       The number of levels as rows in the image    
int *levc::       The number of levels as cols in the image    
int *firstCr::    The first possible C coef at a given level   
int *lastCr::     The last possible C coef at a given level    
int *offsetCr::   Offset from C[0] for certain level's coeffs  
int *firstDr::    The first possible D coef at a given level   
int *lastDr::     The last possible D coef at a given level    
int *offsetDr::   Offset from D[0] for certain level's coeffs  
int *firstCc::    The first possible C coef at a given level   
int *lastCc::     The last possible C coef at a given level    
int *offsetCc::   Offset from C[0] for certain level's coeffs  
int *firstDc::    The first possible D coef at a given level   
int *lastDc::     The last possible D coef at a given level    
int *offsetDc::   Offset from D[0] for certain level's coeffs  
int *type::       The type of wavelet decomposition            
int *bc::         Method of boundary correction                
double *H::       The wavelet filter               
int *LengthH::    The length of the wavelet filter     
int *error::      0=no error, various errors possible      
 *---------------------*/
{
register int i,j;

double *C;      /* temporary store for input data       */
double *D;      /* temporary store for wavelet coefficients */

/* The 1D wavelet reconstruction function   */
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);
/* Basically just do tpwd backwards! */

*error = 0;

if ((C = (double *)malloc(2*(unsigned)*nrow *sizeof(double)))==NULL)    {
    *error = 1;
    return;
    }
if ((D = (double *)malloc((unsigned)*nrow *sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }

/*
 * First do the wavelet reconstruction over all cols in the image for each col
 */

for(j=0; j< *ncol; ++j) {

    /* Copy the row across - sorry there is probably a more efficient
       way to do this in-place, but what the hell. */

    *(C+ (*nrow*2)-2) = ACCESS(image, *ncol, 0, j);

    for(i=1; i< *nrow; ++i) {
        *(D+i-1) = ACCESS(image, *ncol, i, j); 
        } 


    /* Now do the jolly old wavelet RECONSTRUCTION */

    waverecons(C, D, H, LengthH, levc,
        firstCc, lastCc, offsetCc,
        firstDc, lastDc, offsetDc,
        type, bc, error);


    if (*error != 0)
        return;

    /* And put the answers back in the image array  */

    for(i=0; i< *nrow; ++i)
        ACCESS(image, *ncol, i, j) = *(C+i);

    }

free(C);
free(D);

if ((C = (double *)malloc(2*(unsigned)*ncol *sizeof(double)))==NULL)    {
    *error = 1;
    return;
    }
if ((D = (double *)malloc((unsigned)*ncol *sizeof(double)))==NULL)  {
    *error = 2;
    return;
    }

/*
 * Second do the wavelet reconstruction over all rows in the image for each row
 */

for(i=0; i< *nrow; ++i) {

    /* Copy the row across - sorry there is probably a more efficient
       way to do this in-place, but what the hell. */

    *(C+ (*ncol*2)-2) = ACCESS(image, *ncol, i, 0); 

    for(j=1; j< *ncol; ++j) {
        *(D+j-1) = ACCESS(image, *ncol, i, j);
        } 


    /* Now do the jolly old wavelet reconstruction */


    waverecons(C, D, H, LengthH, levr,
        firstCr, lastCr, offsetCr,
        firstDr, lastDr, offsetDr,
        type, bc, error);


    if (*error != 0)
        return;

    /* And put the answers back in the image array  */

    for(j=0; j< *ncol; ++j)
        ACCESS(image, *ncol, i, j) = *(C+j);

    }

free(C);
free(D);
}


#define ZILCHTOL    1.0E-300    /* Zero tolerance for Shannon entropy */
#define STOP        1       /* Code for stopping        */
#define LEFT        2       /* Code for going left      */
#define RIGHT       3       /* Code for going right     */

/* Compute Shannon-Weaver entropy substitute - the l^2 log (l^2) "norm" */
void ShannonEntropy(double *v, int *lengthv, double *zilchtol, double *answer,
	int *error)
{
register int i;
double *vsq;
double sum=0.0;
double SW=0.0;

/* Make private copy of squared coefficients    */

*error = 0;

if ((vsq = (double *)malloc((unsigned)*lengthv*sizeof(double)))==NULL)  {
    *error = 15000;
    return;
    }

for(i=0; i < *lengthv; ++i) {
    *(vsq + i) = *(v+i) * *(v+i);
    sum += *(vsq+i);
    if ( *(vsq+i) == 0.0)
        *(vsq+i) = 1.0;
    SW += *(vsq+i) * log(*(vsq+i));
    }

if (sum < *zilchtol)
    *answer = 0.0;

else
    *answer = -SW;

free(vsq);

return;

}

#define ACCESSU(uvec, fv, lev, j)   *(uvec + *(fv+lev) + j)

void Cmnv(double *wst, double *wstC, int *LengthData, int *nlevels,
	int *upperctrl, double *upperl, int *firstl, int *verbose, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wst::         Table of wavelet packet coefficients     
double *wstC::        Table of scaling function coefficients   
int *LengthData::     Length of original data set          
int *nlevels::        Number of levels in the decomposition    
int *upperctrl::      Vector to record "control" decisions     
double *upperl::      Vector to record minimum entropies       
int *firstl::         Index vector into previous two vectors   
int *verbose::        Print out verbose messages (1=yes, 0=no) 
int *error::          Error condition              
 *---------------------*/
{
register int i,j,k;
register int nll, nul;  /* Number of packets in lower and upper levels  */
register int kl, kr;    /* Daughter packet indices, left and right  */
int PacketLength;   /* Generic packet lengths           */
double *pkt, *pktl, *pktr;  /* Generic packets          */
double *cpkt;       /* Combined packet for level zero computations  */
double mpE, dlE, drE;   /* Entropies for mother and left & right daughters */
double zilchtol;    /* A zero tolerance             */

double *getpacket(double *wst, int nlevels, int level, int index, int *error);
void ShannonEntropy(double *v, int *lengthv, double *zilchtol, double *answer,
	int *error);

*error = 0;
zilchtol = ZILCHTOL;

if (*verbose == 1)
    Rprintf("Cmnv: function entered\n");

nll = (int)*LengthData;
nul = nll >> 1;

/* Go through each level. i refers to the lower level */

for(i=0; i <= *nlevels-1; ++i)  {
    if (*verbose==1)
        Rprintf("Cmnv: Packets. Lower: %d Upper %d\n", nll, nul);

    for(j=0; j<nul; ++j)    {
        if (*verbose==1)
            Rprintf("Upper level index: %d\n", j);
        kl = j << 1;
        kr = kl + 1;

        /* Get mother packet */

        pkt = getpacket(wstC, 1+(int)*nlevels, i+1, j, error);
        PacketLength = 1 << (i+1);

        if (*error != 0)
            return;

        /* Compute mother packet entropy */

        ShannonEntropy(pkt, &PacketLength, &zilchtol, &mpE, error);

        if (*error != 0)
            return;

        /* Don't need mummy packet now, so free her */

        free((void *)pkt);

        /* Now get daughter packets and entropies. The procedure
         * for this depends on whether we're at the bottom level
         * or not
         */

        if (i==0)   {

            pkt = getpacket(wst, 1+(int)*nlevels, i, kl, error);
            PacketLength = 1 << i;

            if (*error != 0)
                return;

            pktl = getpacket(wstC, 1+(int)*nlevels, i, kl, error);

            if (*error != 0)
                return;

            /* Now create a combined packet of both these two */

            if ((cpkt = (double *)malloc((unsigned)2*PacketLength*sizeof(double)))==NULL) {
                *error = 4;
                return;
                }

            for(k = 0; k<PacketLength; ++k) {
                *(cpkt+k) = *(pkt+k);
                *(cpkt+PacketLength+k) = *(pktl+k);
                }

            /* And work out the entropy of this combined packet */

            PacketLength <<= 1;

            ShannonEntropy(cpkt, &PacketLength, &zilchtol, &dlE, error);
            if (*error != 0)
                return;

            /* Don't need combined packets or daughter left pax any
             * more. We can use cpkt again, so don't free it now
             */

            free((void *)pkt);
            free((void *)pktl);

            /* Now do the same for the daughter right packets */

            pkt = getpacket(wst, 1+(int)*nlevels, i, kr, error);
            PacketLength = 1 << i;

            if (*error != 0)
                return;


            pktr = getpacket(wstC, 1+(int)*nlevels, i, kr, error);

            if (*error != 0)
                return;

            for(k = 0; k<PacketLength; ++k) {
                *(cpkt+k) = *(pkt+k);
                *(cpkt+PacketLength+k) = *(pktr+k);
                }

            /* And work out the entropy of this combined packet */

            PacketLength <<= 1;

            ShannonEntropy(cpkt, &PacketLength, &zilchtol, &drE, error);
            if (*error != 0)
                return;

            /* Don't need combined packets or daughter left pax any
             * more
             */

            free((void *)cpkt);
            free((void *)pkt);
            free((void *)pktr);
            }
        else    {   /* We're not at the bottom */

            if (*verbose==1)
              Rprintf("Left Ent C contrib %lf\n",
                ACCESSU(upperl,firstl, i-1, kl));

            pktl = getpacket(wst, 1+(int)*nlevels, i, kl, error);
            PacketLength = 1 << i;

            if (*error != 0)
                return;

            ShannonEntropy(pktl, &PacketLength, &zilchtol, &dlE,
                error);

            if (*error != 0)
                return;

    
            dlE += ACCESSU(upperl,firstl, i-1, kl);

            /* Loose the daughter packet */

            free((void *)pktl);


            /* now the right daughter */

            if (*verbose==1)
              Rprintf("Right Ent C contrib %lf\n",
                ACCESSU(upperl,firstl, i-1, kr));

            pktr = getpacket(wst, 1+(int)*nlevels, i, kr, error);
            PacketLength = 1 << i;

            if (*error != 0)
                return;

            ShannonEntropy(pktr, &PacketLength, &zilchtol, &drE,
                error);

            if (*error != 0)
                return;

    
            drE += ACCESSU(upperl,firstl, i-1, kr);

            /* Loose the daughter packet */

            free((void *)pktr);
            }

        if (*verbose==1)    {
            Rprintf("Mother ent.: %lf\n", mpE);
            Rprintf("Daug. l. ent.: %lf\n", dlE);
            Rprintf("Daug. r. ent.: %lf\n", drE);
            }


        if (mpE < dlE)
          if (mpE < drE)    {
            ACCESSU(upperl, firstl, i, j) = mpE;
            ACCESSU(upperctrl, firstl, i, j) = STOP; 
            }
          else  {
            ACCESSU(upperl, firstl, i, j) = drE;
            ACCESSU(upperctrl, firstl, i, j) = RIGHT; 
            }
        else if (dlE < drE) {
            ACCESSU(upperl, firstl, i, j) = dlE;
            ACCESSU(upperctrl, firstl, i, j) = LEFT; 
            }
          else  {
            ACCESSU(upperl, firstl, i, j) = drE;
            ACCESSU(upperctrl, firstl, i, j) = RIGHT; 
            }

        if (*verbose==1)
          Rprintf("\tSelected %d %lf\n", ACCESSU(upperctrl,firstl,i,j),
            ACCESSU(upperl, firstl, i, j));
        }

    nul >>= 1;
    nll >>= 1;
    }
}

/*
 * Wavelet packet node vector computations
 *
 * (from ~guy/projects/WAVELETS/PACKET/wpCmnv.c)
 */


#define TOP 1
#define BOTTOM  2

void wpCmnv(double *wp, int *LengthData, int *nlevels,
        int *upperctrl, double *upperl, int *firstl, int *verbose, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wp::          Table of wavelet packet coefficients
int *LengthData::     Length of original data set
int *nlevels::        Number of levels in the decomposition
int *upperctrl::      Vector to record "control" decisions
double *upperl::      Vector to record minimum entropies
int *firstl::         Index vector into previous two vectors
int *verbose::        Print out verbose messages (1=yes, 0=no)
int *error::          Error condition
 *---------------------*/
{
register int i,j;
register int nll, nul;  /* Number of packets in lower and upper levels  */
register int kl, kr;    /* Daughter packet indices, left and right  */
int PacketLength;   /* Generic packet lengths           */
double *pkt, *pktl, *pktr;  /* Generic packets          */
double mpE, dE;     /* Entropies for mother and daughters */
double zilchtol;    /* A zero tolerance             */
double tmp;     /* Temporary holder             */

double *getpacket(double *wst, int nlevels, int level, int index, int *error);
void ShannonEntropy(double *v, int *lengthv, double *zilchtol, double *answer,
	int *error);

*error = 0;
zilchtol = ZILCHTOL;

if (*verbose == 1)
    Rprintf("wpCmnv: function entered\n");

nll = (int)*LengthData;
nul = nll >> 1;

/* Go through each level. i refers to the lower level */

for(i=0; i <= *nlevels-1; ++i)  {
    if (*verbose==1)
        Rprintf("wpCmnv: Packets. Lower: %d Upper %d\n", nll, nul);

    for(j=0; j<nul; ++j)    {
        if (*verbose==1)
            Rprintf("Upper level index: %d\n", j);
        kl = j << 1;
        kr = kl + 1;

        /* Get mother packet */

        pkt = getpacket(wp, 1+(int)*nlevels, i+1, j, error);
        PacketLength = 1 << (i+1);

        if (*error != 0)
            return;

        /* Compute mother packet entropy */

        ShannonEntropy(pkt, &PacketLength, &zilchtol, &mpE, error);

        if (*error != 0)
            return;

        /* Don't need mummy packet now, so free her */

        free((void *)pkt);

        /*
         * Now get daughter packets and entropies.
         *
         * If this is at level 0 then it is really the entropies
         *
         * Otherwise the daugther entropy is just the minimum
         * that we decided on before
         * 
         */
        if (i==0)   {

            pktl = getpacket(wp, 1+(int)*nlevels, i, kl, error);
            PacketLength = 1 << i;

            if (*error != 0)
                return;

            pktr = getpacket(wp, 1+(int)*nlevels, i, kr, error);

            if (*error != 0)
                return;

            /* And work out the entropy of the left and right and
             * add them
             */


            ShannonEntropy(pktl, &PacketLength, &zilchtol, &dE, error);
            if (*error != 0)
                return;

            tmp = dE;

            ShannonEntropy(pktr, &PacketLength, &zilchtol, &dE, error);
            if (*error != 0)
                return;

            dE += tmp;


            /* Don't need combined packets or daughter left pax any
             * more. We can use cpkt again, so don't free it now
             */

            free((void *)pktl);
            free((void *)pktr);
        }
        else    {
            dE =  ACCESSU(upperl, firstl, i-1, kl);
            dE += ACCESSU(upperl, firstl, i-1, kr);
        }

        if (*verbose==1)    {
            Rprintf("Mother ent.: %lf\n", mpE);
            Rprintf("Daug. ent.: %lf\n", dE);
            }


        if (mpE < dE)   {
            ACCESSU(upperl, firstl, i, j) = mpE;
            ACCESSU(upperctrl, firstl, i, j) = TOP; 
            }
        else    {
            ACCESSU(upperl, firstl, i, j) = dE;
            ACCESSU(upperctrl, firstl, i, j) = BOTTOM; 
            }

        if (*verbose==1)
          Rprintf("\tSelected %d %lf\n", ACCESSU(upperctrl,firstl,i,j),
            ACCESSU(upperl, firstl, i, j));
        }

    nul >>= 1;
    nll >>= 1;
    }
}

/*
 * WPST -   Stationary wavelet packet algorithm (i.e "The nightmare")
 */

void wpst(double *ansvec, int *lansvec, int *nlev, int *finish_level,
	int *avixstart, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *ansvec::	Vector of length *lansvec that contains the original
			data and will contain the stationary wavelet packet
			coefficients on exit
int *lansvec::		Length of the ansvec vector              
int *nlev::		The number of levels in this transform       
int *finish_level::	The last level to decompose to            
int *avixstart::	A vector of length (*nlev+1). The index ranging from
		 0 to *nlev. The entries in this vector are indices into
		 the ansvec vector indicating the start index for
		 packets for a given level.
		 e.g.  *(avixstart + 0) = 0 (always). So that the first
		 index for level 0 packets in ansvec is 0.
         
		 e.g. *(avixstart+1)=256 (for *nlev=4). So that the first
		 index for level 1 packets in ansvec is 256 etc.
         
double *H::		Filter smoothing coefficients as with all other algs 
int *LengthH::		The number of filter smoothing coefficients      
int *error::		Error code. 0=o.k.                   
			  1 - memory error in creating c_in       
			  2-5 - memory error for c_out, d_out, c_outR, d_outR 
 *---------------------*/
{
register int i,j,k, plev;
int pnpkts, ppktlength;
double *c_in, *c_out, *d_out, *c_outR, *d_outR;

void wpsub(double *c_in, int lc_in, double *c_out, double *d_out,
	double *c_outR, double *d_outR, double *H, int *LengthH);


/*
 * i represents the child level. Go through each child level, filling in
 * coefficients as you go
 */

for(i=(int)*nlev-1; i>=(int)*finish_level; --i) {

    plev = i+1;     /* Parent level             */
    pnpkts = NPKTS((int)plev, *nlev);   /* Number of pkts at p lev */
    ppktlength = PKTLENGTH((int)plev);  /* Length at parent level  */

    /* Create input and output packets  */

    if ((c_in = (double *)malloc((unsigned)(ppktlength*sizeof(double))))==NULL) {
        *error = 1;
        return;
        }

    if ((c_out = (double *)malloc((unsigned)(sizeof(double)*ppktlength/2)))==NULL)  {
        *error = 2;
        return;
        }

    if ((d_out = (double *)malloc((unsigned)(sizeof(double)*ppktlength/2)))==NULL)  {
        *error = 3;
        return;
        }


    if ((c_outR = (double *)malloc((unsigned)(sizeof(double)*ppktlength/2)))==NULL) {
        *error = 4;
        return;
        }

    if ((d_outR = (double *)malloc((unsigned)(sizeof(double)*ppktlength/2)))==NULL) {
        *error = 5;
        return;
        }
        
    for(j=0; j< pnpkts; ++j)    {   /* Go thru each parent pkt */

        /* Copy parent packet to c_in */

        for(k=0; k<ppktlength; ++k)
            *(c_in+k) = ACCWPST(ansvec, plev, avixstart, j, k);

        /* Now compute answer using filter application */ 

        wpsub(c_in, (int)ppktlength, c_out, d_out, c_outR, d_outR,
            H, LengthH);


        /* Now store the answer packets */

        for(k=0; k < ppktlength/2; ++k) {
            ACCWPST(ansvec, i, avixstart, 4*j, k) = *(c_out+k);
            ACCWPST(ansvec, i, avixstart, 4*j+1, k) = *(d_out+k);
            ACCWPST(ansvec, i, avixstart, 4*j+2, k) = *(c_outR+k);
            ACCWPST(ansvec, i, avixstart, 4*j+3, k) = *(d_outR+k);
            }

        /* Next parent packet */
        }

    /* Free packet memory */
    free((void *)c_in);
    free((void *)c_out);
    free((void *)d_out);
    free((void *)c_outR);
    free((void *)d_outR);

    /* Next level */
    }
}

/*
 * WPSUB    -   Subroutine to compute basic step of stationary
 *          wavelet packet algorithm.
 */

void wpsub(double *c_in, int lc_in, double *c_out, double *d_out,
	double *c_outR, double *d_outR, double *H, int *LengthH)
/*---------------------
 * Argument description
 *---------------------
double *c_in::    Data input                       
int lc_in::   	  Length of input                  
double *c_out::   result of low pass applied to input          
double *d_out::   result of high pass applied to input         
double *c_outR::  result of low pass applied to rotated input      
double *d_outR::  result of high pass applied to rotated input     
double *H::       Wavelet filter                   
int *LengthH::    Length of wavelet filter             
 *---------------------*/
{
int lengthout;

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);
void rotater(double *book, int length);

lengthout = lc_in/2;

/* Obtain the father and mother wavelet coefficients of the data */

convolveC(c_in, lc_in, 0, H, (int)*LengthH,
        c_out, 0, (lengthout-1l),
        WAVELET, 1, PERIODIC);

convolveD(c_in, lc_in, 0, H, (int)*LengthH,
        d_out, 0, (lengthout-1l),
        WAVELET, 1, PERIODIC);

/* Obtain the father and mother wavelet coefficients of the rotated data */

rotater(c_in, lc_in);
        
convolveC(c_in, lc_in, 0, H, (int)*LengthH,
        c_outR, 0, (lengthout-1l),
        WAVELET, 1, PERIODIC);

convolveD(c_in, lc_in, 0, H, (int)*LengthH,
        d_outR, 0, (lengthout-1l),
        WAVELET, 1, PERIODIC);

/* That's it */

}


void accessDwpst(double *coefvec, int *lansvec, int *nlev, int *avixstart,
	int *primaryindex, int *nwppkt, int *pklength, int *level,
	double *weave, int *lweave, int *error)
/*---------------------
 * Argument description
 *---------------------
double *coefvec::     Vector storing the stat. wavelet pack. coefs 
int *lansvec::        Length of previous vector            
int *nlev::           The number of levels in the transform    
int *avixstart::      Index into coefvec for starting position
             	      of each level of coefficients        
int *primaryindex::   Packet numbers in the SWPT scheme to take out & interweave
int *nwppkt::         The number of ordinary WP packets at this level
             	      Also the length of the primaryindex vector
int *pklength::       The length of an SWPT packet         
int *level::          The level that we're extracting from     
double *weave::       A vector to store the answer         
int *lweave::         The length of the weave vector       
int *error::          An error code                
 *---------------------*/
{
register int counter, pklcount, i;

*error = 0;

counter = 0;


for(pklcount=0; pklcount< *pklength; ++pklcount)    {
    for(i=0; i<*nwppkt; ++i)    {   /* For each SWPT packet */

        *(weave+counter) = ACCWPST(coefvec, *level, avixstart,
            (*(primaryindex+i)), pklcount);
        ++counter;
        }
    }

}

/*
 * Use the binary representation of *l to build a number in base 4
 */

void c2to4(int *l, int *a)
{
register int ndigits,i ;
int mask;
int multiplier;

*a = 0;

if (*l == 0)
    return;

ndigits = (int)ceil( log((double)*l)/log(2.00));

ndigits++;  /* For exact powers of 2 */

mask = 1;
multiplier=1;

for(i=0; i<ndigits; ++i)    {
    *a += multiplier*((mask& *l)>>i);
    mask <<= 1;
    multiplier *= 4;
    }
    
}


/* Next code is from Arne Kovac */

/* Increasing this value would remove some nearly empty diagonals */

double thr=0.0;

/* The doubledouble structure is used by the makegrid function for
   sorting the data with respect to the x-component. */

struct doubledouble
           {
           double x;
           double y;
           } Doubledouble;

/* The ddcomp function is used by the makegrid function for sorting
   the data as an argument for qsort. */
           
int ddcomp(const void *a, const void *b)        /* MAN: changed to work in qsort below */
  {

struct doubledouble *q1=(struct doubledouble *)a; 
struct doubledouble *q2=(struct doubledouble *)b;  
  int t;
  
  if(q1->x>q2->x)
/*    return 1; */
t=1;
  else
  if(q1->x<q2->x)
/*    return -1;  */
t=-1;
  else
/*    return 0;   */
t=0;    
    return t;
  }

void makegrid(double *x, double *y, int *n, double *gridx, double *gridy,
               int *gridn, double *G, int *Gindex)

               /* This function computes from observations in x und y
                  new data on a grid of length gridn as well as a
                  description of the matrix that maps the original
                  data (or better the ordered original data, so that
                  x[i]<=x[j] when i<=j) to the new grid data. 
               
                  input:
                  x, y   vector of observations of length n 
                  gridn

                  output:
                  gridx, gridy the constructed grid of length gridn */
  {
  struct doubledouble *q;
  int i,li=0,ind;
  

  /* First, sort the data with respect to x. */

  q=(void *)malloc(*n*sizeof(Doubledouble));

  for(i=0;i<*n;i++)
    {
    q[i].x=x[i];
    q[i].y=y[i];
    }
  qsort(q,(size_t)*n,sizeof(struct doubledouble),ddcomp);   /* MAN: used to be (int *)ddcomp */

  /* Now create the new grid data. */

  for(i=0;i<*gridn;i++)
    {
    gridx[i]=(i+0.5)/(*gridn);

    /* Determine the index of the nearest observation left to the grid time 
       point. */

    while((li<*n-1)&&(q[li+1].x<gridx[i]))
      li++;
    if(li==(*n-1)) /* We are at the right end */
      {
      gridy[i]=q[li].y;
      ind=li-1;
      G[i]=0;
      }
    else
    if(q[li].x>=gridx[i]) /* We are at the left end */
      {
      gridy[i]=q[0].y;
      ind=0;
      G[i]=1;
      }
    else
      {
      gridy[i]=q[li].y+(gridx[i]-q[li].x)*(q[li+1].y-q[li].y)
                                                /(q[li+1].x-q[li].x); 
      ind=li;
      G[i]=1-(gridx[i]-q[li].x)/(q[li+1].x-q[li].x);
      }
    Gindex[i]=ind;
    }
  free(q);
  }

/* sigmastruct describes a covariance matrix of size n. diag is an n-vector
   of pointers to double vectors that correspond to the diagonals of the
   matrix. If diag[i]==NULL, then the i-th diagonal is empty. This
   representation is useful for covariance matrices with a band structure. */

struct sigmastruct 
          {
          int n;
          double **diag;
          };

/* createSigma allocates memory for a new covariance matrix of size n. */

int createSigma(struct sigmastruct *Sigma, int n)
  {
  int i;
  
  Sigma->n=n;
  if((Sigma->diag=malloc(n*sizeof(double *)))==NULL)
    {
    return(-1);
    }
  for(i=0;i<n;i++)
    (Sigma->diag)[i]=NULL;
  return(0);
  }
  
/* freeSigma releases the memory used by a sigmastruct element. */

void freeSigma(struct sigmastruct *Sigma)
  {
  int i;

  for(i=0;i<Sigma->n;i++)
    if(Sigma->diag[i]!=NULL)
      free((Sigma->diag)[i]);
  free(Sigma->diag);
  }

/* CleanupSigma removes diagonals that contain only elements < thr. */
   
void cleanupSigma(struct sigmastruct *Sigma)
  {
  int i,j;

  for(i=0;i<Sigma->n;i++)
    if((Sigma->diag)[i]!=NULL)
      {
      j=0;
      while((j<Sigma->n-i)&&(fabs((Sigma->diag)[i][j])<thr))
        j++;
      if(j>=Sigma->n-i)
        {
        free(Sigma->diag[i]);
        Sigma->diag[i]=NULL;
        }
      }
  }
  
/* putSigma changes the entry in the i-th row and j-th column to s and
   allocates memory for the diagonal if necessary. */

int putSigma(struct sigmastruct *Sigma, int i, int j, double s)
  {
  int d=abs(i-j);
  
  if(fabs(s)>0.0000001){        /* MAN: added brace to avoid ambiguity */
  if((i>=Sigma->n)||(j>=Sigma->n))
    {
    return(-1);
    /*
    puts("Error: This element does not exist.");
        */
    }
  else
    {
    if((Sigma->diag)[d]==NULL)
      if((Sigma->diag[d]=calloc(Sigma->n-d,sizeof(double)))==NULL)
        return(-2);
    (Sigma->diag)[d][(i+j-d)/2]=s;
    }
  }         /* MAN: added */
return(0);  /* MAN: added, hopefully won't have bad consequences. */
} 
  
/* allocateSigma allocates memory for diagonals of a covariance matrix,
   specified by the boolean vector d */
   

int allocateSigma(struct sigmastruct *Sigma, int *d)
{
int i;
  
for(i=0;i<Sigma->n;i++)
    if(d[i]==TRUE)
        if((Sigma->diag[i]=calloc(Sigma->n-i,sizeof(double)))==NULL)  {
            *d = (Sigma->n - i)*sizeof(double);
            return(-1);
            }

return(0);
}

/* computec is the function that computes the factors for the variances
   of the wavelet coefficients. 

   Gmatrix, Gindex describe the matrix that maps the original data
                   to the grid data as received by makegrid (s.o.)
   n, gridn        are the numbers of original and grid observations
   H, LengthH      describe the used wavelet filter. The filter coefficients
                   are stored in the vector H, their number in LengthH.
   bc              is either PERIODIC or SYMMETRIC, the boundary correction
                   to be used

   c contains afterwards the coefficients c_{jk}, such that
               Var(w_{jk})=c_{jk}\cdot\sigma^2
*/


void computec(int *n,double *c,int *gridn,double *Gmatrix,int *Gindex,
              double *H, int *LengthH, int *bc, int *error)
  {
  int virtgn,i,j,k,l,d,zaehler=0,gn=*gridn,laststart=0;
  int ii,dd,jj,o1,o2,iiG,iiH,jjG,jjH,gn2,LengthH2=*LengthH/2,dH,dG;
  int *NEEDIT,offset,offset2,band,band1,band2;
  double cellC,cellD,sig,G[20];
  double ProductG[20][20],ProductH[20][20];
  struct sigmastruct Sigma,Sigma2,Sigma3;
  int rc;

  int createSigma(struct sigmastruct *Sigma, int n);
  int putSigma(struct sigmastruct *Sigma, int i, int j, double s);
  int allocateSigma(struct sigmastruct *Sigma, int *d);

  if(*LengthH>20)
    {
    REprintf("Sorry, you can not use this procedure with more than 20 filter coefficients!");
    *error = 1;
    return; 
    }
  if((NEEDIT=malloc(*gridn*sizeof(int)))==NULL) {
    *error = 2;
    *n = *gridn * sizeof(int);  /* Contains number of bytes requested */
    return;
    }

  /* First step: Compute Filter G from Filter H */
  sig=-1.0;
  for(k=0;k<*LengthH;k++)
    {
    G[*LengthH-k-1]=sig*H[k];
    sig=-sig;
    }

  for(k=0;k<*LengthH;k++)
    for(l=0;l<*LengthH;l++)
      {
      ProductG[k][l]=G[k]*G[l];
      ProductH[k][l]=H[k]*H[l];
      }

  /* Second step: Compute the variance/covariance-matrix of the grid data */

  if (createSigma(&Sigma,gn) < 0)   {
        *error = 3;
        *n = (int)gn * sizeof(double);
        return;
        }
  for(i=0;i<*gridn;i++)
    {
    j=laststart;
    while(Gindex[i]-Gindex[j]>=2) 
      j++;
    laststart=j;
    for(;j<=i;j++)
      {
      switch(Gindex[i]-Gindex[j])
        {
        case  1: rc = putSigma(&Sigma,i,j,Gmatrix[i]*(1.0-Gmatrix[j]));

             if (rc < 0)    {
                if (rc == -1)   {
                    *error = 4;
                    return;
                    }

                if (rc == -2)   {
                    *error = 5;
                    *n = (Sigma.n - abs(i-j))*sizeof(double);
                    return;
                    }
                }
                    

                 break;
        case  0: rc= putSigma(&Sigma,i,j,Gmatrix[i]*Gmatrix[j]+(1.0-Gmatrix[i])*(1.0-Gmatrix[j]));
             if (rc < 0)    {
                if (rc == -1)   {
                    *error = 4;
                    return;
                    }

                if (rc == -2)   {
                    *error = 5;
                    *n = (Sigma.n - abs(i-j))*sizeof(double);
                    return;
                    }
                }
                 break;
        }
      }
    }
  
  /* And now the difficult part... */
  
  if(*bc==PERIODIC)
    {
    while(gn>=2) /* Apply the wavelet filters to the covariance matrix Sigma. */
      {
      gn2=gn/2; /* gn and gn2 are the sizes of Sigma and Sigma2 (or Sigma3). */ 

    /* Store the result of the high pass filter in sigma2... */
      if (createSigma(&Sigma2,gn2)<0){
        *error = 3l;
        *n = gn2 * sizeof(double);
        return;
        }

    /* ... and the result of the low pass filter in sigma3. */
      if (createSigma(&Sigma3,gn2)<0)   {
        *error = 3;
        *n = gn2 * sizeof(double);
        return;
        }
      cleanupSigma(&Sigma);       

      /* First we need to know which diagonals in sigma2 and sigma3 will not
         be empty. */

      band1=gn/2;
      band2=gn/2+1;
      while((band1>=0)&&(Sigma.diag[band1]==NULL))
	band1--; /* GPN Bugfix 22 Apr 2022. This was band1- which seems wrong. band1 needs to decrease itself */
      while((band2<=gn-1)&&(Sigma.diag[band2]==NULL))
        band2++;
      if(band1<=gn-band2)
        band=gn-band2;
      else
        band=band1;
        
      if(band+*LengthH>gn)
        for(d=0;d<gn2;d++)
          NEEDIT[d]=TRUE;
      else
        {
        for(d=0;d<gn2;d++)
          NEEDIT[d]=FALSE;
        for(d=0;(d<=(band+*LengthH)/2)&&(d<gn2);d++)
          NEEDIT[d]=TRUE;
        for(d=1;(d<=(band+*LengthH)/2)&&(d<gn2);d++)
          NEEDIT[gn2-d]=TRUE;
        }

      /* We allocate memory only for these diagonals. */

      if (allocateSigma(&Sigma2,NEEDIT) <0) {
        *error = 6;
        *n = *NEEDIT;
        return;
        }
        
      if (allocateSigma(&Sigma3,NEEDIT) < 0)    {
        *error = 6;
        *n = *NEEDIT;
        return;
        }

      /* Every entry of Sigma is involved in the computation of (LengthH/2)^2
         entries of Sigma2 and Sigma3. We find out in which and do the
         necessary computaions. */

      /* Let's start with the main diagonal elements of Sigma. */

      if(Sigma.diag[0]!=NULL)
        for(j=0;j<gn;j++)
          if(fabs(sig=Sigma.diag[0][j])>thr)
            {
            o1=j%2;
            iiH=j/2;
            iiG=(iiH+LengthH2-1)%(gn2);

            for(k=0;k<LengthH2;k++)
              {
              Sigma3.diag[0][iiH]+=sig*ProductH[2*k+o1][2*k+o1];
              Sigma2.diag[0][iiG]+=sig*ProductG[2*k+o1][2*k+o1];
              jjH=((j-2*(k+1)+*gridn)/2)%(gn2);
              jjG=(jjH+LengthH2-1)%(gn2);
              for(l=k+1;l<LengthH2;l++)
                { 
                dH=abs(iiH-jjH);
                dG=abs(iiG-jjG);
                if(dH==0)
                  Sigma3.diag[dH][iiH]+=2*sig*ProductH[2*k+o1][2*l+o1];
                else       
                  Sigma3.diag[dH][(iiH+jjH-dH)/2]+=sig*ProductH[2*k+o1][2*l+o1];
                if(dG==0)
                  Sigma2.diag[dG][iiG]+=2*sig*ProductG[2*k+o1][2*l+o1];
                else         
                  Sigma2.diag[dG][(iiG+jjG-dG)/2]+=sig*ProductG[2*k+o1][2*l+o1];
            
                if((--jjH)<0) jjH=gn2-1;
                if((--jjG)<0) jjG=gn2-1;
                }
              if((--iiH)<0) iiH=gn2-1;
              if((--iiG)<0) iiG=gn2-1;
              }
            }

      /* Now do the rest. */
      for(d=1;d<gn;d++)
        if(Sigma.diag[d]!=NULL)
          for(j=0;j<gn-d;j++)
            if(fabs(sig=Sigma.diag[d][j])>thr)
              {
              o1=(d+j)%2;
              iiH=((d+j+*gridn)/2)%(gn2);
              iiG=(iiH+LengthH2-1)%(gn2);

              for(k=0;k<LengthH2;k++)
                {
                jjH=((j+*gridn)/2)%(gn2);
                jjG=(jjH+LengthH2-1)%(gn2);
                o2=j%2;
                for(l=0;l<LengthH2;l++)
                  {
                  dH=abs(iiH-jjH);
                  dG=abs(iiG-jjG);
                  if(dH==0)
                    Sigma3.diag[dH][iiH]+=2*sig*ProductH[o1][o2];
                  else       
                    Sigma3.diag[dH][(iiH+jjH-dH)/2]+=sig*ProductH[o1][o2];
                  if(dG==0)
                    Sigma2.diag[dG][iiG]+=2*sig*ProductG[o1][o2];
                  else         
                    Sigma2.diag[dG][(iiG+jjG-dG)/2]+=sig*ProductG[o1][o2];
                  if((--jjH)<0) jjH=gn2-1;
                  if((--jjG)<0) jjG=gn2-1;
                  o2+=2;
                  }
                if((--iiH)<0) iiH=gn2-1;
                if((--iiG)<0) iiG=gn2-1;
                o1+=2;
                }
              }

      /* Now Sigma2 und Sigma3 should contain the right values. Store
         only the main diagonal of Sigma2 in c. */
      for(j=0;j<gn2;j++)
        c[zaehler++]=Sigma2.diag[0][j];

      /* Sigma and Sigma2 are now redundant. */
      freeSigma(&Sigma);
      freeSigma(&Sigma2);

      /* Sigma3 becomes our new Sigma. */
      memcpy(&Sigma,&Sigma3,sizeof(Sigma));
      gn=gn2;
      }
    } /* end of the periodic case */
  else
  if(*bc==SYMMETRIC) /* The symmetric case is even more difficult... */
    {
    virtgn=gn; /* virtgn is used to determine when we can finish and takes
                  exactly the same values as gn in the periodic case. */
    offset=0; /* Symmetric boundary conditions provide extra coefficients
                 on both sides. offset gives the number of such extra
                 coeffs on the left-hand side. */
    while(virtgn>=2)
      {
      /* First of all, we want to know how many diagonals and extra coeffs
         Sigma2 and Sigma3 have. */

      if(offset%2==0)
        gn2=(gn+1)/2+LengthH2-1;
      else
        gn2=(gn+2)/2+LengthH2-1;
      offset2=(offset+*LengthH-1)/2;

      /* Now allocate memory for them. */
      
      if (createSigma(&Sigma2,gn2) <0)  {
        *error = 3;
        *n = gn2 * sizeof(double);
        return;
        }
        
      if (createSigma(&Sigma3,gn2) <0)  {
        *error = 3;
        *n = gn2*sizeof(double);
        return;
        }
      cleanupSigma(&Sigma);
      
      /* Again, we need to know which diagonals in sigma2 and sigma3 will not
         be empty. */

      band=gn-1;
      while((band>=0)&&(Sigma.diag[band]==NULL))
        band--;
      if(band+2*(*LengthH)>gn)
        for(d=0;d<gn2;d++)
          NEEDIT[d]=TRUE;
      else
        {
        for(d=0;d<gn2;d++)
          NEEDIT[d]=FALSE;
        for(d=0;(d<=(band+*LengthH)/2)&&(d<gn2);d++)
          NEEDIT[d]=TRUE;
        }

      /* We allocate memory only for these diagonals. */

      if (allocateSigma(&Sigma2,NEEDIT) < 0)    {
        *error = 6;
        *n = *NEEDIT;
        return;
        }
      if (allocateSigma(&Sigma3,NEEDIT) < 0)    {
        *error = 6;
        *n = *NEEDIT;
        return;
        }

      /* Every entry of Sigma is involved in the computation of (LengthH/2)^2
         entries of Sigma2 and Sigma3. We find out in which and do the
         necessary computaions. */

      /* Let's start with the values that are not effected by boundary 
         correction.  First the main diagonal. */

      if(Sigma.diag[0]!=NULL)
        for(j=0;j<gn;j++)
          if(fabs(sig=Sigma.diag[0][j])>thr)
            {
            o1=(j+offset)%2;
            iiH=(j+offset%2)/2+LengthH2-1;
            iiG=iiH;

            for(k=0;k<LengthH2;k++)
              {
              Sigma3.diag[0][iiH]+=sig*ProductH[2*k+o1][2*k+o1];
              Sigma2.diag[0][iiG]+=sig*ProductG[2*k+o1][2*k+o1];
              jjH=(j+offset%2)/2+LengthH2-1-(k+1);
              jjG=jjH;
              for(l=k+1;l<LengthH2;l++)
                { 
                if((iiH<0)||(jjH<0))
                  REprintf("computec: ERROR in C code");
                dH=abs(iiH-jjH);
                dG=abs(iiG-jjG);
                if(dH==0)
                  Sigma3.diag[dH][iiH]+=2*sig*ProductH[2*k+o1][2*l+o1];
                else       
                  Sigma3.diag[dH][(iiH+jjH-dH)/2]+=sig*ProductH[2*k+o1][2*l+o1];
                if(dG==0)
                  Sigma2.diag[dG][iiG]+=2*sig*ProductG[2*k+o1][2*l+o1];
                else       
                  Sigma2.diag[dG][(iiG+jjG-dG)/2]+=sig*ProductG[2*k+o1][2*l+o1];
                jjH--;
                jjG--;
                }
              iiH--;
              iiG--;
              }
            }

      /* Now do the not-effected values on the other diagonals. */

      for(d=1;d<gn;d++)
        if(Sigma.diag[d]!=NULL)
          /*for(j=*LengthH-2+offset%2;j<gn-d-(*LengthH-2+offset%2);j++)*/
          for(j=0;j<gn-d;j++)
            if(fabs(sig=Sigma.diag[d][j])>thr)
              {
              o1=(d+j+offset)%2;
              iiH=(d+j+offset%2)/2+LengthH2-1;
              iiG=iiH;

              for(k=0;k<LengthH2;k++)
                {
                jjH=(j+offset%2)/2+LengthH2-1;
                jjG=jjH;
                o2=(j+offset)%2;
                for(l=0;l<LengthH2;l++)
                  {
                  dH=abs(iiH-jjH);
                  dG=abs(iiG-jjG);
                  if(dH==0)
                    Sigma3.diag[dH][iiH]+=2*sig*ProductH[o1][o2];
                  else       
                    Sigma3.diag[dH][(iiH+jjH-dH)/2]+=sig*ProductH[o1][o2];
                  if(dG==0)
                    Sigma2.diag[dG][iiG]+=2*sig*ProductG[o1][o2];
                  else         
                    Sigma2.diag[dG][(iiG+jjG-dG)/2]+=sig*ProductG[o1][o2];
                  jjH--;
                  jjG--;
                  o2+=2;
                  }
                iiH--;
                iiG--;
                o1+=2;
                }
              }

      /* Finally, do the boundary entries. */

      for(d=0;d<gn2;d++)
        if(Sigma2.diag[d]!=NULL)
          {
          for(j=0;(j<gn2-d)&&(j<LengthH2-1+offset%2);j++)
            {
            cellC=0;
            cellD=0;
            i=d+j;
            iiG=2*i-(*LengthH-2)-offset%2;
            for(k=0;k<*LengthH;k++,iiG++)
              {
              jjG=2*j-(*LengthH-2)-offset%2;
              for(l=0;l<*LengthH;l++,jjG++)
                {
                ii=iiG;
                jj=jjG;
                if(ii<0)
                  ii=-ii-1;
                if(jj<0)
                  jj=-jj-1;
                if(ii>=gn)
                  {
                  ii=2*gn-ii-1;
                  }
                if(jj>=gn)
                  {
                  jj=2*gn-jj-1;
                  }
                dd=abs(ii-jj);
                if(Sigma.diag[dd]!=NULL)
                  {
                  sig=(Sigma.diag)[dd][(ii+jj-dd)/2];
                  cellC+=sig*ProductH[k][l];
                  cellD+=sig*ProductG[k][l];
                  }
                }
              }
            Sigma2.diag[d][j]=cellD;
            Sigma3.diag[d][j]=cellC;
            }
          for(j=gn2-d-1;(j>=0)&&(j>=gn2-d-LengthH2+1-offset%2);j--)
            {
            cellC=0;
            cellD=0;
            i=d+j;
            iiG=2*i-(*LengthH-2)-offset%2;
            for(k=0;k<*LengthH;k++,iiG++)
              {
              jjG=2*j-(*LengthH-2)-offset%2;
              for(l=0;l<*LengthH;l++,jjG++)
                {
                ii=iiG;
                jj=jjG;
                if(ii<0)
                  ii=-ii-1;
                if(jj<0)
                  jj=-jj-1;
                if(ii>=gn)
                  {
                  ii=2*gn-ii-1;
                  }
                if(jj>=gn)
                  {
                  jj=2*gn-jj-1;
                  }
                dd=abs(ii-jj);
                if(Sigma.diag[dd]!=NULL)
                  {
                  sig=(Sigma.diag)[dd][(ii+jj-dd)/2];
                  cellC+=sig*ProductH[k][l];
                  cellD+=sig*ProductG[k][l];
                  }
                }
              }
            Sigma2.diag[d][j]=cellD;
            Sigma3.diag[d][j]=cellC;
            }
          }

      /* This looks now pretty the same as in the periodic case. */

      for(j=0;j<gn2;j++)
        c[zaehler++]=Sigma2.diag[0][j];
      freeSigma(&Sigma);
      freeSigma(&Sigma2);
      memcpy(&Sigma,&Sigma3,sizeof(Sigma));

      gn=gn2;
      offset=offset2;
      virtgn=virtgn/2;
      }
    } /* end of the symmetric case */

  /* We are finished with the computation of c and can release any memory. */

  freeSigma(&Sigma);
  free(NEEDIT);
  }

/* Back to GPN code */

/*
 * Error codes for the rainmat suite.
 *
 * 0 - No error
 *
 * Errors here:
 *
 * rainmat: 100 memory error for creation of w
 *      101 memory error for creation of subvectors of w
 *          J returns the number of the subvector
 *          that couldn't be created. 
 *
 * wlpart   110 memory error for creation of TheData 
 *      111 memory error for creation of ixvec
 *
 * rotateleft   120 memory error for creation of tmp
 *
 * rainmatPARENT 130    memory error for creation of lvec
 *
 * mkcoef   140 memory error for creation of ixvec
 *      141 memory error for creation of TheData
 *      142 memory error for creation of lcoefvec
 *      143 memory error for creation of tmpcfvec
 *
 * rainmatPARTIAL 150   memory error for creation of lvec
 *
 * PsiJonly:    160 wout is not int enough to store answer
 *          (change a constant in the S calling function,
 *          the variable lwout contains the number of doubles
 *          required)
 *      161 memory error for creation of w
 *      162 memory error for creation of subvectors of w
 *          J returns the number of the subvector
 *          that couldn't be created. 
 *
 * 
 *
 * Other error codes can occur from other parts of WaveThresh
 *
 */
        

/*
 * Compute the matrix that obtains S from I
 *
 */

void rainmat(int *J, int *donej, double **coefvec, int *lvec,
	double *fmat, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::	     The desired maximum level (positive)     
int *donej::         The first j columns already filled       
double **coefvec::   The \psi_{jk} stacked into one vector    
int *lvec::          A vector of lengths of each \psi_j vector in coefvec.
		     The jth element is the length of the jth \psi_j in coefvec
double *fmat::       This vector will contain the answer. This is
                     the lower triangular portion of the J*J matrix,
                     and therefore is of length J(J-1)/2 
int *error::         Error code                   
 *---------------------*/
{

/* First we compute the w. One for each j           */

double **w;
register int j,k,m,l;
double sum;
int lj,ll;


if ((w = (double **)malloc((unsigned)*J*sizeof(double *)))==NULL)   {
    *error = 100;
    return;
    }

/* Now populate each of the *w */

for(j=0; j<*J; ++j) {
    if ((*(w+j) = (double *)malloc((unsigned)(*(lvec+j)*2-1)*sizeof(double)))==NULL)    {
        *error = 101;
        *J = (int)j;
        return;
        }
    }

/* Now compute each of the wjk */

for(j=0; j< *J; ++j)    {
    lj = *(lvec+j);
    for(k = 1-lj; k <= lj-1; ++k)   {
        sum = 0.0;
        for(m = max(0, k); m <= min(lj-1, lj-1+k); ++m) {
            sum += *((*(coefvec+j))+m) *
                *((*(coefvec+j))+m-k);
            }
        ACCESSW(w, j, k-1+lj) = sum;
        }
    }

/* Now compute the F */

for(j=0; j<*J; ++j) {
    lj = *(lvec+j);
    for(l=j; l<*J; ++l) {
        if (l >= *donej)    {
            ll = *(lvec+l);
            sum = 0.0;
            for(k=max(1-ll, 1-lj); k <= min(lj-1, ll-1); ++k) {
                sum += ACCESSW(w, j, k-1+lj) *
                    ACCESSW(w, l, (-k)-1+ll);
                }
            *(fmat+*J*l+j) = *(fmat+*J*j+l) = sum;
            }
        }
    }


/* Now free the w */
for(j=0; j<*J; ++j) {
    free((void *)*(w+j));
    }
free((void *)w);
}



void wlpart(int *J, int *BigJ, double *H, int *LengthH, int *error)
{
register int KeepGoing;
register int somefull;
register int allnonzero;
register int i;
register int j;
double *TheData;
int ndata;
double *C, *D;
int *firstC, *lastC, *offsetC, *firstD, *lastD, *offsetD;
int LengthC, LengthD, levels;
int type,bc;
int *ixvec;


void simpleWT(double *TheData, int *ndata, double *H, int *LengthH,
    double **C, int *LengthC, double **D, int *LengthD, int *levels,
    int **firstC, int **lastC, int **offsetC,
    int **firstD, int **lastD, int **offsetD,
    int *type, int *bc, int *error);
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);


*error=0;

*BigJ = *J + 1;
KeepGoing = TRUE;

while(KeepGoing)    {

    /* Rprintf("Entered loop BigJ is %ld\n", *BigJ); */

    ndata = (int)0x01 << *BigJ;

    /* Rprintf("ndata is %ld\n", ndata);*/

    /*
     * Basically a dummy wavelet transform to set up first/last stuff
     */
    if ((TheData = (double *)malloc((unsigned)ndata*sizeof(double)))==NULL) {
        *error = 110;
        return;
        }

    for(i=0; i<ndata; ++i)
        *(TheData+i) = 0.0;

    /*
     * Do the wavelet transform
     */
    
    simpleWT(TheData, &ndata, H, LengthH,
        &C, &LengthC, &D, &LengthD, &levels,
        &firstC, &lastC, &offsetC,
        &firstD, &lastD, &offsetD,
        &type, &bc, error);

    if (*error != 0)
        return;

    /*
     * Create ixvec
     */

    if ((ixvec = (int *)malloc((unsigned)*BigJ*sizeof(int)))==NULL){
        *error = 111;
        return;
        }

    for(i=0; i< *BigJ; ++i)
        *(ixvec+i) =(0x01 << (*BigJ -1 - i)); 

    for(i=1; i< *BigJ; ++i)
        *(ixvec+i) = *(ixvec+i-1) + *(ixvec+i);

    for(i=0; i< *BigJ; ++i)
        --*(ixvec+i);

    somefull = FALSE;

    for(i=0; i<*J; ++i) {

        for(j=0; j<LengthD; ++j)
            *(D+j) = 0;

        *(D+ *(ixvec+i)) = 1;

        waverecons(C, D, H, LengthH, &levels,
            firstC, lastC, offsetC, firstD, lastD, offsetD,
            &type, &bc, error);

        if (*error != 0)
            return;

        allnonzero = TRUE;

        for(j=0; j<ndata; ++j)  {
            if (*(C+j) == 0.0)  {
                allnonzero = FALSE;
                break;
                }
            }
        somefull = somefull || allnonzero;


        if (somefull)
            break;

    }
    /* Free memory */

    free((void *)C);
    free((void *)D);
    free((void *)firstC); free((void *)lastC);
    free((void *)offsetC);
    free((void *)firstD); free((void *)lastD);
    free((void *)offsetD);
    free((void *)ixvec);
    free((void *)TheData);

    if (!somefull)
        KeepGoing = FALSE;
    else
        *BigJ = *BigJ + 1;
    }
}

/*
 * Return the index number of the last zero in a vector
 * Return -1 if no such exists
 */

int idlastzero(double *v, int *nv)
{
register int i;

for(i= *nv-1; i>=0; --i)
        if (*(v+i) == 0.0)
                break;

return(i);
}

/*
 * Cyclically rotate a vector n places to the left
 * (a C replacement for guyrot)
 */

void rotateleft(double *v, int *nv, int *n, int *error)
{
register int i;
double *tmp;    /* Storage for the ones the fall off the left   */

*error = 0;

*n = *n % *nv;

if (*n == 0)    /* No rotation required */
    return;

if ((tmp = (double *)malloc((unsigned)(*n)*sizeof(double))) == NULL)    {
    *error = 120;
    return;
    }

for(i=0; i< *n; ++i)    
    *(tmp+i) = *(v+i);

for(i=0; i< *nv - *n; ++i)
    *(v+i) = *(v+i+*n);

for(i=0; i< *n; ++i)
    *(v+ i + *nv - *n) = *(tmp+i);
    
free((void *)tmp);
}

void rainmatPARENT(int *J, double *H, int *LengthH, double *fmat,
	double *tol, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::     	 The dimension of the problem             
double *H::      The wavelet filter coefficients          
int *LengthH::   The number of wavelet filter coefficients        
double *fmat::   The answer                       
double *tol::    Elements smaller than this will be deleted       
int *error::     Error code. Nonzero is an error          
 *---------------------*/
{
register int i;
int BigJ;  /* The level we must go to to be able to compute
         * coefficients without error
         */
int donej; /* Only for partial                 */

double **coefvec;   /* These are the \Psi_j (\tau)          */
int *lvec;     /* Vector of length *J contains the length  
             * of each vector in coefvec
             */

void wlpart(int *J, int *BigJ, double *H, int *LengthH, int *error);
void mkcoef(int *J, int BigJ, double *H, int *LengthH, double ***coefvec,
	int *lvec, double *tol, int *error);
void rainmat(int *J, int *donej, double **coefvec, int *lvec,
	double *fmat, int *error);
void haarmat(int *J, int *donej, double *fmat, int *error);

donej = 0;

if (*LengthH == 2)  /* Haar - can compute exactly */
    {
    haarmat(J, &donej, fmat, error);
    return;
    }

/* whichlevel */

wlpart(J, &BigJ, H, LengthH, error);

if (*error != 0)
    return;

/* mkcoef */

if ((lvec = (int *)malloc((unsigned)*J*sizeof(int)))==NULL)   {
    *error = 130;
    return;
    }

for(i=0; i<*J; ++i)
    *(lvec+i) = 0;

mkcoef(J, BigJ, H, LengthH, &coefvec, lvec, tol, error); 

if (*error != 0)
    return;

/* rainmat */ 
donej = 0;

rainmat(J, &donej, coefvec, lvec, fmat, error);

if (*error != 0)
    return;

free((void *)lvec);

for(i=0; i<*J; ++i)
    free((void *)*(coefvec+i));

free((void *)coefvec);
}

/* Make \Psi_j(\tau) components */

void mkcoef(int *J, int BigJ, double *H, int *LengthH, double ***coefvec,
	int *lvec, double *tol, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::             Dimension of the problem             
int BigJ::           The maximum depth that we have to go to      
double *H::          Wavelet filter coefficients              
int *LengthH::       Number of wavelet filter coefficients        
double ***coefvec::  Coefficients of \Psi_j(\tau)          
int *lvec::          Vector of length *J that will contain length of each
		     component of coefvec 
double *tol::        Elements smaller than this will be deleted       
int *error::         Error code                       
 *---------------------*/
{
register int i,j;
register int large_ones;
int ndata;
int *ixvec;        /* Index vector for inserting 1s into blank WT  */
double **lcoefvec;  /* Local version of coefvec         */
double *tmpcfvec;   /* Temporary vector             */

/* Things needed for the simpleWT */
double *TheData;
double *C, *D;
int *firstC, *lastC, *offsetC, *firstD, *lastD, *offsetD;
int LengthC, LengthD, levels;
int type,bc;
int n_to_rotate;


void simpleWT(double *TheData, int *ndata, double *H, int *LengthH,
    double **C, int *LengthC, double **D, int *LengthD, int *levels,
    int **firstC, int **lastC, int **offsetC,
    int **firstD, int **lastD, int **offsetD,
    int *type, int *bc, int *error);
int idlastzero(double *v, int *nv);
void rotateleft(double *v, int *nv, int *n, int *error);
void waverecons(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error);

ndata = (int)0x01 << BigJ;

/*
 * Create ixvec
 */

if ((ixvec = (int *)malloc((unsigned)BigJ*sizeof(int)))==NULL){
    *error = 140;
    return;
    }

for(i=0; i< BigJ; ++i)
    *(ixvec+i) =(0x01 << (BigJ -1 - i)); 

for(i=1; i< BigJ; ++i)
    *(ixvec+i) = *(ixvec+i-1) + *(ixvec+i);

for(i=0; i< BigJ; ++i)
    --*(ixvec+i);

/*
 * Basically a dummy wavelet transform to set up first/last stuff
 */
if ((TheData = (double *)malloc((unsigned)ndata*sizeof(double)))==NULL) {
    *error = 141;
    return;
    }

for(i=0; i<ndata; ++i)
    *(TheData+i) = 0.0;

/*
 * Do the wavelet transform
 */

simpleWT(TheData, &ndata, H, LengthH,
    &C, &LengthC, &D, &LengthD, &levels,
    &firstC, &lastC, &offsetC,
    &firstD, &lastD, &offsetD,
    &type, &bc, error);

if (*error != 0)
    return;

if ((lcoefvec = (double **)malloc((unsigned)*J*sizeof(double *)))==NULL){
    *error = 142;
    return;
    }

for(i=1; i<= *J; ++i)   {
    for(j=0; j<LengthD; ++j)
        *(D+j) = 0.0;

    *(D+ *(ixvec+i-1)) = 1;

    waverecons(C, D, H, LengthH, &levels,
        firstC, lastC, offsetC, firstD, lastD, offsetD,
        &type, &bc, error);

    if (*error != 0)
        return;

    /* Now copy reconstruction into TheData  (vec in S) */

    for(j=0; j<ndata; ++j)
        *(TheData+j) = *(C+j);

    n_to_rotate = (int)idlastzero(TheData, &ndata);

    if (n_to_rotate < 0)
        n_to_rotate = 0;

    rotateleft(TheData, &ndata, &n_to_rotate, error);

    if (*error != 0)
        return;

    large_ones = 0;

    for(j=0; j<ndata; ++j)
        if (fabs(*(TheData+j)) > *tol)
            ++large_ones;

    /* Now get memory for the large ones */

    if ((tmpcfvec = (double *)malloc((unsigned)large_ones*sizeof(double)))==NULL)   {
        *error = 143;
        return;
        }

    large_ones = 0;

    for(j=0; j<ndata; ++j)
        if (fabs(*(TheData+j)) > *tol)
            *(tmpcfvec+large_ones++) = *(TheData+j);


    /* Install this vector into the array */

    *(lcoefvec+i-1) = tmpcfvec;
    *(lvec+i-1) = (int)large_ones;
    }

/* Install the lcoefvec into the coefvec */

*coefvec = lcoefvec;

free((void *)ixvec);
free((void *)TheData);
}

void rainmatOLD(int *J, double *coefvec, int *ixvec, int *lvec,
	double *fmat, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::           The desired maximum level (positive)     
double *coefvec::  The \psi_{jk} stacked into one vector    
int *ixvec::
int *lvec::        A vector of lengths of each \psi_j vector in coefvec.
		   The jth element is the length of the jth \psi_j in coefvec
double *fmat::     This vector will contain the answer. This is the lower
		   triangular portion of the J*J matrix, and therefore is of
		   length J(J-1)/2 
int *error::       Error code
                	1-  Generating **w
                	2+j Memory error on 2+j th one
 *---------------------*/
{

/* First we compute the w. One for each j           */

double **w;
register int j,k,m,l,cnt;
double sum;
int lj,ll;

if ((w = (double **)malloc((unsigned)*J*sizeof(double *)))==NULL)   {
    *error = 1;
    return;
    }

/* Now populate each of the *w */

for(j=0; j<*J; ++j) {
    if ((*(w+j) = (double *)malloc((unsigned)(*(lvec+j)*2-1)*sizeof(double)))==NULL)    {
        *error = (int)(2+j);
        return;
        }
    }

/* Now compute each of the wjk */

for(j=0; j< *J; ++j)    {
    lj = *(lvec+j);
    for(k = 1-lj; k <= lj-1; ++k)   {
        sum = 0.0;
        for(m = max(0, k); m <= min(lj-1, lj-1+k); ++m) {
            sum += *(coefvec+*(ixvec+j)+m) *
                *(coefvec+*(ixvec+j)+m-k);
            }
        ACCESSW(w, j, k-1+lj) = sum;
        }
    }

/* Now compute the F */

cnt = 0;
for(j=0; j<*J; ++j) {
    lj = *(lvec+j);
    for(l=j; l<*J; ++l) {
        ll = *(lvec+l);
        sum = 0.0;
        for(k=max(1-ll, 1-lj); k <= min(lj-1, ll-1); ++k)   {
            sum += ACCESSW(w, j, k-1+lj) *
                ACCESSW(w, l, (-k)-1+ll);
            }
        
        *(fmat+*J*l+j) = *(fmat+*J*j+l) = sum;
        ++cnt;
        }
    }

/* Now free the w */
for(j=0; j<*J; ++j) {
    free((void *)*(w+j));
    }
free((void *)w);
}

/* rainmatPARTIAL - partial matrix filling              */

void rainmatPARTIAL(int *J, int *donej, double *H, int *LengthH,
	double *fmat, double *tol, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::         The dimension of the problem             
int *donej::     The first j dimensions are already filled        
double *H::      The wavelet filter coefficients          
int *LengthH::   The number of wavelet filter coefficients        
double *fmat::   The answer                       
double *tol::    Elements smaller than this will be deleted       
int *error::     Error code. Nonzero is an error          
 *---------------------*/
{
register int i;
int BigJ;  /* The level we must go to to be able to compute
         * coefficients without error
         */

double **coefvec;   /* These are the \Psi_j (\tau)          */
int *lvec;     /* Vector of length *J contains the length  
             * of each vector in coefvec
             */

void wlpart(int *J, int *BigJ, double *H, int *LengthH, int *error);
void mkcoef(int *J, int BigJ, double *H, int *LengthH, double ***coefvec,
	int *lvec, double *tol, int *error);
void rainmat(int *J, int *donej, double **coefvec, int *lvec,
	double *fmat, int *error);
void haarmat(int *J, int *donej, double *fmat, int *error);

if (*LengthH == 2)  /* Haar - can compute exactly */
    {
    haarmat(J, donej, fmat, error);
    return;
    }

/* whichlevel */

wlpart(J, &BigJ, H, LengthH, error);

if (*error != 0)
    return;

/* mkcoef */

if ((lvec = (int *)malloc((unsigned)*J*sizeof(int)))==NULL)   {
    *error = 150;
    return;
    }

for(i=0; i<*J; ++i)
    *(lvec+i) = 0;

mkcoef(J, BigJ, H, LengthH, &coefvec, lvec, tol, error); 

if (*error != 0)
    return;

/* rainmat */ 

rainmat(J, donej, coefvec, lvec, fmat, error);

if (*error != 0)
    return;

free((void *)lvec);

for(i=0; i<*J; ++i)
    free((void *)*(coefvec+i));

free((void *)coefvec);
}

void PsiJ(int *J, double *H, int *LengthH, double *tol, double *wout,
	int *lwout, int *rlvec, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::          The dimension of the problem             
double *H::       The wavelet filter coefficients          
int *LengthH::    The number of wavelet filter coefficients        
double *tol::     Elements smaller than this will be deleted       
double *wout::    Answers for \Psi_j(\tau)             
int *lwout::      Length of previous array             
int *rlvec::      Vector of length J contains lengths of \psi_j    
int *error::      Error code. Nonzero is an error          
 *---------------------*/
{
register int i;
int BigJ;  /* The level we must go to to be able to compute
         * coefficients without error
         */
double **coefvec;   /* These are the \psi_j (\tau)          */
int *lvec;     /* Vector of length *J contains the length  
             * of each vector in coefvec
             */

void wlpart(int *J, int *BigJ, double *H, int *LengthH, int *error);
void mkcoef(int *J, int BigJ, double *H, int *LengthH, double ***coefvec,
	int *lvec, double *tol, int *error);
void PsiJonly(int *J, double **coefvec, int *lvec, double *wout,
	int *lwout, int *error);

/* whichlevel */

wlpart(J, &BigJ, H, LengthH, error);

if (*error != 0)
    return;

/* mkcoef */

if ((lvec = (int *)malloc((unsigned)*J*sizeof(int)))==NULL)   {
    *error = 130;
    return;
    }

for(i=0; i<*J; ++i)
    *(lvec+i) = 0;

mkcoef(J, BigJ, H, LengthH, &coefvec, lvec, tol, error); 

if (*error != 0)
    return;


PsiJonly(J, coefvec, lvec, wout, lwout, error);

if (*error != 0)
    return;

for(i=0; i<*J; ++i)
    *(rlvec + i) = *(lvec+i);

free((void *)lvec);

for(i=0; i<*J; ++i)
    free((void *)*(coefvec+i));

free((void *)coefvec);
}

void PsiJonly(int *J, double **coefvec, int *lvec, double *wout,
	int *lwout, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::              The desired maximum level (positive)     
double **coefvec::    The \psi_{jk} stacked into one vector    
int *lvec::           A vector of lengths of each \psi_j vector in
                      coefvec. The jth element is the length of the
                      jth \psi_j in coefvec
double *wout::        Output contains the \Psi_j(\tau)     
int *lwout::          Length of this vector. If it is not long enough an
		      error code is returned     
int *error::          Error code                   
 *---------------------*/
{

/* First we compute the w. One for each j           */

double **w;
register int j,k,m;
double sum;
int totall;
int lj,cnt;

/* Check output vector is long enough to store answer */

totall = 0;
for(j=0; j < *J; ++j)
    totall += *(lvec+j)*2l - 1l; 

if (totall > *lwout)    {
    *error = 160;
    *lwout = totall;
    return;
    }


if ((w = (double **)malloc((unsigned)*J*sizeof(double *)))==NULL)   {
    *error = 161;
    return;
    }

/* Now populate each of the *w */

for(j=0; j<*J; ++j) {
    if ((*(w+j) = (double *)malloc((unsigned)(*(lvec+j)*2-1)*sizeof(double)))==NULL)    {
        *error = 162;
        *J = (int)j;
        return;
        }
    }

/* Now compute each of the wjk */

for(j=0; j< *J; ++j)    {
    lj = *(lvec+j);
    for(k = 1-lj; k <= lj-1; ++k)   {
        sum = 0.0;
        for(m = max(0, k); m <= min(lj-1, lj-1+k); ++m) {
            sum += *((*(coefvec+j))+m) *
                *((*(coefvec+j))+m-k);
            }
        ACCESSW(w, j, k-1+lj) = sum;
        }
    }

/* Store the w */

cnt = 0;

for(j=0; j < *J; ++j)   {
    lj = *(lvec+j);
    for(k = 1-lj; k <= lj-1; ++k)   {
        *(wout+cnt) = ACCESSW(w, j, k-1+lj);
        ++cnt;
        }
    }


/* Now free the w */
for(j=0; j<*J; ++j) {
    free((void *)*(w+j));
    }
free((void *)w);
}

/* haarmat - Computes matrix exactly using formula      */

void haarmat(int *J, int *donej, double *fmat, int *error)
/*---------------------
 * Argument description
 *---------------------
int *J::         The desired maximum level (positive)     
int *donej::     The first j columns already filled       
double *fmat::   This vector will contain the answer. This is the lower
		 triangular portion of the J*J matrix, and therefore is of
 		 length J(J-1)/2 
int *error::     Error code                   
 *---------------------*/
{
register int j,l;
double a;
double twoj, twol, two2j, two2jmo;

for(j=0; j<*J; ++j) {
    for(l=j; l<*J; ++l) {
        if (l >= *donej)    {

            if (l==j)   {
                twoj = pow(2.0, ((double)j+1)); 
                two2j = twoj*twoj;
                a = (two2j + 5.0)/(3.0*twoj);
                }

            else    {
                two2jmo = pow(2.0, (double)(2*j+1));
                twol = pow(2.0, ((double)l+1));
                a = (two2jmo + 1.0)/twol;
                }
            *(fmat+*J*l+j) = *(fmat+*J*j+l) = a;
            }
        }
    }

}

/*
 * Now follows the code from swt2d.c
 */


/*
 * Perform whole of SWT2D after initialising
 */

void SWT2Dall(double *m, int *nm, double *am, int *J, double *H,
	int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *m::    The input data                   
int *nm::      The dimension of the square matrix m         
double *am::   The *big* answer 3D array                
int *J::       The level at which to store the initial information  
double *H::    The smoothing filter             
int *LengthH:: The length of the smoothing filter       
int *error::   Error code 0=ok, anything else is memory error   
 *---------------------*/
{
int D1, D12;    /* Dimensions of am array               */
int nm2, nm4;   /* nm divided by 2 then 4               */
void initSWT2D(double *m, int *nm, double *am, int *J, double *H,
	int *LengthH, int *error);
void SWT2Drec(double *am, int D1, int D12, int x, int y, int TWOsl,
	int sl, int J, double *H, int *LengthH, int *error);

*error = 0;

initSWT2D(m, nm, am, J, H, LengthH, error);

if (*error != 0)
    return;

/*
 * Now for each level use the previous level as the coefficients to
 * do a 2D wavelet transform for the next level
 *
 * Produce level J-2 from J-1 (which was done in initSWT2D)
 * Produce level J-3 from J-2 ...
 * ...
 * Produce level 0 from 1
 * Go home!
 */
D12 = (*J)*(*nm * 2);
D1 = (*J);

nm2 = *nm/2;
nm4 = nm2/2;

SWT2Drec(am, D1, D12, 0, 0, nm2, nm4, *J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, *nm, 0, nm2, nm4, *J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, 0, *nm, nm2, nm4, *J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, *nm, *nm, nm2, nm4, *J-1, H, LengthH, error);
if (*error != 0)
    return;
}

void SmallStore(double *am, int D1, int D12, int J, int sl, int x, int y,
	int ix, int jy,
	double *hhout, double *hgout, double *ghout, double *ggout, int nm)
/*---------------------
 * Argument description
 *---------------------
double *am::      The *big* matrix to store everything in      
int D1::          First dimension of am                
int D12::         First and second dimensions of am multiplied     
int J::           The level to fill                    
int sl::          Side length of small packets             
int x::           The origin x coordinate              
int y::           The origin y coordinate              
int ix::          The smaller matrix i offset              
int jy::          The smaller matrix j offset              
double *hhout::   The new smoothed matrix              
double *hgout::   The new horizontal detail matrix         
double *ghout::   The new vertical detail matrix           
double *ggout::   The new diagonal detail matrix           
int nm::          Size of the hhout, hgout, ghout, ggout       
 *---------------------*/
{
register int i,j;

for(i=0; i< sl; ++i)
  for(j=0; j< sl; ++j)  {
    ACCESS3D(am, D1, D12, J, x+i, y+j) = ACCESS(hhout, nm, ix+i, jy+j);
    ACCESS3D(am, D1, D12, J, x+i, sl+y+j) = ACCESS(hgout, nm, ix+i, jy+j);
    ACCESS3D(am, D1, D12, J, sl+x+i, y+j) = ACCESS(ghout, nm, ix+i, jy+j);
    ACCESS3D(am, D1, D12, J, sl+x+i, sl+y+j) = ACCESS(ggout, nm, ix+i, jy+j);
    }
}

/* initialise the answer matrix */

void initSWT2D(double *m, int *nm, double *am, int *J, double *H,
	int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *m::     The input data                   
int *nm::       The dimension of the square matrix m         
double *am::    The *big* answer 3D array                
int *J::        The level at which to store the initial information  
double *H::     The smoothing filter             
int *LengthH::	The length of the smoothing filter       
int *error::    Error code 0=ok, anything else is memory error   
 *---------------------*/
{
int mlength;    /* Length of vector representing matrix         */
int D1, D12;    /* 1st and Second dimension of answer matrix        */
double *hhout, *hgout, *ghout, *ggout;  /* Intermediate stores      */
int nm2;    /* Half of *nm                      */

/* Carries out a step of the SWT2D algorithm        */
void SWT2D(double *m, int *nm,
	double *hhout, double *hgout, double *ghout, double *ggout,
	double *H, int *LengthH, int *error);

*error = 0;

mlength = *nm * *nm;

/* First create some space for hhout, hgout, ghout and ggout.       */

if ((hhout = (double *)malloc((size_t)(mlength* sizeof(double))))==NULL){
    *error = 7;
    return;
    }

if ((hgout = (double *)malloc((size_t)(mlength* sizeof(double))))==NULL){
    *error = 8;
    return;
    }

if ((ghout = (double *)malloc((size_t)(mlength* sizeof(double))))==NULL){
    *error = 9;
    return;
    }

if ((ggout = (double *)malloc((size_t)(mlength* sizeof(double))))==NULL){
    *error = 10;
    return;
    }

/* Apply the 2D SWT to the initial data and store the HH, GH, HG, GG
 * matrices in their appropriate place in the big matrix
 */

SWT2D(m, nm, hhout, hgout, ghout, ggout, H, LengthH, error);

#ifdef PRINTON
Rprintf("First hhout matrix\n");
{
int i,j;
for(i=0; i<*nm; ++i)    {
    Rprintf("[%d, ] ", i);
    for(j=0; j<*nm; ++j)
      Rprintf("%lf ", ACCESS(hhout, *nm, i,j));
    Rprintf("\n");
    }
}
#endif
    

if (*error != 0)
    return;

/*
 * Now copy each of the results hhout, hgout, ghout and ggout to the answer
 * matrix am
 */


D12 = (*J)*(*nm * 2);
D1 = (*J);

nm2 = *nm / 2;

SmallStore(am, D1, D12, *J-1, nm2, 0l, 0l, 0l, 0l,
    hhout, hgout, ghout, ggout, *nm);
SmallStore(am, D1, D12, *J-1, nm2, *nm, 0l, nm2, 0l,
    hhout, hgout, ghout, ggout, *nm);
SmallStore(am, D1, D12, *J-1, nm2, 0l, *nm, 0l, nm2,
    hhout, hgout, ghout, ggout, *nm);
SmallStore(am, D1, D12, *J-1, nm2, *nm, *nm, nm2, nm2,
    hhout, hgout, ghout, ggout, *nm);


free((void *)hhout);
free((void *)hgout);
free((void *)ghout);
free((void *)ggout);
}

    

void SWT2Drec(double *am, int D1, int D12, int x, int y, int TWOsl,
	int sl, int J, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *am::    The big storage array                
int D1::        First dimension of am                
int D12::       First and second dimensions of am multiplied     
int x::         X origin coordinate of smoothed data         
int y::         Y origin coordinate of smoothed data         
int TWOsl::     Side length of smoothed data             
int sl::        Side length of result packets (2*sl = TWOsl)     
int J::         Level we accessing from (and putting into j-1)   
double *H::     The smoothing filter             
int *LengthH::	The length of the smoothing filter       
int *error::    Error code                       
 *---------------------*/
{
register int i,j;
double *m;  /* Somewhere to put the smoothed data           */
int mlength;    /* The length of this matrix                */
double *hhout, *hgout, *ghout, *ggout; /* Smoothed, hori, verti & diag  */
int sl2;    /* sl divided by 2                  */

void SmallStore(double *am, int D1, int D12, int J, int sl, int x, int y,
	int ix, int jy,
	double *hhout, double *hgout, double *ghout, double *ggout, int nm);
void SWT2D(double *m, int *nm,
	double *hhout, double *hgout, double *ghout, double *ggout,
	double *H, int *LengthH, int *error);
void SWT2Drec(double *am, int D1, int D12, int x, int y, int TWOsl,
	int sl, int J, double *H, int *LengthH, int *error);

*error = 0;

#ifdef PRINTON
Rprintf("SWT2Drec: x=%ld, y=%ld, TWOsl=%ld, sl=%ld, J=%ld\n",x,y,TWOsl,sl,J);
#endif

mlength = TWOsl * TWOsl;

/* Create space for TWOsl * TWOsl matrix m*/

if ((m = (double *)malloc((size_t)mlength*sizeof(double)))==NULL){
    *error = 11;
    return;
    }

/* Fill matrix from am from x,y at origin at level j*/

for(i=0; i<TWOsl; ++i)
 for(j=0; j<TWOsl; ++j)
  ACCESS(m, TWOsl, i, j) = ACCESS3D(am, D1, D12, J, x+i, y+j);

/* Create space for hhout, etc */

if ((hhout = (double *)malloc((size_t)mlength*sizeof(double)))==NULL){
    *error = 12;
    return;
    }
if ((hgout = (double *)malloc((size_t)mlength*sizeof(double)))==NULL){
    *error = 13;
    return;
    }

if ((ghout = (double *)malloc((size_t)mlength*sizeof(double)))==NULL){
    *error = 14;
    return;
    }

if ((ggout = (double *)malloc((size_t)mlength*sizeof(double)))==NULL){
    *error = 15;
    return;
    }

/* Process using SWT2D and into hhout, hgout, ghout, hhout */

SWT2D(m, &TWOsl, hhout, hgout, ghout, ggout, H, LengthH, error);

if (*error != 0)
    return;

#ifdef PRINTON
Rprintf("HH matrix for level %ld\n", J); 
for(i=0; i<TWOsl; ++i)  {
 Rprintf("[%d,] ",i);
 for(j=0; j<TWOsl; ++j)
  Rprintf("%lf ", ACCESS(hhout, TWOsl, i, j));
 Rprintf("\n");
 }
#endif

/* Free m */

free((void *)m);

/* Small Store the hhout etc into am at level j-1 */

SmallStore(am, D1, D12, J-1, sl, x, y, 0l, 0l, hhout, hgout, ghout, ggout,
    TWOsl);
SmallStore(am, D1, D12, J-1, sl, x+TWOsl, y, sl, 0l, hhout, hgout, ghout, ggout,
    TWOsl);
SmallStore(am, D1, D12, J-1, sl, x, y+TWOsl, 0l, sl, hhout, hgout, ghout, ggout,
    TWOsl);
SmallStore(am, D1, D12, J-1, sl, x+TWOsl, y+TWOsl, sl, sl, hhout, hgout, ghout,
    ggout, TWOsl);

/* Free hhout, etc */

free((void *)hhout); free((void *)hgout);
free((void *)ghout); free((void *)ggout);

/* Call SWTrec again 4 times if J != 0 , otherwise return */

if (J == 1)
    return;

sl2 = sl/2;

SWT2Drec(am, D1, D12, x, y, sl, sl2, J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, x+TWOsl, y, sl, sl2, J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, x, y+TWOsl, sl, sl2, J-1, H, LengthH, error);
if (*error != 0)
    return;
SWT2Drec(am, D1, D12, x+TWOsl, y+TWOsl, sl, sl2, J-1, H, LengthH, error);
if (*error != 0)
    return;
}


/* Performs "head" part of 2D packet SWT algorithm 
 *
 * This is just the decomposition step. The full algorithm is implemented
 * in S
 */

void SWT2D(double *m, int *nm,
	double *hhout, double *hgout, double *ghout, double *ggout,
	double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *m::       The matrix to decompose          
int *nm::         The dimension of the square matrix m     
double *hhout::   The smoothed-smoothed matrix of dimension m  
double *ghout::   The detail-smoothed matrix of dimension m    
double *hgout::   The smoothed-detail matrix of dimension m    
double *ggout::   The detail-detail matrix of dimension m  
double *H::       The smoothing filter             
int *LengthH::    The length of the smoothing filter       
int *error::      Error code 0=ok, 1=memory error      
 *---------------------*/
{
int mlength;        /* The number of items in h and g       */
double *h;      /* Intermediate smoothed matrix         */
double *g;      /* Intermediate detail matrix           */

/* Apply H and G smoothers across rows      */
void SWT2DROWblock(double *m, int *nm, double *hout, double *gout,
	double *H, int LengthH, int *error);
/* Apply H and G smoothers across cols      */
void SWT2DCOLblock(double *m, int *nm, double *hout, double *gout,
	double *H, int LengthH, int *error);

*error = 0;

mlength = *nm * *nm;


/* First need to create intermediate matrices h and g   */

if ((h = (double *)malloc((size_t)(mlength* sizeof(double))))==NULL){
    *error = 3;
    return;
    }
if ((g = (double *)malloc((size_t)(mlength * sizeof(double))))==NULL){
    *error = 4;
    return;
    }

/* Fill up the h and g matrices by "across rows" operator       */

SWT2DROWblock(m, nm, h, g, H, (int)*LengthH, error);

if (*error != 0)
    return;

/* Now for each of h and g apply the "across cols" operator     */

SWT2DCOLblock(h, nm, hhout, hgout, H, (int)*LengthH, error);

if (*error != 0)
    return;

SWT2DCOLblock(g, nm, ghout, ggout, H, (int)*LengthH, error);

if (*error != 0)
    return;

/* That's it */

free((void *)h);
free((void *)g);
}

/*
 * Performs part of 2D packet SWT algorithm. 
 *
 * Takes matrix m and applies the G and H filters across the rows of
 * the matrix. Then this is done again but to the rotated versions.
 */ 

void SWT2DROWblock(double *m, int *nm, double *hout, double *gout,
	double *H, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *m::       Input matrix                 
int *nm::         Dimension of square input matrix     
double *hout::    Smoothed answer matrix (same dim as m)   
double *gout::    Detail answer matrix (same dim as m)     
double *H::       The smoothing filter             
int LengthH::     The length of the smoothing filter       
int *error::      Error code. 0=ok, 1=memory error     
 *---------------------*/
{
register int i,j;
int lengthout;      /* Length of out matrices (half of *nm)     */
double *col;        /* A container for columns of matrices      */
double *out;        /* A vector for storing the processed coeffs    */

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);
void rotater(double *book, int length);

*error = 0;

/* Get some memory for a column container called col    */

if ((col = (double *)malloc((unsigned)*nm*sizeof(double)))==NULL)   {
    *error = 1;
    return;
    } 

/* Get some memory for results of a filtering op    */ 

lengthout = (int)*nm/2;

if ((out = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL) {
    *error = 2;
    return;
    } 

/* Process one column at a time */

for(j=0; j< *nm; ++j)   {

    /* Copy column j */

    for(i=0; i< *nm; ++i)
        *(col+i) = ACCESS(m, *nm, i, j);

    
    /* Apply DWT-smooth to j */

    convolveC(col, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this column */

    for(i=0; i < lengthout; ++i)
        ACCESS(hout, *nm, i, j) = *(out+i);

    /* Apply DWT-detail to j */ 

    convolveD(col, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this column */

    for(i=0; i < lengthout; ++i)
        ACCESS(gout, *nm, i, j) = *(out+i);


    /* Rotate column j */

    rotater(col, (int)*nm);

    /* Apply DWT-smooth to rotated j */

    convolveC(col, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this column */

    for(i=0; i < lengthout; ++i)
        ACCESS(hout, *nm, i+lengthout, j) = *(out+i);

    /* Apply DWT-detail to j */ 

    convolveD(col, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this column */

    for(i=0; i < lengthout; ++i)
        ACCESS(gout, *nm, i+lengthout, j) = *(out+i);
    }

/* That's it */
free((void *)col);
free((void *)out);
return;
}

/*
 * Performs part of 2D packet SWT algorithm. 
 *
 * Takes matrix m and applies the G and H filters across the cols of
 * the matrix. Then this is done again but to the rotated versions.
 */ 


void SWT2DCOLblock(double *m, int *nm, double *hout, double *gout,
	double *H, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *m::       Input matrix                 
int *nm::         Dimension of square input matrix     
double *hout::    Smoothed answer matrix (same dim as m)   
double *gout::    Detail answer matrix (same dim as m)     
double *H::       The smoothing filter             
int LengthH::     The length of the smoothing filter       
int *error::      Error code. 0=ok, 1=memory error     
 *---------------------*/
{
register int i,j;
int lengthout;      /* Length of out matrices (half of *nm)     */
double *row;        /* A container for columns of matrices      */
double *out;        /* A vector for storing the processed coeffs    */

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);
void rotater(double *book, int length);

*error = 0;

/* Get some memory for a column container called row    */

if ((row = (double *)malloc((unsigned)*nm*sizeof(double)))==NULL)   {
    *error = 5;
    return;
    } 

/* Get some memory for results of a filtering op    */ 

lengthout = (int)*nm/2;

if ((out = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL) {
    *error = 6;
    return;
    } 

/* Process one row at a time */

for(i=0; i< *nm; ++i)   {

    /* Copy row i */

    for(j=0; j< *nm; ++j)
        *(row+j) = ACCESS(m, *nm, i, j);

    
    /* Apply DWT-smooth to i */

    convolveC(row, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this row */

    for(j=0; j < lengthout; ++j)
        ACCESS(hout, *nm, i, j) = *(out+j);

    /* Apply DWT-detail to i */ 

    convolveD(row, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this row */

    for(j=0; j < lengthout; ++j)
        ACCESS(gout, *nm, i, j) = *(out+j);


    /* Rotate row i */

    rotater(row, (int)*nm);

    /* Apply DWT-smooth to rotated i */

    convolveC(row, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this row */

    for(j=0; j < lengthout; ++j)
        ACCESS(hout, *nm, i, j+lengthout) = *(out+j);

    /* Apply DWT-detail to rotated i */ 

    convolveD(row, (int)*nm, 0, H, LengthH, out,
        0, lengthout-1, WAVELET, 1, PERIODIC);  

    /* Store this row */

    for(j=0; j < lengthout; ++j)
        ACCESS(gout, *nm, i, j+lengthout) = *(out+j);
    }

/* That's it */
free((void *)row);
free((void *)out);
return;
}

void ixtoco(int *level, int *maxlevel, int *index, int *x, int *y)
{
register int i;
int lic;
int mf;

mf = 0x01 << *level;

for(i=*level; i<= *maxlevel; ++i)        {
        lic = *index % 10;
        *index /= 10;
        *x += ((lic & 0x01) *  mf) << 1;
        *y += ((lic & 0x02) *  mf) ;
        mf <<= 1;
        }
}


/* A reconstruction at level levj   */

void SWTRecon(double *am, int D1, int D12, int levj, double *out,
	int x, int y, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *am::  The big storage array                
int D1::      First dimension of am                
int D12::     First and second dimensions of am multiplied     
int levj::    The level to reconstruct             
double *out:: The matrix where the reconstruction is put       
int x::       The x coordinate of the origin           
int y::       The y coordinate of the origin           
double *H::   The filter                       
int *LengthH::The length of the filter             
int *error::  Error code 0=ok                  
 *---------------------*/
{
register int i,j;
int sl; /* Side length of matrix holding details        */
int matele; /* Number of elements in matrix             */
int bc=PERIODIC;    /* Periodic boundary conditions         */
double *hhout, *hgout, *ghout, *ggout;

/* Get the smooth coefficients at partclr level */
void SWTGetSmooth(double *am, int D1, int D12, double *TheSmooth, int levj,
	int x, int y, int sl, double *H, int *LengthH, int *error);

void ImageReconstructStep(double *ImCC, double *ImCD, double *ImDC,
	double *ImDD,
	int LengthCin, int firstCin,
	int LengthDin, int firstDin,
	double *H, int LengthH,
	int LengthCout, int firstCout, int lastCout,
	double *ImOut, int *bc, int *error);
void tpose(double *m, int l);

*error = 0;

sl = 1 << levj; /* Dim of detail arrays is 2^levj   */
matele = sl*sl;     /* Number of elements in matrix     */

/* Get memory for horizontal, vertical and diagonal detail and smooth   */

if ((hhout = (double *)malloc((size_t)matele*sizeof(double)))==NULL)    {
    *error = 16;
    return;
    }

if ((hgout = (double *)malloc((size_t)matele*sizeof(double)))==NULL)    {
    *error = 17;
    return;
    }

if ((ghout = (double *)malloc((size_t)matele*sizeof(double)))==NULL)    {
    *error = 18;
    return;
    }

if ((ggout = (double *)malloc((size_t)matele*sizeof(double)))==NULL)    {
    *error = 19;
    return;
    }

/* Get hold of the horizontal, vertical and diagonal detail     */

for(i=0; i<sl; ++i)
  for(j=0; j<sl; ++j)   {
    ACCESS(hgout, sl, i, j) = ACCESS3D(am, D1, D12, levj, x+i, y+sl+j);
    ACCESS(ghout, sl, i, j) = ACCESS3D(am, D1, D12, levj, x+sl+i, y+j);
    ACCESS(ggout, sl, i, j) = ACCESS3D(am, D1, D12, levj, x+sl+i, y+sl+j);
    }

/*
 * Get hold of the smooth at this level 
 *
 * If levj==0 then this is just the final smoothed scaling function coefficients
 *
 * If levj>0 then we have to recursively get hold of the smooth at the j-1
 * level
 * 
 */

if (levj == 0)
    ACCESS(hhout, sl, 0, 0) = ACCESS3D(am, D1, D12, levj, x, y);
else    {
    SWTGetSmooth(am, D1, D12, hhout, levj, x, y, sl, H, LengthH, error);
    if (*error != 0)
        return;
    }

/*
 * Use S,H,V, and D to reconstruct at level levj, x, y 
 * and put it into out.
 */

#ifdef PRINTON
Rprintf("This is ggout\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(ggout, sl, i, j));
  Rprintf("\n");
  }
Rprintf("This is hgout\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(hgout, sl, i, j));
  Rprintf("\n");
  }
Rprintf("This is ghout\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(ghout, sl, i, j));
  Rprintf("\n");
  }
#endif

/*
tpose(hhout, sl); - this doesn't need to be done because its already been done
*/
tpose(hgout, sl);
tpose(ghout, sl);
tpose(ggout, sl);

ImageReconstructStep(hhout, hgout, ghout, ggout, sl, 0, sl, 0,
    H, *LengthH, 2*sl, 0, 2*sl-1, out, &bc, error);

#ifdef PRINTON
Rprintf("Immediate IRS\n");
for(i=0; i<2*sl; ++i)   {
    Rprintf("[%d, ] ", i);
    for(j=0; j<2*sl; ++j)   {
        Rprintf("%lf ", ACCESS(out, 2*sl, i,j));
        }
    Rprintf("\n");
    }
#endif

if (*error != 0)
    return;

/*
 * Free memory if necessary
 */

free((void *)hhout);
free((void *)hgout);
free((void *)ghout);
free((void *)ggout);
}

void SAvBasis(double *am, int *D1, int *D12, double *TheSmooth, int *levj,
	double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *am::          The big storage array                
int *D1::             First dimension of am                
int *D12::            First and second dimensions of am multiplied     
double *TheSmooth::   The returned smooth              
int *levj::           The level you want it at         
double *H::           The filter                   
int *LengthH::        Length of the filter             
int *error::          Error code. O=ok             
 *---------------------*/
{
void tpose(double *m, int l);
/* Get the smooth coefficients at partclr level */
void SWTGetSmooth(double *am, int D1, int D12, double *TheSmooth, int levj,
	int x, int y, int sl, double *H, int *LengthH, int *error);

*error = 0;

SWTGetSmooth(am, *D1, *D12, TheSmooth, *levj, 0, 0, 1l<<*levj,
    H, LengthH, error);

if (*error != 0)
    return;

tpose(TheSmooth, 1<<*levj);

}

void SWTGetSmooth(double *am, int D1, int D12, double *TheSmooth, int levj,
	int x, int y, int sl, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *am::          The big storage array                
int D1::              First dimension of am                
int D12::             First and second dimensions of am multiplied     
double *TheSmooth::   The returned smooth              
int levj::            The level you want it at         
int x::               The x ordinate for the origin        
int y::               The y ordinate for the origin        
int sl::              Sidelength of TheSmooth array        
double *H::           The filter                   
int *LengthH::        Length of the filter             
int *error::          Error code. O=ok             
 *---------------------*/
{
register int i,j;
double *holder, *holder2;   /* Some storage space           */
double *rc0, *rc1, *rc2, *rc3;
            /* The four reconstructions         */
int matele;     /* Number of elements in the recon matrices */

void tpose(double *m, int l);
void rotateback(double *book, int length);

/* Get a reconstruction at level,x,y        */
void SWTRecon(double *am, int D1, int D12, int levj, double *out,
	int x, int y, double *H, int *LengthH, int *error);

/*
 * Get memory for SWTRecon
 */

matele = sl*sl;

if ((rc0 = (double *)malloc((size_t)matele*sizeof(double)))==NULL)  {
    *error = 20;
    return;
    } 

if ((rc1 = (double *)malloc((size_t)matele*sizeof(double)))==NULL)  {
    *error = 21;
    return;
    } 

if ((rc2 = (double *)malloc((size_t)matele*sizeof(double)))==NULL)  {
    *error = 22;
    return;
    } 

if ((rc3 = (double *)malloc((size_t)matele*sizeof(double)))==NULL)  {
    *error = 24;
    return;
    } 


/*
 * Perform SWTRecon for each of the 4 indices
 *
 * and maybe do some rotation
 */
SWTRecon(am, D1, D12, levj-1, rc0, x, y, H, LengthH, error);
if (*error != 0) return;
SWTRecon(am, D1, D12, levj-1, rc1, x+sl, y, H, LengthH, error);
if (*error != 0) return;
SWTRecon(am, D1, D12, levj-1, rc2, x, y+sl, H, LengthH, error);
if (*error != 0) return;
SWTRecon(am, D1, D12, levj-1, rc3, x+sl, y+sl, H, LengthH, error);
if (*error != 0) return;

/* Rotation -
 *
 * rc0  -   none
 * rc1  -   rotate rows
 * rc2  -   rotate columns
 * rc3  -   rotate rows and columns
 *
 */ 

if ((holder = (double *)malloc((size_t)sl*sizeof(double)))==NULL)   {
    *error = 25;
    return;
    }
if ((holder2 = (double *)malloc((size_t)sl*sizeof(double)))==NULL)  {
    *error = 26;
    return;
    }


/* Rotate row by row */

for(i=0; i<sl; ++i) {

    /* Get row i */
    for(j = 0; j<sl; ++j)   {
        *(holder+j) = ACCESS(rc1,sl,i,j);
        *(holder2+j) = ACCESS(rc3,sl,i,j);
        }

    rotateback(holder, (int)sl);
    rotateback(holder2, (int)sl);

    for(j=0; j<sl; ++j) {
        ACCESS(rc1,sl,i,j) = *(holder+j);
        ACCESS(rc3,sl,i,j) = *(holder2+j);
        }
    }

/* Rotate col by col */

for(j=0; j<sl; ++j) {

    /* Get col j */
    for(i = 0; i<sl; ++i)   {
        *(holder+i) = ACCESS(rc2,sl,i,j);
        *(holder2+i) = ACCESS(rc3,sl,i,j);
        }

    rotateback(holder, (int)sl);
    rotateback(holder2, (int)sl);

    for(i=0; i<sl; ++i) {
        ACCESS(rc2,sl,i,j) = *(holder+i);
        ACCESS(rc3,sl,i,j) = *(holder2+i);
        }
    }
#ifdef PRINTON
Rprintf("This is for level %ld\n", levj);
Rprintf("r0\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(rc0, sl, i, j));
  Rprintf("\n");
  }
Rprintf("r1\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ",i );
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(rc1, sl, i, j));
  Rprintf("\n");
  }
Rprintf("r2\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(rc2, sl, i, j));
  Rprintf("\n");
  }
Rprintf("r3\n");
for(i=0; i<sl; ++i) {
  Rprintf("[%d,] ", i);
  for(j=0; j<sl; ++j)
    Rprintf("%lf ", ACCESS(rc3, sl, i, j));
  Rprintf("\n");
  }
#endif

/*
 * Free holder memory
 */

free((void *)holder);
free((void *)holder2);


/*
 * Work out the average and put it into TheSmooth
 */

for(i=0; i<sl; ++i)
  for(j=0; j<sl; ++j)   
    ACCESS(TheSmooth, sl, i,j) = (ACCESS(rc0,sl,i,j)+
                  ACCESS(rc1,sl,i,j)+
                  ACCESS(rc2,sl,i,j)+
                  ACCESS(rc3,sl,i,j))/4.0;

#ifdef PRINTON
Rprintf("Transposed TheSmooth\n");
/*
tpose(TheSmooth, sl);
*/

for(i=0; i<sl; ++i) {
    Rprintf("[%d, ] ", i);
    for(j=0; j<sl; ++j) {
        Rprintf("%lf ", ACCESS(TheSmooth, sl, i,j));
        }
    Rprintf("\n");
    }
#endif


/*
 * Free memory
 */

free((void *)rc0);
free((void *)rc1);
free((void *)rc2);
free((void *)rc3);
}

void tpose(double *m, int l)
{
double tmp;
register int i,j;

for(i=0; i<l; ++i)
    for(j=0; j<i; ++j)  {
        tmp = ACCESS(m, l, i, j);
        ACCESS(m, l, i, j) = ACCESS(m, l, j, i);
        ACCESS(m, l, j, i) = tmp;
        }
}

void getpacketwst2D(double *am, int *D1, int *D12, int *maxlevel,
	int *level, int *index, int *type, double *out, int *sl)
/*---------------------
 * Argument description
 *---------------------
double *am::     The big storage array                
int *D1::        First dimension of am                
int *D12::       First and second dimensions of am multiplied     
int *maxlevel::  The maximum level (in C numbering)            
int *level::     Which level of coefficients you require      
int *index::     The index number of the packet           
int *type::      The type of coefficients S,H,V or D (TYPES etc.) 
double *out::    A square matrix of dimension sl           
int *sl::        Dimension of out square matrix           
 *---------------------*/
{
register int i,j;
int x,y;    /* The coordinates into am              */

/* Convert index to coordinates             */
void ixtoco(int *level, int *maxlevel, int *index, int *x, int *y);
void tpose(double *m, int l);

x=y=0;
ixtoco(level, maxlevel, index, &x, &y);

switch(*type)   {
    case TYPES:
        break;

    case TYPEH:
        y += *sl;
        break;

    case TYPEV:
        x += *sl;
        break;

    case TYPED:
        x += *sl;
        y += *sl;
        break;
    };

for(i=0; i<*sl; ++i)
  for(j=0; j<*sl; ++j)
    ACCESS(out, *sl, i, j) = ACCESS3D(am, *D1, *D12, *level, x+i, y+j);

tpose(out, *sl);

}
void putpacketwst2D(double *am, int *D1, int *D12, int *maxlevel,
	int *level, int *index, int *type, double *in, int *sl)
/*---------------------
 * Argument description
 *---------------------
double *am::     The big storage array                
int *D1::        First dimension of am                
int *D12::       First and second dimensions of am multiplied     
int *maxlevel::  The maximum level (in C numbering)            
int *level::     Which level of coefficients you require      
int *index::     The index number of the packet           
int *type::      The type of coefficients S,H,V or D (TYPES etc.) 
double *in::     A square matrix of dimension sl           
int *sl::        Dimension of out square matrix           
 *---------------------*/
{
register int i,j;
int x,y;    /* The coordinates into am              */
/* Convert index to coordinates             */
void ixtoco(int *level, int *maxlevel, int *index, int *x, int *y);

x=y=0;
ixtoco(level, maxlevel, index, &x, &y);

tpose(in, *sl);

switch(*type)   {
    case TYPES:
        break;

    case TYPEH:
        y += *sl;
        break;

    case TYPEV:
        x += *sl;
        break;

    case TYPED:
        x += *sl;
        y += *sl;
        break;
    };


for(i=0; i<*sl; ++i)
  for(j=0; j<*sl; ++j)
    ACCESS3D(am, *D1, *D12, *level, x+i, y+j) = ACCESS(in, *sl, i, j);

}

/*
 * The following code was first written in C++ by Markus Monnerjahn,
 * Universitat * Kaiserslautern. The code was rewritten in C
 * (and a few errors corrected) * by Piotr Fryzlewicz, Wroclaw University
 * whilst he was visiting University of * Bristol in 1998-9. The code is
 * designed to integrate with S-Plus
 */

/* The following was in Fryzlewicz's ``Filters.h'' */

#define Nmax 8

typedef struct {
int Length;
double H[2 * Nmax];
double G[2 * Nmax];
double HLeft[Nmax][3 * Nmax - 1];
double GLeft[Nmax][3 * Nmax - 1];
double HRight[Nmax][3 * Nmax - 1];
double GRight[Nmax][3 * Nmax - 1];
double PreLeft[Nmax][Nmax];
double PreInvLeft[Nmax][Nmax];
double PreRight[Nmax][Nmax];
double PreInvRight[Nmax][Nmax];
} Filter;

Filter GetFilt(int);

double Sum(double* vect, int length);

/* The following was in Fryzlewicz's ``WavInt.h'' */

#define NORMAL 0
#define INVERSE 1

#define NO_PRECOND 0
#define PRECOND 1

void Precondition(int Scale, int Direction, Filter F, double* Vect);

void TransStep(int Scale, Filter F, double* Vect);

void InvTransStep(int Scale, Filter F, double* Vect);

void Transform(int MinScale, int Direction, int FilterNumber, double* Vect, int 
Size, int Precond,
int* FilterHistory);

/* The following was in Fryzlewicz's ``Coef.h'' */

const double Interior[] = {
0.70710678, 0.70710678, 0.482962913145, 0.836516303738, 
0.224143868042, -0.129409522551, 0.332670552950, 0.806891509311, 
0.459877502118, -0.135011020010, -0.085441273882, 0.035226291882, 
0.045570345896, -0.0178247014417, -0.140317624179, 0.421234534204, 
1.13665824341, 0.703739068656, -0.0419109651251, -0.107148901418, 
0.0276321529578, -0.0298424998687, -0.247951362613, 0.0234789231361, 
0.89658164838, 1.02305296689, 0.281990696854, -0.0553441861166, 
0.0417468644215, 0.0386547959548, -0.0110318675094, 0.00249992209279, 
0.06325056266, -0.0297837512985, -0.102724969862, 0.477904371333, 
1.11389278393, 0.694457972958, -0.0683231215866, -0.166863215412, 
0.00493661237185, 0.0217847003266, 0.014521394762, 0.00567134268574, 
-0.152463871896, -0.198056706807, 0.408183939725, 1.08578270981, 
0.758162601964, 0.0246656594886, -0.070078291222, 0.0960147679355, 
0.043155452582, -0.0178704316511, -0.0014812259146, 0.0037926585342, 
0.00267279339281, -0.000428394300246, -0.0211456865284, 0.00538638875377, 
0.0694904659113, -0.0384935212634, -0.0734625087609, 0.515398670374, 
1.09910663054, 0.68074534719, -0.0866536154058, -0.202648655286, 
0.0107586117505, 0.0448236230437, -0.000766690896228, -0.0047834585115 
};

const double Left[] = {
0.70710678, 0.70710678, 0.70710678, -0.70710678, 
0.6033325119E+00, -0.7965435169E+00, 0.6908955318E+00, 0.5463927140E+00, 
-0.3983129977E+00, -0.2587922483E+00, 0.3751746045E-01, 0.1003722456E-01, 
0.4573276599E+00, 0.1223510431E+00, 0.8500881025E+00, 0.2274281117E+00, 
0.2238203570E+00, -0.8366029212E+00, -0.1292227434E+00, 0.4830129218E+00, 
0.3888997639E+00, 0.5837809695E+00, -0.8820782813E-01, 0.7936188402E+00, 
-0.8478413085E+00, 0.1609551773E+00, 0.3494874367E+00, -0.5884171124E-01, 
-0.6211483178E+00, -0.3493401824E+00, 0.5225273932E+00, 0.2989205572E+00, 
-0.2000080031E+00, -0.3283013439E+00, 0.3378673486E+00, -0.3322637132E+00, 
-0.3997707705E+00, 0.6982497200E+00, 0.1648201297E+00, -0.2878789995E+00, 
-0.9587863831E-02, 0.1015059001E-02, 0.3712255319E-03, -0.3930133175E-04, 
0.3260097101E+00, -0.3451437111E-01, 0.8016481645E+00, -0.8486981030E-01, 
0.4720552620E+00, 0.1337307023E+00, -0.1400420809E+00, 0.4604064538E+00, 
-0.8542510010E-01, -0.8068932216E+00, 0.3521962365E-01, 0.3326712589E+00, 
0.9097539258E+00, -0.7573970762E-01, 0.4041658894E+00, 0.3254391718E+00, 
0.8904031866E-01, -0.6843490805E+00, -0.1198419201E-01, 0.6200442107E+00, 
-0.3042908414E-01, -0.1885851398E+00, -0.2728514077E+00, 0.1665959592E+00, 
0.5090815232E+00, -0.4847843089E+00, 0.6236424434E+00, 0.3564635425E+00, 
0.4628400863E+00, 0.4839896156E+00, 0.2467476417E+00, -0.6057543651E+00, 
-0.1766953329E-01, 0.3451833289E-01, -0.4517364549E-01, 0.8824901639E-01, 
0.1261179286E+00, 0.2082535326E+00, -0.2308557268E+00, -0.4018227932E+00, 
-0.5279923525E-01, -0.6872148762E-01, 0.2192651713E+00, 0.3302135113E+00, 
0.4634807211E+00, 0.5580212962E+00, 0.7001197140E+00, -0.5994974134E+00, 
0.4120325790E+00, -0.6909199198E-01, -0.2622276250E-01, 0.2785356997E-01, 
-0.6704069413E-01, 0.7120999037E-01, -0.2907980427E-01, 0.6548500701E-01, 
0.5992807229E-01, -0.1349524295E+00, 0.6176427778E-02, -0.1390873929E-01, 
-0.4021099904E-01, 0.9055141946E-01, -0.3952587013E-01, 0.8900857304E-01, 
-0.5259906257E-01, 0.3733444476E+00, 0.3289494480E+00, -0.8404653708E+00, 
0.7966378967E+00, 0.3156849362E+00, 0.4901130336E+00, 0.1202976509E+00, 
-0.2943287768E-01, -0.1307020280E-01, -0.7524762313E-01, -0.3341507090E-01, 
0.9302490657E+00, 0.2803062751E-01, 0.3488878121E+00, -0.1899596351E+00, 
0.1098578445E+00, 0.5126446690E+00, 0.2701025958E-01, -0.6992703507E+00, 
0.7897329608E-02, 0.4475860840E+00, 0.7300852033E-02, -0.1049415950E+00, 
-0.3099501188E+00, -0.6736388320E-01, 0.5987214008E+00, 0.3342473230E+00, 
0.6208411495E+00, -0.5086311754E+00, 0.3720164550E+00, 0.1104945888E-01, 
0.1451060379E+00, 0.6367318477E+00, 0.6263262144E-02, -0.4656214411E+00, 
0.1687072021E-01, 0.3800200446E-01, 0.1562115518E-01, 0.3518730685E-01, 
0.1138401479E+00, 0.1327317986E+00, -0.3641831811E+00, -0.4544562923E+00, 
0.1228916167E-01, 0.2604140629E+00, 0.4117590210E+00, 0.4050120843E+00, 
0.5954973448E+00, -0.8108043041E-02, 0.5457944280E+00, -0.6031278756E+00, 
0.1756416185E+00, -0.3803908839E-01, -0.1751671276E-01, 0.4149398330E+00, 
0.2376581050E-01, -0.5938666471E-01, 0.2200554624E-01, -0.5498806769E-01, 
-0.4098816981E-01, -0.2828796750E-01, 0.1431055982E+00, 0.1003358539E+00, 
-0.6777711666E-01, -0.8695683128E-01, -0.1740980466E+00, -0.7322615732E-01, 
-0.9559132536E-01, 0.9631505948E-01, 0.1512981220E+00, 0.2847632852E+00, 
0.6357421629E+00, -0.7284839004E+00, 0.6823643465E+00, 0.5708436334E+00, 
0.1963091141E+00, -0.5757775517E-02, -0.3301996982E-01, -0.1685142803E+00, 
0.2833562039E-01, 0.2152759089E-01, 0.2623688365E-01, 0.1993310503E-01, 
0.7965029532E-02, 0.1109055091E-01, -0.2796122002E-01, -0.3893335648E-01, 
0.1706238600E-01, 0.2375776006E-01, 0.2928655475E-01, 0.4077875984E-01, 
0.2857806991E-02, 0.3979226173E-02, -0.4686793277E-01, -0.6525916723E-01, 
-0.1702677471E+00, -0.3158187796E-01, 0.3028147873E-01, -0.1797740855E+00, 
0.6351733467E+00, 0.7214516731E+00, 0.7207789094E+00, -0.6372838903E+00, 
0.1995469939E+00, 0.1632013715E-01, -0.3864256419E-01, 0.1755600827E+00, 
0.2947288315E-01, -0.2116688221E-01, 0.2728991267E-01, -0.1959911299E-01, 
0.9231184460E+00, -0.2321705918E-01, 0.3781064509E+00, 0.1171366871E+00, 
0.6815972396E-01, -0.3314072447E+00, -0.1063392935E-01, 0.6083387113E+00, 
-0.9724819798E-02, -0.6313905276E+00, 0.1297887318E-02, 0.3213644020E+00, 
0.5723757426E-02, -0.6248310180E-01, -0.2868924125E+00, 0.6092307458E-01, 
0.5940469387E+00, -0.2496185025E+00, 0.6609604336E+00, 0.4764189916E+00, 
0.3412988002E+00, -0.3565792957E+00, 0.9716877076E-01, -0.2893047969E+00, 
-0.1860855577E-01, 0.6336061692E+00, -0.4013700006E-01, -0.3069767199E+00, 
0.1537307481E-02, 0.4772815467E-02, 0.6783960390E-02, 0.2106188350E-01, 
0.1472870932E+00, -0.1146309215E+00, -0.3468945123E+00, 0.3715476805E+00, 
0.6221043773E-01, -0.4207564723E+00, 0.4938022108E+00, -0.1034481165E+00, 
0.5849105895E+00, 0.4067108104E+00, 0.4598788833E+00, 0.2977351919E+00, 
0.2274224529E+00, -0.6148599157E+00, -0.2788953680E-01, 0.6853115007E-01, 
-0.6367943024E-01, 0.1538100267E+00, 0.2177711568E-02, -0.5448874481E-02, 
0.9609989671E-02, -0.2404525386E-01, -0.8324138021E-01, 0.1150865862E+00, 
0.2035630921E+00, -0.2924528794E+00, -0.8289131955E-01, 0.1456717540E+00, 
-0.1939904105E+00, 0.2552553980E+00, -0.1908457968E-01, -0.1279750816E-02, 
0.2636403540E+00, -0.3051944288E+00, 0.5123094150E+00, -0.4352292068E+00, 
0.6538041187E+00, 0.7155631228E+00, 0.3695601186E+00, -0.1145542755E+00, 
-0.3993953251E-01, -0.4300360725E-01, -0.9479455728E-01, -0.6987749084E-01, 
0.2986558645E-02, 0.4395928052E-02, 0.1317933842E-01, 0.1939872286E-01, 
0.3164853364E-01, 0.7165752442E-01, -0.8214471961E-01, -0.1869771835E+00, 
0.4327551059E-01, 0.1021125523E+00, 0.6830653932E-01, 0.1533850430E+00, 
-0.4840048148E-02, -0.1330465569E-01, -0.7320834090E-01, -0.1616269775E+00, 
-0.5786918156E-01, -0.1222261416E+00, 0.2055855208E-01, -0.2142989488E+00, 
0.4001944206E+00, 0.8103733356E+00, 0.7644983569E+00, -0.4058841524E+00, 
0.4653802512E+00, -0.1135045856E+00, -0.4686421163E-01, 0.2621633969E-01, 
-0.1137136979E+00, 0.6071714788E-01, 0.3413297317E-02, -0.2015600401E-02, 
0.1506248690E-01, -0.8894611814E-02, -0.5550804971E-02, 0.1092693197E-01, 
0.1471087926E-01, -0.2895882270E-01, -0.8859117912E-02, 0.1743944876E-01, 
-0.1158923365E-01, 0.2281376635E-01, 0.1569832259E-02, -0.3090263591E-02, 
0.1157111171E-01, -0.2277809280E-01, 0.7445526684E-02, -0.1465675051E-01, 
0.2902493903E-01, -0.8673422325E-01, -0.2999940975E-01, 0.6576188609E-01, 
-0.6576340468E-01, 0.4762105224E+00, 0.3438258786E+00, -0.7968622388E+00, 
0.7867130346E+00, 0.3399096026E+00, 0.4895844920E+00, 0.7597309556E-01, 
-0.4825468670E-01, -0.2118117647E-01, -0.1177993251E+00, -0.4513809001E-01, 
0.3488475665E-02, 0.1772121221E-02, 0.1539424027E-01, 0.7820166308E-02, 
0.9426568713E+00, 0.2111640418E-01, 0.2691475418E+00, -0.3198838656E+00, 
0.1685897688E+00, 0.5837712653E+00, 0.9379442214E-01, -0.5101698301E+00, 
0.4001600171E-01, 0.4662343677E+00, 0.1057520816E-01, -0.2666472827E+00, 
-0.1804517915E-02, 0.8683768041E-01, 0.4620417362E-02, -0.1225652919E-01, 
-0.3237338194E+00, -0.3337831390E-01, 0.6746802222E+00, 0.3515528841E+00, 
0.5985769136E+00, -0.4361153142E+00, 0.2331427899E+00, 0.2303567833E-01, 
0.7615188825E-02, 0.3942229572E+00, 0.8787238208E-01, -0.5996924995E+00, 
0.1327822780E+00, 0.3966958788E+00, -0.4208296665E-01, -0.1101295389E+00, 
-0.4763766998E-02, -0.2001792798E-02, 0.1219755973E-01, 0.5125562862E-02, 
-0.2713311808E-01, -0.6003418428E-01, -0.2899110477E+00, 0.4010347311E+00, 
0.6661404639E-01, -0.2255897141E+00, 0.7251175151E+00, -0.3612133157E+00, 
0.5750729858E+00, 0.4028777516E+00, 0.4362186445E-01, 0.2083956160E+00, 
-0.1857980605E+00, -0.6015619432E+00, 0.1201045634E+00, 0.2500040680E+00, 
0.5011061230E-01, 0.1224078560E+00, -0.2991531277E-01, -0.9490678016E-01, 
-0.1573393605E-02, -0.3494310966E-02, 0.4028652635E-02, 0.8947135057E-02, 
0.2428297386E-01, -0.9808774833E-02, -0.1012534619E-01, 0.2839107988E-01, 
-0.1589124403E+00, 0.2026374721E-02, -0.1400872185E+00, 0.3959072389E-02, 
0.3332116140E+00, -0.1570876545E-01, 0.7053065273E+00, -0.9207675872E-01, 
0.5793528972E+00, -0.5918868498E-02, 0.3896614469E-02, 0.5375286822E+00, 
-0.5955471765E-01, -0.7665446854E+00, 0.7479188863E-01, 0.2868941871E+00, 
0.3081374376E-01, 0.1404254891E+00, -0.1340206604E-01, -0.1078898335E+00, 
-0.1047350551E-02, -0.4024427396E-02, 0.2681726647E-02, 0.1030449087E-01, 
0.2396415492E-02, 0.1254105463E-02, -0.2199441967E-01, -0.1041703156E-01, 
0.3831273498E-01, 0.1308623534E-01, 0.8416826933E-02, 0.7187262673E-02, 
-0.1404553727E+00, -0.3166691268E-01, -0.9975512873E-01, -0.1089545516E-01, 
0.2671255526E+00, 0.5917635073E-01, 0.7725679637E+00, 0.4986270055E-01, 
0.5400947777E+00, 0.1873760767E-01, 0.1395136806E-01, -0.5366643035E+00, 
-0.4947167749E-01, 0.7678225002E+00, 0.6809969285E-01, -0.2885552645E+00, 
0.3044230228E-01, -0.1401119926E+00, -0.1259630803E-01, 0.1078385285E+00, 
-0.1045026849E-02, 0.4012411756E-02, 0.2675776840E-02, -0.1027372499E-01, 
-0.9493922980E-03, -0.4169893239E-04, 0.8445740060E-02, 0.4623987716E-03, 
-0.8462283108E-02, -0.6221335658E-03, -0.3100522289E-02, 0.1572706697E-03, 
0.2269052507E-01, 0.2991444009E-02, -0.9675176425E-02, 0.4406571767E-03, 
-0.1013545517E+00, -0.1220793447E-01, -0.1405275956E+00, -0.3060471910E-01, 
0.2877816730E+00, 0.6783517122E-01, 0.7686169120E+00, 0.4959733828E-01, 
0.5358841401E+00, 0.1742614880E-01, 0.1743190445E-01, -0.5361004996E+00, 
-0.4949055461E-01, 0.7677690648E+00, 0.6785067390E-01, -0.2886325963E+00, 
0.3050565786E-01, -0.1400474616E+00, -0.1262873733E-01, 0.1078086397E+00, 
-0.1047102354E-02, 0.4010247499E-02, 0.2681091142E-02, -0.1026818344E-01, 
0.1364120169E-03, 0.3562754358E-04, -0.1153014522E-02, -0.3011395265E-03, 
0.1050196147E-02, 0.2742856482E-03, 0.6396328807E-03, 0.1670567176E-03, 
-0.1940525665E-02, -0.5068179562E-03, 0.1962825019E-02, 0.5126423314E-03, 
0.9431646581E-02, 0.2463316874E-02, 0.4027633417E-02, 0.1051919469E-02, 
-0.1076854888E+00, -0.1260408191E-01, -0.1401281648E+00, -0.3053643220E-01, 
0.2886799721E+00, 0.6790539526E-01, 0.7677614604E+00, 0.4955094521E-01, 
0.5360881498E+00, 0.1743687297E-01, 0.1744975917E-01, -0.5360995353E+00, 
-0.4954996630E-01, 0.7677648691E+00, 0.6789074612E-01, -0.2886301380E+00, 
0.3051526119E-01, -0.1400473018E+00, -0.1263602514E-01, 0.1078082911E+00, 
-0.1047379021E-02, 0.4010246532E-02, 0.2681799545E-02, -0.1026818096E-01, 
0.9281136260E+00, -0.7274859222E-02, 0.3668618286E+00, 0.4042756257E-01, 
0.6266997447E-01, -0.1352706631E+00, -0.7647264163E-02, 0.3196302046E+00, 
-0.5088004810E-02, -0.5501818613E+00, 0.1171898539E-02, 0.6205221664E+00, 
0.1656055473E-02, -0.4113728030E+00, -0.1929333000E-03, 0.1430304841E+00, 
-0.1203653707E-02, -0.2044325554E-01, -0.2926634432E+00, 0.2188690374E-01, 
0.6373078980E+00, -0.1065695633E+00, 0.6452435440E+00, 0.2843724182E+00, 
0.2965161471E+00, -0.4589840425E+00, 0.5906254343E-01, 0.3646732234E+00, 
-0.1366965641E-01, 0.1447078292E+00, -0.1323318522E-01, -0.5672211027E+00, 
0.2213708145E-02, 0.4500458402E+00, 0.8306912600E-02, -0.1356216791E+00, 
-0.1892199826E-03, 0.8474645355E-03, -0.1180561737E-02, 0.5287413026E-02, 
0.1537394280E+00, 0.4622773093E-01, -0.3933877381E+00, -0.1945832994E+00, 
0.1368645026E+00, 0.3992968143E+00, 0.6064470308E+00, -0.3762144385E+00, 
0.5685338138E+00, -0.7211460588E-01, 0.3166230149E+00, 0.4092764167E+00, 
0.9637901837E-01, 0.5023569070E-01, -0.2215743488E-01, -0.5779911399E+00, 
-0.4105959919E-01, 0.3834915553E+00, 0.3062891878E-02, -0.1151747954E-01, 
0.1191510029E-01, -0.4999008246E-01, -0.2473506634E-03, 0.7518452020E-03, 
-0.1543244666E-02, 0.4690834793E-02, -0.9143117079E-01, -0.7950461567E-01, 
0.2571933236E+00, 0.2823466127E+00, -0.1858137201E+00, -0.4114017970E+00, 
-0.2318000501E+00, 0.1023384550E+00, 0.1598859051E+00, 0.3519831583E+00, 
0.4871009258E+00, -0.6162537699E-01, 0.5734443542E+00, -0.4034315192E+00, 
0.4450618899E+00, -0.1544031905E+00, 0.2033090498E+00, 0.6034144392E+00, 
-0.3331081673E-01, -0.1069950598E+00, -0.6859620275E-01, -0.2019995950E+00, 
0.4406096842E-02, 0.1468973030E-01, 0.1763488599E-01, 0.5758759998E-01, 
-0.3388224811E-03, -0.1171092606E-02, -0.2113946167E-02, -0.7306559816E-02, 
0.5709470789E-01, 0.6367040385E-01, -0.1665850443E+00, -0.1917735089E+00, 
0.1454047216E+00, 0.1859227361E+00, 0.1015538648E+00, 0.8760929009E-01, 
-0.1362852159E+00, -0.1698024241E+00, -0.1769345243E+00, -0.1650717698E+00, 
0.8596479652E-02, 0.4044496899E-01, 0.3062987789E+00, 0.2651503207E+00, 
0.5431752460E+00, 0.3382384364E+00, 0.6140955155E+00, -0.7711604374E+00, 
0.3284966144E+00, 0.2754111726E+00, -0.4723534035E-01, 0.4710140765E-01, 
-0.1039762360E+00, 0.2958356475E-01, 0.6053819361E-02, -0.8186932177E-02, 
0.2473853830E-01, -0.2780308340E-01, -0.4480360998E-03, 0.8002315701E-03, 
-0.2795340476E-02, 0.4992722014E-02, -0.2771885178E-01, 0.5905246892E-01, 
0.8389262409E-01, -0.1800708576E+00, -0.8147463473E-01, 0.1795432638E+00, 
-0.3866640174E-01, 0.7542773595E-01, 0.7464550842E-01, -0.1632828547E+00, 
0.6920099668E-01, -0.1404076028E+00, -0.2231653624E-01, 0.5304071283E-01, 
-0.1020129060E+00, 0.2049015597E+00, -0.6275852804E-01, 0.1129722066E+00, 
0.8683887629E-01, 0.6297984646E-01, 0.4540741128E+00, -0.7409643915E+00, 
0.7328223558E+00, 0.5042713730E+00, 0.4319437621E+00, 0.7588728111E-01, 
-0.5734179466E-01, -0.4120932567E-01, -0.1314896433E+00, -0.8313424682E-01, 
0.7195965493E-02, 0.5540498594E-02, 0.2979946085E-01, 0.2315912504E-01, 
-0.5190323961E-03, -0.3922287088E-03, -0.3238293222E-02, -0.2447152777E-02, 
0.8529185218E-02, 0.1789547417E-01, -0.2624362023E-01, -0.5515123460E-01, 
0.2697489980E-01, 0.5699017157E-01, 0.9748173003E-02, 0.1992423503E-01, 
-0.2431726904E-01, -0.5131208900E-01, -0.1901319060E-01, -0.3947317762E-01, 
0.8695522469E-02, 0.1852759787E-01, 0.2729552258E-01, 0.5658785929E-01, 
0.1231155208E-01, 0.2487421369E-01, 0.9187597846E-02, 0.5732429572E-01, 
-0.4745631970E-01, -0.1020173091E+00, -0.2924026466E-01, -0.4251128534E+00, 
0.3831732042E+00, 0.8001954561E+00, 0.7723764353E+00, -0.3763869007E+00, 
0.4743837775E+00, -0.6582574501E-01, -0.6082976218E-01, 0.2842730285E-01, 
-0.1418145038E+00, 0.5331847291E-01, 0.7565491348E-02, -0.3918460367E-02, 
0.3149214931E-01, -0.1548884491E-01, -0.5401012058E-03, 0.3080051475E-03, 
-0.3369743560E-02, 0.1921673848E-02, -0.1148179452E-02, 0.2054346689E-02, 
0.3567347972E-02, -0.6382774933E-02, -0.3784785333E-02, 0.6771819551E-02, 
-0.1105680377E-02, 0.1978304192E-02, 0.3387180654E-02, -0.6060419089E-02, 
0.2395725348E-02, -0.4286482122E-02, -0.1281174092E-02, 0.2292309675E-02, 
-0.3407989030E-02, 0.6097652254E-02, -0.1283416571E-02, 0.2296316381E-02, 
-0.9841020551E-02, 0.2254192382E-01, 0.6341196709E-02, -0.1213654065E-01, 
0.4619539219E-01, -0.1379967819E+00, -0.2943416082E-01, 0.6522286075E-01, 
-0.5094613096E-01, 0.4794644197E+00, 0.3655480429E+00, -0.7789628537E+00, 
0.7770499808E+00, 0.3646968778E+00, 0.4810695095E+00, 0.5250853492E-01, 
-0.6126039684E-01, -0.2724349082E-01, -0.1432436897E+00, -0.4923361338E-01, 
0.7606483137E-02, 0.3810616823E-02, 0.3168981524E-01, 0.1496216731E-01, 
-0.5420982042E-03, -0.3029798461E-03, -0.3382203026E-02, -0.1890320507E-02
};

const double LeftPre[] = {
0.3248940489E+00, 0.3077926491E+01, 0.3715801512E-01, -0.1142045672E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1001445405E+01, 0.9985566812E+00, 
0.1007941579E+00, 0.9921209927E+01, -0.5929310244E+00, 0.2752403722E+02, 
-0.1509649764E-01, -0.6946796955E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.2137256653E+00, 0.4678895250E+01, 0.3068542369E-01, -0.1435467050E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1000189333E+01, 0.9998107029E+00, 0.2489911114E+01, 0.4016207624E+00, 
-0.2752988536E+01, 0.6592239446E+00, 0.1687841447E+01, -0.1870967456E+00, 
-0.4022221173E+00, 0.3352374610E-01, 0.0000000000E+00, 0.0000000000E+00, 
0.1677210550E+01, 0.5962280646E+00, -0.7075375437E+00, 0.3732739492E+00, 
0.1763544288E+00, -0.8158414568E-01, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1130145142E+01, 0.8848420994E+00, 
-0.6162121550E-01, 0.5415219932E-01, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1006885156E+01, 0.9931619247E+00, 0.3743818971E+01, 0.2671069322E+00, 
-0.6861305063E+01, 0.1047409836E+01, 0.7183286173E+01, -0.6066998614E+00, 
-0.3711751230E+01, 0.2622361574E+00, 0.7595392314E+00, -0.5336282441E-01, 
0.0000000000E+00, 0.0000000000E+00, 0.1749746932E+01, 0.5715112179E+00, 
-0.1112386159E+01, 0.5118280722E+00, 0.5611031648E+00, -0.2021081048E+00, 
-0.1136294342E+00, 0.3953339052E-01, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1242099062E+01, 0.8050887651E+00, 
-0.2151591879E+00, 0.1662752175E+00, 0.4303743724E-01, -0.3186979674E-01, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1041780293E+01, 0.9598952930E+00, 
-0.1641111826E-01, 0.1572808312E-01, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1001581378E+01, 0.9984211186E+00, 
0.2691258238E+01, 0.3715734097E+00, -0.4691390835E+01, 0.7677024985E+00, 
0.5650878809E+01, -0.1862336274E+00, -0.4003610222E+01, 0.4546833284E-01, 
0.1548107666E+01, -0.9293368857E-02, -0.2523403787E+00, 0.1131936169E-02, 
0.0000000000E+00, 0.0000000000E+00, 0.2270666166E+01, 0.4403993925E+00, 
-0.2346218035E+01, 0.6446164864E+00, 0.1686964492E+01, -0.2228021018E+00, 
-0.6557730032E+00, 0.6332382245E-01, 0.1069443363E+00, -0.8983535001E-02, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1602926731E+01, 0.6238588332E+00, -0.7485464376E+00, 0.3995416825E+00, 
0.3030447729E+00, -0.1345491657E+00, -0.5083020320E-01, 0.2195487513E-01, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1168807480E+01, 0.8555728955E+00, 
-0.1287884925E+00, 0.1077442661E+00, 0.2180567976E-01, -0.1782040660E-01, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1022680347E+01, 0.9778226434E+00, -0.7652542041E-02, 0.7478034810E-02, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1000641088E+01, 0.9993593230E+00, 
0.5804269091E+00, 0.1722869812E+01, -0.7323447716E+00, 0.5738683253E+01, 
0.2824499599E+01, -0.1245382394E+02, -0.2917456460E+01, 0.7964976943E+01, 
0.1789930270E+01, -0.4649550414E+01, -0.6097284388E+00, 0.1520565899E+01, 
0.8862324360E-01, -0.2151255312E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.2198648442E+00, 0.4548248738E+01, 0.5966667693E+00, -0.4076684813E+01, 
-0.9937306950E-01, 0.1227270380E+01, 0.5522044783E-01, -0.6678395557E+00, 
-0.1473788829E-01, 0.1967618134E+00, 0.1693029724E-02, -0.2546662797E-01, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.6656852331E+00, 0.1502211481E+01, 0.1901881588E+00, -0.2856934692E+00, 
-0.1038800577E+00, 0.1544684019E+00, 0.3255433820E-01, -0.4820986746E-01, 
-0.4462431097E-02, 0.6604466353E-02, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1000032785E+01, 0.9999672157E+00, -0.4315341245E-02, 0.4305484704E-02, 
0.1542073127E-02, -0.1534917109E-02, -0.2143390832E-03, 0.2132631237E-03, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1002256439E+01, 0.9977486406E+00, -0.1554134553E-02, 0.1550217349E-02, 
0.2204062611E-03, -0.2197888501E-03, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1000269825E+01, 0.9997302475E+00, -0.7738658261E-04, 0.7736527399E-04, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1000005602E+01, 0.9999943980E+00, 0.2810123818E+01, 0.3558562059E+00, 
-0.6559223197E+01, 0.8146829794E+00, 0.1166265739E+02, -0.1366545660E+00, 
-0.1370803617E+02, -0.2512978838E-01, 0.1057161896E+02, 0.3435512306E-01, 
-0.5161014832E+01, -0.1989050428E-01, 0.1450032681E+01, 0.6413814665E-02, 
-0.1790147104E+00, -0.8795283072E-03, 0.0000000000E+00, 0.0000000000E+00, 
0.2865090274E+01, 0.3490291420E+00, -0.4726642907E+01, 0.7527005851E+00, 
0.5658991531E+01, -0.1741151174E+00, -0.4391120800E+01, 0.2189754588E-01, 
0.2150489727E+01, 0.7628539402E-02, -0.6053327281E+00, -0.4888651292E-02, 
0.7483245389E-01, 0.8230905056E-03, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.2191756126E+01, 0.4562551409E+00, 
-0.2253964453E+01, 0.6427198947E+00, 0.1808275069E+01, -0.2393601344E+00, 
-0.9018774144E+00, 0.8134531543E-01, 0.2565106244E+00, -0.1922603801E-01, 
-0.3190347407E-01, 0.2116725807E-02, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1600048291E+01, 0.6249811369E+00, -0.8330291569E+00, 0.4302638030E+00, 
0.4321721196E+00, -0.1746850867E+00, -0.1262794856E+00, 0.4820621391E-01, 
0.1601926411E-01, -0.6049455037E-02, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1210019309E+01, 0.8264330929E+00, 
-0.2036886238E+00, 0.1611623644E+00, 0.6037258981E-01, -0.4551697161E-01, 
-0.7732501809E-02, 0.5803277507E-02, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1044505769E+01, 0.9573905957E+00, -0.2609642521E-01, 0.2489086650E-01, 
0.3366051566E-02, -0.3198214876E-02, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1003760640E+01, 0.9962534497E+00, 
-0.9726417449E-03, 0.9689366947E-03, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1000062955E+01, 0.9999370493E+00
};

const double Right[] = {
0.70710678, 0.70710678, 0.70710678, -0.70710678, 
0.8705087534E+00, -0.2575129195E+00, 0.4348969980E+00, 0.8014229620E+00, 
0.2303890438E+00, -0.5398225007E+00, -0.1942334074E+00, 0.3717189665E+00, 
0.1901514184E+00, -0.3639069596E+00, 0.3749553316E+00, -0.7175799994E+00, 
0.7675566693E+00, 0.4010695194E+00, 0.4431490496E+00, 0.2315575950E+00, 
0.9096849943E+00, 0.7221948764E-01, 0.3823606559E+00, -0.4265622188E+00, 
0.1509872153E+00, 0.8042331410E+00, 0.5896101069E-01, -0.4074776976E+00, 
-0.2904078511E+00, -0.1535052307E+00, 0.4189992290E+00, 0.5223942379E+00, 
0.4969643721E+00, -0.9819800088E-01, 0.4907578307E+00, -0.7678795743E+00, 
0.4643627674E+00, 0.2985152397E+00, 0.1914505442E+00, 0.1230738317E+00, 
0.8183541840E-01, 0.2294775484E+00, -0.1587582156E+00, -0.4451794444E+00, 
-0.9124735623E-01, -0.2558698912E+00, 0.6042558204E-03, 0.1694414802E-02, 
0.7702933610E-01, 0.7598761510E+00, 0.5200601778E+00, 0.1391503267E+00, 
0.7642591993E+00, -0.2725472352E+00, 0.3150938230E+00, -0.1123675716E+00, 
0.9154705188E+00, -0.1982779906E+00, 0.3919142810E+00, 0.6040677868E+00, 
0.5947771124E-01, -0.6495297603E+00, -0.2519180851E-01, 0.4050309541E+00, 
0.6437934569E-01, -0.9924194741E-01, -0.2191626469E+00, -0.2726273506E+00, 
0.4488001781E+00, 0.5092867484E+00, 0.7540005084E+00, 0.6811856626E-01, 
0.3937758157E+00, -0.6735345752E+00, -0.1581338944E+00, 0.4499340827E+00, 
-0.1614201190E-01, 0.2719065540E-01, 0.4126840881E-01, -0.6951519362E-01, 
0.1290078289E-01, 0.4581959340E-02, -0.1390716006E+00, 0.3062032365E-01, 
0.2921367950E-01, 0.1388737145E-01, 0.4606168537E+00, -0.9504834999E-01, 
0.8164119742E+00, -0.3015815012E+00, 0.2986473346E+00, 0.8033636317E+00, 
-0.1027663536E+00, -0.4968427065E+00, -0.1257488211E-01, -0.2963204370E-01, 
0.3214874197E-01, 0.7575680778E-01, -0.6775603652E-02, 0.2880305124E-02, 
0.1913244129E-01, -0.8133189531E-02, -0.1770918425E-01, 0.7528163799E-02, 
-0.6765916174E-01, 0.2876186983E-01, -0.3023588481E-01, 0.1285325684E-01, 
0.4977920821E+00, -0.9920191439E-01, 0.8039495996E+00, -0.2977899841E+00, 
0.2977111011E+00, 0.8037936841E+00, -0.9910804055E-01, -0.4976470524E+00, 
-0.1259895190E-01, -0.2963766018E-01, 0.3221027840E-01, 0.7577116678E-01, 
0.6629994791E+00, 0.2704377396E+00, 0.5590853114E+00, -0.2099389697E+00, 
-0.3638895712E+00, 0.6232723378E+00, -0.3190292420E+00, -0.6204577645E+00, 
-0.8577516141E-01, 0.3239177989E+00, 0.7938922949E-01, -0.6655879286E-01, 
-0.1615976602E+00, -0.4460509080E+00, 0.7429922432E+00, 0.2531551137E+00, 
0.5771521705E+00, -0.3811235714E+00, 0.2773360273E+00, -0.2939476383E+00, 
0.1530007177E-01, 0.6480763210E+00, -0.1063942878E+00, -0.2856775435E+00, 
-0.1216769629E-01, -0.4448843004E-01, 0.1126647052E-01, 0.4119330180E-01, 
0.2576681838E+00, -0.4091805529E+00, -0.3995181391E-01, 0.1427046841E+00, 
0.1290658921E-01, 0.2088087907E-01, 0.4163708782E+00, -0.3899540483E+00, 
0.6469374504E+00, -0.2052536834E+00, 0.5608421202E+00, 0.7190940497E+00, 
0.4995700940E-02, -0.3116714966E+00, -0.1569482750E+00, -0.1763286609E-01, 
-0.2009485730E-01, -0.4039743821E-01, 0.1860648984E-01, 0.3740531781E-01, 
-0.8955158953E-01, -0.1414833249E+00, 0.3895616532E-01, 0.6067058589E-01, 
0.4751157158E-02, -0.1094361815E-01, -0.6596188737E-01, -0.1401642518E+00, 
0.1520881337E-01, -0.7844525367E-01, 0.2466644601E+00, 0.5612750679E-01, 
0.7199998123E+00, 0.5913073605E+00, 0.6131240289E+00, -0.7443064557E+00, 
0.1537319456E-01, 0.1990571237E+00, -0.1721947543E+00, 0.4730213481E-01, 
-0.2083858871E-01, 0.3055544136E-01, 0.1929513523E-01, -0.2829228895E-01, 
0.1546583496E-01, 0.1104876084E-01, -0.6559247482E-02, -0.4685912973E-02, 
0.2728093541E-02, 0.1948944437E-02, 0.1830645747E-01, 0.1307809574E-01, 
0.1708296448E-01, 0.1220403484E-01, 0.2156026058E-01, 0.1540260595E-01, 
-0.3692524465E-01, -0.1737014260E+00, 0.2045334663E+00, -0.1298838460E-01, 
0.7233635779E+00, 0.6336000391E+00, 0.6327590443E+00, -0.7243137518E+00, 
0.1649999072E-01, 0.1992776018E+00, -0.1751674471E+00, 0.3929252813E-01, 
-0.2109311508E-01, 0.2952572161E-01, 0.1953080958E-01, -0.2733883754E-01, 
0.9236675275E+00, -0.5565971353E-01, 0.3793343885E+00, 0.2242498085E+00, 
0.4986459778E-01, -0.5135424238E+00, -0.1688947163E-01, 0.6321916480E+00, 
-0.3311806998E-02, -0.4901838615E+00, 0.2804629640E-02, 0.2042322283E+00, 
-0.1237553624E-01, -0.3508415169E-01, -0.2831955072E+00, 0.1057159278E+00, 
0.6209232596E+00, -0.3394364574E+00, 0.6494435367E+00, 0.4787854504E+00, 
0.3279016732E+00, -0.7187645472E-01, 0.4613465423E-01, -0.5042780615E+00, 
-0.2201105042E-01, 0.5796966704E+00, 0.4772996268E-01, -0.2205355783E+00, 
0.1811083865E-02, -0.3313371482E-02, -0.7992103956E-02, 0.1462152573E-01, 
0.1415889950E+00, 0.1743290418E+00, -0.3264117354E+00, -0.4255531446E+00, 
0.5898718517E-01, 0.2638047100E+00, 0.5662575494E+00, 0.3616438204E+00, 
0.6511436245E+00, -0.2644671055E+00, 0.3434709841E+00, -0.4908547063E+00, 
-0.6897565912E-01, 0.4962107578E+00, -0.2139852264E-01, 0.6828427694E-01, 
0.4514540543E-01, -0.1606899475E+00, 0.1807019829E-02, -0.5156673789E-02, 
-0.7974169836E-02, 0.2275580595E-01, -0.4080929146E-01, 0.3244709516E-02, 
0.1251977353E+00, 0.5111253513E-03, -0.9208945297E-01, 0.4381492218E-02, 
-0.1815441604E+00, -0.3035847066E-01, 0.7492893632E-01, -0.3796933443E-01, 
0.4921890166E+00, 0.7322186910E-01, 0.7598885402E+00, 0.3436818708E+00, 
0.3298535864E+00, -0.7867428461E+00, -0.6885203537E-01, 0.4886190401E+00, 
-0.2056224689E-01, 0.4823811158E-01, 0.4357038095E-01, -0.1176377166E+00, 
0.1729458762E-02, -0.3491719268E-02, -0.7631901801E-02, 0.1540855391E-01, 
0.1473338385E-01, -0.5211034254E-02, -0.3879928306E-01, 0.1379054345E-01, 
0.3222685441E-01, -0.1172287487E-01, 0.4181351746E-01, -0.1651839644E-01, 
-0.4470160402E-01, 0.1484603300E-01, -0.1143592801E+00, 0.4328769202E-01, 
-0.3079264597E-01, 0.1491923547E-01, 0.4953303341E+00, -0.7352728683E-01, 
0.7831063037E+00, -0.3353441076E+00, 0.3365261769E+00, 0.7880177709E+00, 
-0.7161554542E-01, -0.4916362852E+00, -0.2096645980E-01, -0.4832809329E-01, 
0.4445995152E-01, 0.1180259590E+00, 0.1762244193E-02, 0.3492043263E-02, 
-0.7776580122E-02, -0.1540998367E-01, -0.2279778462E-02, 0.1154376053E-02, 
0.5934267127E-02, -0.3004842779E-02, -0.4653920620E-02, 0.2356533579E-02, 
-0.4697326198E-02, 0.2378512204E-02, 0.7905498252E-02, -0.4002984474E-02, 
0.1478426761E-01, -0.7486080223E-02, 0.6356609850E-03, -0.3218698441E-03, 
-0.1184137622E+00, 0.4493383491E-01, -0.4719524338E-01, 0.2049260653E-01, 
0.4914056494E+00, -0.7279363945E-01, 0.7873827096E+00, -0.3377604925E+00, 
0.3378777169E+00, 0.7876649023E+00, -0.7256793624E-01, -0.4910841873E+00, 
-0.2105690525E-01, -0.4831332253E-01, 0.4471373463E-01, 0.1179951690E+00, 
0.1767573251E-02, 0.3490782271E-02, -0.7800096641E-02, -0.1540441905E-01, 
0.9091449027E+00, 0.1058599696E-01, 0.4122956702E+00, -0.4261024318E-01, 
0.8463996964E-02, 0.1626466305E+00, -0.4985469571E-01, -0.3983434261E+00, 
-0.2885243693E-01, 0.6359366721E+00, -0.8468475604E-02, -0.5811383289E+00, 
0.9290995299E-03, 0.2621339682E+00, 0.2377975298E-02, -0.4571068446E-01, 
-0.2699308042E+00, -0.3369311938E-01, 0.6287804833E+00, 0.1173158297E+00, 
0.6425124921E+00, -0.3345684650E+00, 0.3295936378E+00, 0.4993594177E+00, 
0.9680040842E-01, -0.2270999192E+00, -0.1505090528E-02, -0.4092454872E+00, 
-0.2576996713E-01, 0.5883989035E+00, -0.1663793375E-01, -0.2415648560E+00, 
0.7886583178E-03, 0.4918478872E-02, 0.2019348750E-02, 0.1259369707E-01, 
0.2157402688E+00, -0.6868038150E-01, -0.3660961628E+00, 0.2059784468E+00, 
0.9979560410E-01, -0.4260199250E+00, 0.5343306371E+00, 0.3030769028E+00, 
0.5651956498E+00, 0.2528028219E+00, 0.3984651316E+00, -0.2995453713E+00, 
0.2006608171E+00, -0.4153788373E+00, 0.3790996947E-01, 0.5783273640E+00, 
-0.5279896883E-01, -0.1200871578E+00, -0.3509616877E-01, -0.6583871196E-01, 
0.1600620645E-02, 0.3864118919E-02, 0.4098367097E-02, 0.9894023233E-02, 
-0.1219999250E+00, 0.1199301213E+00, 0.2497494925E+00, -0.3021493088E+00, 
-0.1489636838E+00, 0.4092549285E+00, -0.2314941689E+00, 0.2642183651E-01, 
0.4306729538E-01, -0.3532455770E+00, 0.3349988297E+00, -0.2065478346E+00, 
0.5170449138E+00, 0.2227171845E+00, 0.5314908191E+00, 0.5240786066E+00, 
0.3984867666E+00, -0.2234158448E+00, 0.1134176909E+00, -0.3876543947E+00, 
-0.8664578083E-01, 0.1720771534E+00, -0.6206801690E-01, 0.6600629896E-01, 
0.2555164772E-02, -0.5990160538E-02, 0.6542464174E-02, -0.1533772348E-01, 
0.6643773700E-01, 0.4081778499E-01, -0.1425076588E+00, -0.7665915987E-01, 
0.1150002072E+00, 0.2089854289E-01, 0.1110261320E+00, 0.1140640809E+00, 
-0.5456405522E-01, 0.1937409327E-01, -0.1403852745E+00, -0.1187989690E+00, 
-0.9002913124E-01, -0.1838791619E+00, 0.9927107272E-01, -0.1133363121E+00, 
0.3252695570E+00, 0.7805202816E+00, 0.5955739150E+00, -0.5430653297E+00, 
0.6177120462E+00, 0.6990288263E-01, 0.2102093627E+00, -0.6927193818E-01, 
-0.1207077310E+00, 0.6202160209E-01, -0.9029360085E-01, 0.1264388999E-01, 
0.3498469159E-02, -0.2337275009E-02, 0.8957782052E-02, -0.5984560439E-02, 
-0.2472698999E-01, 0.1089534551E+00, 0.5451727819E-01, -0.2418936349E+00, 
-0.4913268141E-01, 0.2243741754E+00, -0.3543504434E-01, 0.1488992192E+00, 
0.2981657384E-01, -0.1382661605E+00, 0.5443081267E-01, -0.2312770306E+00, 
0.3080352023E-01, -0.1080114084E+00, -0.2311127518E-01, 0.1488311123E+00, 
-0.2663146891E-01, 0.2352279421E+00, -0.3986350788E-01, 0.4020354092E+00, 
0.1012109980E+00, -0.7029809958E+00, 0.5654368499E+00, -0.6175503239E-02, 
0.7419995864E+00, 0.1485821868E+00, 0.2730680885E+00, 0.1174257523E+00, 
-0.1373359727E+00, -0.5325615968E-01, -0.1050857718E+00, -0.2772249506E-01, 
0.3942767529E-02, 0.1737252935E-02, 0.1009540190E-01, 0.4448212190E-02, 
0.4060965940E-02, 0.1545016950E-01, -0.9082349871E-02, -0.3455430249E-01, 
0.8675066430E-02, 0.3300477019E-01, 0.5263347825E-02, 0.2002469773E-01, 
-0.5426164014E-02, -0.2064414148E-01, -0.8281419894E-02, -0.3150712024E-01, 
-0.2928559727E-02, -0.1114186730E-01, 0.7407888982E-02, 0.2818372382E-01, 
0.1392100926E-02, -0.5377738331E-01, 0.2865053765E-01, 0.1320738411E+00, 
0.5527738563E-01, 0.2388167237E+00, -0.5250033174E-01, -0.7743478585E+00, 
0.2337385282E-01, 0.5554662583E+00, 0.5390029057E+00, -0.7133893666E-02, 
0.7666684946E+00, -0.5462109919E-01, 0.2878132266E+00, -0.7104460368E-01, 
-0.1399692009E+00, 0.3084906050E-01, -0.1077108034E+00, 0.1309957900E-01, 
0.4008607630E-02, -0.1053633670E-02, 0.1026398458E-01, -0.2697814488E-02, 
0.9280961080E+00, -0.1623137602E-01, 0.3685482772E+00, 0.7700181909E-01, 
0.5169019500E-01, -0.2216538064E+00, -0.1072856966E-01, 0.4518914876E+00, 
-0.3547264354E-02, -0.5968036841E+00, 0.1881724221E-02, 0.5361053678E+00, 
0.2253426546E-03, -0.2972611580E+00, -0.4074549423E-03, 0.9152878068E-01, 
0.2542121693E-02, -0.1206728019E-01, -0.2880715337E+00, 0.3847998474E-01, 
0.6422316271E+00, -0.1590369267E+00, 0.6540393284E+00, 0.3562112107E+00, 
0.2730170577E+00, -0.4595861878E+00, 0.4300012686E-01, 0.1637283171E+00, 
-0.1532507378E-01, 0.3534175472E+00, -0.3095082211E-02, -0.5814453781E+00, 
0.3081063250E-02, 0.3686477291E+00, -0.1192117285E-01, -0.9703824753E-01, 
-0.2510321213E-03, -0.5800646631E-03, 0.1566213637E-02, 0.3619079426E-02, 
0.1554699662E+00, -0.7121179707E-01, -0.3928508188E+00, 0.2523023478E+00, 
0.1655006354E+00, -0.4171714520E+00, 0.6186900270E+00, 0.2451267006E+00, 
0.5757733437E+00, 0.2725111911E+00, 0.2765573214E+00, -0.3518314346E+00, 
0.4215137090E-01, -0.2344734683E+00, -0.2022141803E-01, 0.5903810045E+00, 
0.3170527734E-01, -0.3062645760E+00, 0.2969565406E-02, -0.9326190637E-02, 
-0.1133537854E-01, 0.4097674398E-01, -0.2472629557E-03, 0.5916913076E-03, 
0.1542697449E-02, -0.3691619185E-02, -0.9273812077E-01, 0.1179909741E+00, 
0.2525480154E+00, -0.3454318691E+00, -0.1903081753E+00, 0.3669763444E+00, 
-0.2325263183E+00, 0.8571956937E-01, 0.1951123754E+00, -0.3733637204E+00, 
0.5747752998E+00, -0.1532441521E+00, 0.6031889613E+00, 0.3576146478E+00, 
0.3197987067E+00, 0.3379298220E+00, -0.2318110629E-01, -0.4826776512E+00, 
-0.2335276594E-01, -0.1247229279E+00, 0.3904642048E-01, 0.2577717071E+00, 
0.3358324305E-02, 0.1649336963E-01, -0.1296942483E-01, -0.6672983676E-01, 
-0.2744734004E-03, -0.1243663709E-02, 0.1712466040E-02, 0.7759337931E-02, 
0.4096472123E-01, 0.9544232866E-02, -0.1240767851E+00, -0.2688549916E-01, 
0.1278114403E+00, 0.3303656538E-01, 0.7138489316E-01, -0.1099762148E-01, 
-0.1584871307E+00, -0.3331736763E-01, -0.1727866355E+00, 0.3186523735E-01, 
0.1143640280E+00, 0.7533252701E-01, 0.5129116028E+00, -0.5202591621E-01, 
0.7137757767E+00, -0.3899644872E+00, 0.3414317954E+00, 0.7710342652E+00, 
-0.3881329228E-01, -0.4656907083E+00, -0.2526829049E-01, -0.6031146457E-01, 
0.4424255303E-01, 0.1391049995E+00, 0.3575822034E-02, 0.7543265904E-02, 
-0.1394488775E-01, -0.3124469163E-01, -0.2875903358E-03, -0.5438414561E-03, 
0.1794303867E-02, 0.3393079342E-02, -0.1790191282E-01, -0.5044968595E-02, 
0.5295161530E-01, 0.1488052556E-01, -0.5704547538E-01, -0.1581581210E-01, 
-0.1894676577E-01, -0.6949272400E-02, 0.6968256270E-01, 0.2057861785E-01, 
0.4547075753E-01, 0.1765306809E-01, -0.7337937479E-01, -0.1862118927E-01, 
-0.1444748191E+00, -0.4891212508E-01, -0.1409191494E-01, -0.1375364073E-01, 
0.4948704363E+00, 0.5400109528E-01, 0.7618895951E+00, 0.3573536308E+00, 
0.3591100053E+00, -0.7782672053E+00, -0.4810103061E-01, 0.4837479787E+00, 
-0.2674149513E-01, 0.6134832074E-01, 0.4783092078E-01, -0.1436058269E+00, 
0.3754819526E-02, -0.7612327455E-02, -0.1470868404E-01, 0.3171165160E-01, 
-0.2997260590E-03, 0.5426010792E-03, 0.1870019815E-02, -0.3385340511E-02, 
0.4842573875E-02, 0.2223331313E-02, -0.1435818114E-01, -0.6588142232E-02, 
0.1564773617E-01, 0.7160556501E-02, 0.3635704163E-02, 0.1734300817E-02, 
-0.1807616145E-01, -0.8390109930E-02, -0.7951972463E-02, -0.3916508325E-02, 
0.2177252927E-01, 0.9923280681E-02, 0.3159516732E-01, 0.1489081491E-01, 
-0.5393073080E-02, -0.2164190023E-02, -0.1454356948E+00, -0.5005924913E-01, 
-0.5484912621E-01, -0.2423995567E-01, 0.4834662634E+00, 0.5271872910E-01, 
0.7753953496E+00, 0.3632737173E+00, 0.3639514649E+00, -0.7773694792E+00, 
-0.5141318290E-01, 0.4816288882E+00, -0.2717238099E-01, 0.6128840935E-01, 
0.4899043649E-01, -0.1433368916E+00, 0.3804055936E-02, -0.7609148826E-02, 
-0.1492982098E-01, 0.3170285605E-01, -0.3026846009E-03, 0.5422216523E-03, 
0.1888478444E-02, -0.3382973230E-02, -0.5886912723E-03, 0.3289309542E-03, 
0.1750384584E-02, -0.9780261437E-03, -0.1931250902E-02, 0.1079085355E-02, 
-0.3621632497E-03, 0.2023585594E-03, 0.2085884179E-02, -0.1165486300E-02, 
0.6409242844E-03, -0.3581158117E-03, -0.2736289413E-02, 0.1528899888E-02, 
-0.3368917714E-02, 0.1882379797E-02, 0.1038171472E-02, -0.5800785496E-03, 
0.3193936333E-01, -0.1508848895E-01, 0.6823471193E-02, -0.3370611049E-02, 
-0.1435014917E+00, 0.4925184726E-01, -0.6098305316E-01, 0.2705640950E-01, 
0.4814477634E+00, -0.5199195925E-01, 0.7771077176E+00, -0.3643916717E+00, 
0.3644295324E+00, 0.7771916813E+00, -0.5192289139E-01, -0.4813707824E+00, 
-0.2721779258E-01, -0.6127396434E-01, 0.4913249792E-01, 0.1432965635E+00, 
0.3808651449E-02, 0.7607535473E-02, -0.1495173150E-01, -0.3169534309E-01, 
-0.3029170562E-03, -0.5421338829E-03, 0.1889928755E-02, 0.3382425628E-02
};

const double RightPre[] = {
0.1089843053E+01, 0.9175633109E+00, 0.0000000000E+00, 0.0000000000E+00, 
-0.8008132342E+00, 0.3505220326E+00, 0.2096292884E+01, 0.4770325785E+00, 
0.1055782528E+01, 0.9471647557E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.4658798302E+00, 0.2539461793E+00, 
0.1737631796E+01, 0.5754959149E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.2417783151E+01, -0.1059694498E+00, -0.6663367942E+01, 0.6796520006E+00, 
0.5642212524E+01, 0.1772354366E+00, 0.1000398078E+01, 0.9996020804E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.2241154296E-02, 0.2235092825E-02, 
0.1002312956E+01, 0.9976923812E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.1844504727E-01, 0.2335082980E-01, 
0.9170462790E-01, -0.1171759021E+00, 0.7808176166E+00, 0.1280708809E+01, 
0.0000000000E+00, 0.0000000000E+00, -0.7373304912E-02, -0.2846460011E-02, 
-0.9310068502E-03, 0.9005357890E-01, 0.3767386370E+00, -0.9639839213E+00, 
0.5005192311E+00, 0.1997925230E+01, 0.1000413355E+01, 0.9995868157E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
-0.4572991430E-02, 0.4514077458E-02, 0.1012632589E+01, 0.9875250029E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.9186034191E-02, -0.8522720044E-02, 
-0.5149118250E-01, 0.4842233071E-01, 0.1050111166E+01, 0.9522801325E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
-0.8653411835E-01, 0.4412830877E-01, 0.4462922645E+00, -0.2280710782E+00, 
-0.9332961770E+00, 0.5124746278E+00, 0.1734250554E+01, 0.5766179504E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.6651201997E-01, 0.2197171709E+00, 
0.3272287057E+00, -0.1071279051E+01, -0.6835308959E+00, 0.2403809220E+01, 
0.1157116986E+00, -0.2710986078E+00, 0.2461150317E+00, 0.4063140690E+01, 
0.1000078420E+01, 0.9999215862E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
-0.9850976530E-03, 0.9819736622E-03, 0.1003102675E+01, 0.9969069215E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.2633468462E-02, -0.2560290665E-02, -0.1664138984E-01, 0.1623091263E-01, 
0.1022118540E+01, 0.9783601031E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
-0.1429038400E-02, 0.1388671224E-02, 0.9368484963E-02, -0.9150173851E-02, 
-0.2392821028E-01, 0.2393097496E-01, 0.9782470759E+00, 0.1022236636E+01, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1195277199E-01, -0.1168953970E-01, -0.6711312742E-01, 0.6542835502E-01, 
0.1479375368E+00, -0.1465646909E+00, -0.7026196422E-01, 0.7358661282E-01, 
0.9760519092E+00, 0.1024535673E+01, 0.0000000000E+00, 0.0000000000E+00, 
-0.1248032832E-01, 0.2141640226E-01, 0.7502310126E-01, -0.1275021840E+00, 
-0.1804736950E+00, 0.3095551764E+00, 0.1461500318E+00, -0.2545824950E+00, 
0.1504596966E+00, -0.2445553979E+00, 0.6303329548E+00, 0.1586463142E+01, 
0.1000408431E+01, 0.9995917354E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.5293682422E-02, 0.5202484794E-02, 
0.1017114208E+01, 0.9831737600E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1591493727E-01, -0.1341009335E-01, -0.1031512206E+00, 0.8847321480E-01, 
0.1146285615E+01, 0.8723829270E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.3646661500E-01, 0.1680436654E-01, 
0.2406286152E+00, -0.1136375469E+00, -0.6581590373E+00, 0.3658356512E+00, 
0.1569466249E+01, 0.6371592893E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.8620589496E-01, -0.7264200771E-02, -0.5839504132E+00, 0.5536106388E-01, 
0.1663842388E+01, -0.2090053065E+00, -0.2536273536E+01, 0.6450028470E+00, 
0.2505431179E+01, 0.3991328951E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.3341753604E+00, -0.3790478410E-02, 
0.2298323873E+01, 0.1873956341E-01, -0.6712323867E+01, -0.2295905912E-01, 
0.1070192046E+02, -0.9239613787E-01, -0.9843372614E+01, 0.7726440127E+00, 
0.5084895171E+01, 0.1966608881E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.1653022182E+00, 0.3026314320E-01, -0.1367863363E+01, -0.1720350030E+00, 
0.4921932660E+01, 0.4114532707E+00, -0.1002831667E+02, -0.6057060145E+00, 
0.1256338987E+02, 0.5848923814E+00, -0.9758759367E+01, 0.4444537862E+00, 
0.4318033378E+01, 0.2315869083E+00, 0.1000011417E+01, 0.9999885829E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.1885821248E-03, 0.1884331061E-03, 
0.1000779405E+01, 0.9992212023E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.6980138374E-03, -0.6895681425E-03, 
-0.5775869059E-02, 0.5710508693E-02, 0.1010657918E+01, 0.9894544755E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.1587420489E-02, 0.1473040141E-02, 
0.1315517251E-01, -0.1221739122E-01, -0.4838107469E-01, 0.4544822540E-01, 
0.1053305612E+01, 0.9493920742E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.2302522936E-02, -0.1832761923E-02, 
-0.1919511889E-01, 0.1530411869E-01, 0.7113262876E-01, -0.5742370539E-01, 
-0.1541810491E+00, 0.1326319561E+00, 0.1103642518E+01, 0.9060904992E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.6297296878E-02, 0.3996745292E-02, 
0.5008238228E-01, -0.3139508197E-01, -0.1744341513E+00, 0.1085741195E+00, 
0.3478392588E+00, -0.2247406296E+00, -0.4139824420E+00, 0.3061849704E+00, 
0.1225094612E+01, 0.8162634872E+00, 0.0000000000E+00, 0.0000000000E+00, 
0.0000000000E+00, 0.0000000000E+00, 0.1205751485E-01, -0.5357488581E-02, 
-0.9775783312E-01, 0.4355462708E-01, 0.3469638042E+00, -0.1557775571E+00, 
-0.7026319735E+00, 0.3312968007E+00, 0.8627444834E+00, -0.4752264599E+00, 
-0.6802055242E+00, 0.4601200251E+00, 0.1206700215E+01, 0.8287062417E+00, 
0.0000000000E+00, 0.0000000000E+00, -0.1638857751E-01, 0.1062802918E-01, 
0.1318974327E+00, -0.8500762499E-01, -0.4651165383E+00, 0.2996885100E+00, 
0.9374435384E+00, -0.6305890566E+00, -0.1161167601E+01, 0.9184100559E+00, 
0.8600562361E+00, -0.8062512474E+00, -0.2230327919E+00, 0.2486081011E+00, 
0.7434539179E+00, 0.1345073280E+01
};

/* The following is in Fryzlewicz's ``WavInt.c'' */

void Precondition(int Scale, int Direction, Filter F, double* Vect) {

double* templ;
double* tempr;

if (F.Length >= 3) {

int length, step, i, j;

length = F.Length / 2;
step = (int) pow(2.0, Scale);

templ = (double*) malloc(length * sizeof(double));
tempr = (double*) malloc(length * sizeof(double)); 

for (i = 0; i < length; i++) {
templ[i] = tempr[i] = 0;
switch (Direction) {
case NORMAL:
for (j = 0; j < length; j++) {
templ[i] += Vect[j] * F.PreLeft[i][j];
tempr[i] += Vect[step - length + j] * F.PreRight[i][j];
}
break;
case INVERSE:
for (j = 0; j < length; j++) {
templ[i] += Vect[j] * F.PreInvLeft[i][j];
tempr[i] += Vect[step - length + j] * F.PreInvRight[i][j];
}
}
}

for (i = 0; i < length; i++) {
Vect[i] = templ[i];
Vect[step - length + i] = tempr[i];
}

free(templ);
free(tempr);

}
}

void TransStep(int Scale, Filter F, double* Vect) {

double* temp;
int length, halflength, N, pos, i, j, p;

length = (int) pow(2.0, Scale);
halflength = length / 2;
N = F.Length / 2;

temp = (double*) malloc(length * sizeof(double));

if (N > 1) {

pos = 0;

for (i = 0; i < N; i++) {
p = 2 * i;
temp[pos] = temp[pos + halflength] = 0;
for (j = 0; j <= N + p; j++) {
temp[pos] += Vect[j] * F.HLeft[i][j];
temp[pos + halflength] += Vect[j] * F.GLeft[i][j];
}
pos++;
}

for (i = N; i < halflength - N; i++) {
p = 2 * i - N + 1;
temp[pos] = temp[pos + halflength] = 0;
for (j = 0; j < 2 * N; j++) {
temp[pos] += Vect[p + j] * F.H[j];
temp[pos + halflength] += Vect[p + j] * F.G[j];
}
pos++;
}

for (i = N - 1; i >= 0; i--) {
p = 2 * i;
temp[pos] = temp[pos + halflength] = 0;
for (j = 0; j <= N + p; j++) {
temp[pos] += Vect[length - j - 1] * F.HRight[i][j];
temp[pos + halflength] += Vect[length - j - 1] * F.GRight[i][j];
}
pos++;
}

}

else

for (i = 0; i < halflength; i++) {
p = 2 * i;
temp[i] = temp[i + halflength] = 0;
for (j = 0; j < 2 * N; j++) {
temp[i] += Vect[p + j] * F.H[j];
temp[i + halflength] += Vect[p + j] * F.G[j];
}
}

for(i = 0; i < length; i++) Vect[i] = temp[i];

free(temp);

}

void InvTransStep(int Scale, Filter F, double* Vect) {

double* temp;
int length, doublelength, N, pos, i, j, p;

length = (int) pow(2.0, Scale);
doublelength = 2 * length;
N = F.Length / 2;

temp = (double*) malloc(doublelength * sizeof(double));
for (i = 0; i < doublelength; i++) temp[i] = 0;

if (N > 1) {

pos = 0;

for (i = 0; i < N; i++) {
p = 2 * i;
for (j = 0; j <= N + p; j++) {
temp[j] += Vect[pos] * F.HLeft[i][j];
temp[j] += Vect[pos + length] * F.GLeft[i][j];
}
pos++;
}

for (i = N; i < length - N; i++) {
p = 2 * i - N + 1;
for (j = 0; j < 2 * N; j++) {
temp[p + j] += Vect[pos] * F.H[j];
temp[p + j] += Vect[pos + length] * F.G[j];
}
pos++;
}

for (i = N - 1; i >= 0; i--) {
p = 2 * i;
for (j = 0; j <= N + p; j++) {
temp[doublelength - j - 1] += Vect[pos] * F.HRight[i][j];
temp[doublelength - j - 1] += Vect[pos + length] * F.GRight[i][j];
}
pos++;
}

} 

else

for (i = 0; i < length; i++) {
p = 2 * i;
for (j = 0; j < 2; j++) {
temp[p + j] += Vect[i] * F.H[j];
temp[p + j] += Vect[i + length] * F.G[j];
}
}

for (i = 0; i < doublelength; i++) Vect[i] = temp[i];

free(temp);

}

void Trans(int MinScale, int Direction, int FilterNumber, double* Vect, int Size, int
Precond, int* FilterHistory) {

int scale, maxscale, N, NPrev, NNext;
Filter temp, temp1;

maxscale = (int)(log(Size)/log(2));

if (MinScale >= maxscale) {
Rprintf("MinScale must be less than log2(Size).\nNo transformation performed.\n");
return;
}

if (FilterNumber < 1 || FilterNumber > Nmax) {
Rprintf("Filter no %d not implemented.\nNo transformation performed.\n", FilterNumber);
return;
}

N = FilterNumber;

if (Direction == NORMAL)
    for (scale = maxscale; scale > MinScale; scale--) {
        NPrev = N;
        while (((int)pow(2.0, scale)) < 8*N && N != 1) N--;
        FilterHistory[maxscale - scale] = N;
        temp = GetFilt(N);
        if (Precond){           /* MAN: added brace for unambiguity */
            if (scale == maxscale) Precondition(scale, NORMAL, temp, Vect);
            else if (N != NPrev) {
                temp1 = GetFilt(NPrev);
                Precondition(scale, INVERSE, temp1, Vect);
                Precondition(scale, NORMAL, temp, Vect);
            }
        }    /* MAN: added brace for unambiguity */
        TransStep(scale, temp, Vect);
    }

else {
    while (((int)pow(2.0, MinScale+1)) < 8*N && N != 1) N--;
    for (scale = MinScale; scale < maxscale; scale++) {
        N = FilterHistory[maxscale - scale - 1];
        if (scale < maxscale - 1) NNext = FilterHistory[maxscale - scale - 2];
            else NNext = N;     
        temp = GetFilt(N);
        InvTransStep(scale, temp, Vect);
        if (Precond){       /* MAN: added for unambiguity */
            if (scale + 1 == maxscale) Precondition(maxscale, INVERSE, temp,
Vect);
            else if (N != NNext) {
                temp1 = GetFilt(NNext);
                Precondition(scale+1, INVERSE, temp, Vect);
                Precondition(scale+1, NORMAL, temp1, Vect);
            }
        }       /* MAN: added for unambiguity */    
    }
}

}

/* The following is in Fryzlewicz's WavIntC.c */

void dec(double* data, int* size, int* filternumber, int* minscale, int* precond,       /* MAN: added missing void fn type */
int* filterhistory) {

Trans(*minscale, NORMAL, *filternumber, data, *size, *precond, filterhistory);

}

void rec(double* data, int* size, int* filterhistory, int* currentscale, int* precond) {    /* MAN: added missing void fn type */

Trans(*currentscale, INVERSE, filterhistory[0], data, *size, *precond, filterhistory);

}

/* The following in Fryzlewicz's ``Filters.c'' */

Filter GetFilt(int N) {

Filter temp;
int i, j, len, offset, offset1;
double NormH, NormHR, NormHL, NormGL, NormGR;

temp.Length = 0;
for (i = 0; i < 2 * Nmax; i++) temp.H[i] = temp.G[i] = 0;
for (i = 0; i < Nmax; i++)
for (j = 0; j < 3 * Nmax -1; j++)
temp.HLeft[i][j] = temp.GLeft[i][j] = temp.HRight[i][j] = temp.GRight[i][j] = 0;
for (i = 0; i < Nmax; i++)
for (j = 0; j < Nmax; j++)
temp.PreLeft[i][j] = temp.PreInvLeft[i][j] = temp.PreRight[i][j] = temp.PreInvRight[i][j] = 0;

if (N < 1 || N > Nmax) {
Rprintf("Filter no %d not implemented.", N);
return temp;
}

temp.Length = 2 * N;

/* Interior */

offset = 0;
len = 2 * N;

for (i = 1; i < N; i++) offset += 2 * i;

for (i = 0; i < len; i++) temp.H[i] = Interior[i + offset];
NormH = Sum(temp.H, len);
for (i = 0; i < len; i++) temp.H[i] = temp.H[i] / NormH * sqrt(2.0);

for (i = 0; i < len; i++) temp.G[i] = (-2 * (i % 2) + 1) * temp.H[len - i - 1];

/* Left and Right */

offset = offset1 = 0;

for (i = 1; i < N; i++) offset += 4 * i * i;

for (i = 0; i < N; i++) {
len = N + 2 * i + 1;
NormHL = 0.0;
NormGL = 0.0;
NormHR = 0.0;
NormGR = 0.0;
for (j = 0; j < len; j++) {
temp.HLeft[i][j] = Left[offset + offset1 + 2 * j];
NormHL += pow(temp.HLeft[i][j], 2.0);
temp.GLeft[i][j] = Left[offset + offset1 + 2 * j + 1];
NormGL += pow(temp.GLeft[i][j], 2.0);
temp.HRight[i][j] = Right[offset + offset1 + 2 * j];
NormHR += pow(temp.HRight[i][j], 2.0);
temp.GRight[i][j] = Right[offset + offset1 + 2 * j + 1];
NormGR += pow(temp.GRight[i][j], 2.0);
}
for (j = 0; j < len; j++) {
    temp.HLeft[i][j] /= sqrt(NormHL);
    temp.GLeft[i][j] /= sqrt(NormGL);
    temp.HRight[i][j] /= sqrt(NormHR);
    temp.GRight[i][j] /= sqrt(NormGR);
}
offset1 += 2 * len;
}

/* Preconditioning Matrices: Left and Right */

if (N > 1) {

offset = 0;

for (i = 2; i < N; i++) offset += 2 * i * i;

for (i = 0; i < N; i++)
for (j = 0; j < N; j++) {
offset1 = 2 * N * i + 2 * j;
temp.PreLeft[i][j] = LeftPre[offset + offset1];
temp.PreInvLeft[i][j] = LeftPre[offset + offset1 + 1];
temp.PreRight[i][j] = RightPre[offset + offset1];
temp.PreInvRight[i][j] = RightPre[offset + offset1 + 1];
}
}
return temp;
}

double Sum(double* vect, int length) {
    double ssum;
    int i;
    ssum = 0.0;
    for (i = 0; i < length; i++) ssum += vect[i];
    return ssum;
}

/*
 * ThreeD wavelets suite: three-dimensional DWT and inverse
 *
 * A generic 3D array has nr rows, nc cols and ns pixels in the sides.
 * We use the letters r,c and s to index each type.
 */


/* Macro to access 3D array */

#define ACCESSW3D(array, nr, nc, r, c, s) *(array + (nr)*((c) + (s)*(nc))+(r))

/*
 * CreateArray3D:   Create a 3D array of doubles
 *
 *  Arguments.
 *
 *      nr: number of rows      (integer)
 *      nc: number of columns   (integer)
 *      ns: number of sides     (integer)
 *      error:  error code. 0 is o.k., 3001 is memory error.
 *
 *  Returns:    NULL on error or pointer to requested array.
 *
 */


double *CreateArray3D(int nr, int nc, int ns, int *error)
{
double *array;

*error = 0;


if ((array = (double *)malloc((unsigned)(nr*nc*ns)*sizeof(double)))==NULL){
    *error = 3001;
    return(NULL);
    }

else
    return(array);
}

/*
 * DestroyArray3D:  Release memory associated with an array 
 *
 *  Arguments.
 *
 *      array:  pointer to 3D array
 *      error:  error code. 0=O.k. 3002 means NULL pointer was passed.
 *
 *  Returns:    NULL on error or pointer to requested array.
 *
 */

void DestroyArray3D(double *array, int *error)
{


*error = 0;

if (array == NULL)  {
    *error = 3002;
    return;
    }
else
    free((void *)array);
}

/*
 * wd3Dstep - the guts of the 3D DWT algorithm. This algorithm
 *      could be made more efficient by less memory allocation
 *      but I ain't got time to do it
 */


void wd3Dstep(double *Carray, int *truesize, int *size, double *H,
	int *LengthH, int *error)
/*---------------------
 * Argument description
double *Carray::  Input 3D array. All dimensions are size      
int *truesize::   The true dimensions of the Carray            
int *size::       Number of rows, columns and sides (power of 2)   
                  For this invocation of the routine only      
double *H::       Wavelet filter coefficients              
int *LengthH::    Number of wavelet filter coefficients        
int *error::      Error code.  0=O.k.                  
         Memory errors 3003 to 3017               
 *---------------------*/
{
register int r,c,s; /* Counters for rows, cols and sides        */
double *Ha,*Ga; /* Will be storage for first application of filters */

double *HH, *GH, *HG, *GG;  /* Storage for second application   */

double *HHH,*GHH,*HGH,*GGH,*HHG,*GHG,*HGG,*GGG; /* Third application    */

int ndata;          /* Length of TheData            */

int halfsize;

double *c_in, *c_out, *d_out;

/* Creates a 3D array           */
double *CreateArray3D(int nr, int nc, int ns, int *error);

void convolveC(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *c_out, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void convolveD(double *c_in, int LengthCin, int firstCin,
	double *H, int LengthH,
	double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);



*error = 0;

halfsize = *size/2;


/*
 * Get memory for first application
 */

if ((Ha = CreateArray3D((int)halfsize, (int)*size, (int)*size, error))==NULL){
    return;
    }
if ((Ga = CreateArray3D((int)halfsize, (int)*size, (int)*size, error))==NULL){
    return;
    }

/* Get some storage for c_in, c_out, d_out */

ndata = *size;
if ((c_in = (double *)malloc((unsigned)ndata*sizeof(double)))==NULL){
    *error = 3003;
    return;
    }
if ((c_out = (double *)malloc((unsigned)(int)halfsize*sizeof(double)))==NULL){
    *error = 3004;
    return;
    }
if ((d_out = (double *)malloc((unsigned)(int)halfsize*sizeof(double)))==NULL){
    *error = 3005;
    return;
    }

/*
 * Now perform wavelet transform across rows for each column and side
 */

for(c=0; c< *size; ++c)
    for(s=0; s < *size; ++s)    {

        /* Load up c_in array */

        for(r=0; r < *size; ++r)    {
            *(c_in+r) = ACCESSW3D(Carray, (int)*truesize,
                (int)*truesize, r, c, s);
            /*
            Rprintf("Carray[%d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(Carray,
                (int)*truesize, (int)*truesize,
                                r, c, s));
            */
            }


        /* Now do convolveC and convolveD on c_in */

        convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
            c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

        convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
            d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);



        /* Now store C in Ha and D in Ga */

        for(r=0; r < (int)halfsize; ++r)    {
            ACCESSW3D(Ha, (int)halfsize, (int)*size, r, c, s) =
                *(c_out+r);

            /*
            Rprintf("Ha[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(Ha, (int)halfsize,
                    (int)*size, r,c,s));
            */
                }
        for(r=0; r < (int)halfsize; ++r)    {
            ACCESSW3D(Ga, (int)halfsize, (int)*size, r, c, s) =
                *(d_out+r);
            /*
            Rprintf("Ga[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(Ga, (int)halfsize,
                    (int)*size, r,c,s));
            */
            }

        /* Go round and do it again */
        }


/* Now create memory for second application */

if ((HH = CreateArray3D((int)halfsize, (int)halfsize, (int)*size, error))==NULL){
    *error = 3006;
    return;
    }
if ((GH = CreateArray3D((int)halfsize, (int)halfsize, (int)*size, error))==NULL){
    *error = 3007;
    return;
    }
if ((HG = CreateArray3D((int)halfsize, (int)halfsize, (int)*size, error))==NULL){
    *error = 3008;
    return;
    }
if ((GG = CreateArray3D((int)halfsize, (int)halfsize, (int)*size, error))==NULL){
    *error = 3009;
    return;
    }

/* Ha to HH and GH */

/*
 * Now perform convolution steps over cols on H for each row and side. 
 */

for(r=0; r < halfsize; ++r)
    for(s=0; s < *size; ++s)    {

        /* Load up c_in array */

        for(c=0; c < *size; ++c)
            *(c_in+c) = ACCESSW3D(Ha, (int)halfsize, (int)*size,
                r, c, s);



        /* Now do convolveC and convolveD on c_in */

        convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
            c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

        convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
            d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);


        /* Now store C in HH and D in GH */

        for(c=0; c < (int)halfsize; ++c)    {
            ACCESSW3D(HH, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+c);
            /*
            Rprintf("HH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
                }
        for(c=0; c < (int)halfsize; ++c)    {

            ACCESSW3D(GH, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+c);
            /*
            Rprintf("GH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(GH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
            }

        /* Go round and do it again */
        }

/* Ga to HG and GG */

/*
 * Now perform convolution over cols on G for each row and side. 
 */

for(r=0; r < halfsize; ++r)
    for(s=0; s < *size; ++s)    {

        /* Load up c_in array */

        for(c=0; c < *size; ++c)
            *(c_in+c) = ACCESSW3D(Ga, (int)halfsize, (int)*size,
                r, c, s);


        /* Now do convolveC and convolveD on c_in */ 

        convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
                        c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

                convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
                        d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);


        /* Now store C in HG and D in GG */

        for(c=0; c < (int)halfsize; ++c)    {
            ACCESSW3D(HG, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+c);

            /*
            Rprintf("HG[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HG, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
                }
        for(c=0; c < (int)halfsize; ++c)    {
            ACCESSW3D(GG, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+c);
            /*
            Rprintf("GG[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HG, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
            }

        /* Go round and do it again */
        }

/* Now we've used Ha and Ga and so we can free them */

free((void *)Ha);
free((void *)Ga);

/* THIRD LEVEL APPLICATION  */

/* Now create memory for third application */

if ((HHH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3010;
    return;
    }
if ((GHH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3011;
    return;
    }
if ((HGH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3012;
    return;
    }
if ((GGH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3013;
    return;
    }
if ((HHG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3014;
    return;
        }
if ((GHG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3015;
    return;
    }
if ((HGG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3016;
    return;
    }
if ((GGG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3017;
    return;
    }

/* HH to HHH and GHH */

/*
 * Now perform wavelet transform over sides on HH for each row and col. 
 */

for(r=0; r < halfsize; ++r)
    for(c=0; c < halfsize; ++c) {

        /* Load up c_in array */

        for(s=0; s < *size; ++s)
            *(c_in+s) = ACCESSW3D(HH, (int)halfsize,
                (int)halfsize, r, c, s);


        /* Now do convolveC and convolveD on c_in */

                convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
                        c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

                convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
                        d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);


        /* Now store C in HHH and D in GHH */

        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(HHH, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
            /*
            Rprintf("HHH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HHH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
                }
        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(GHH, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+s);
            }

        /* Go round and do it again */
        }

/* GH to HGH and GGH */

/*
 * Now perform wavelet transform over sides on GH for each row and col. 
 */

for(r=0; r < halfsize; ++r)
    for(c=0; c < halfsize; ++c) {

        /* Load up c_in array */

        for(s=0; s < *size; ++s)
            *(c_in+s) = ACCESSW3D(GH, (int)halfsize,
                (int)halfsize, r, c, s);


        /* Now do convolveC and convolveD on c_in */

                convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
                        c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

                convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
                        d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

        /* Now store C in HGH and D in GGH */

        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(HGH, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
            /*
            Rprintf("HGH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HGH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
                }
        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(GGH, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+s);
            /*
            Rprintf("GGH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(GGH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
            }

        /* Go round and do it again */
        }

/* HG to HHG and GHG */

/*
 * Now perform wavelet transform over sides on HG for each row and col. 
 */

for(r=0; r < halfsize; ++r)
    for(c=0; c < halfsize; ++c) {

        /* Load up c_in array */

        for(s=0; s < *size; ++s)
            *(c_in+s) = ACCESSW3D(HG, (int)halfsize,
                (int)halfsize, r, c, s);


        /* Now do convolveC and convolveD on c_in */

                convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
                        c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

                convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
                        d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

        /* Now store C in HHG and D in GHG */

        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(HHG, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
                }
        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(GHG, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+s);
            }

        /* Go round and do it again */
        }

/* GG to HGG and GGG */

/*
 * Now perform wavelet transform over sides on GG for each row and col. 
 */

for(r=0; r < halfsize; ++r)
    for(c=0; c < halfsize; ++c) {

        /* Load up c_in array */

        for(s=0; s < *size; ++s)
            *(c_in+s) = ACCESSW3D(GG, (int)halfsize,
                (int)halfsize, r, c, s);


        /* Now do convolveC and convolveD on c_in */

                convolveC(c_in, (int)*size, 0, H, (int)*LengthH,
                        c_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

                convolveD(c_in, (int)*size, 0, H, (int)*LengthH,
                        d_out, 0, (int)halfsize-1, WAVELET, 1, PERIODIC);

        /* Now store C in HGG and D in GGG */

        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(HGG, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
                }
        for(s=0; s < (int)halfsize; ++s)    {
            ACCESSW3D(GGG, (int)halfsize, (int)halfsize, r, c, s) =
                *(d_out+s);
            }

        /* Go round and do it again */
        }

/* Now we can get rid of the second level memory */

free((void *)HH);
free((void *)GH);
free((void *)HG);
free((void *)GG);


/* Now store the answers in the C array             */

/* HHH */

for(r=0; r < (int)halfsize; ++r)
  for(c=0; c < (int)halfsize; ++c)
    for(s=0; s < (int)halfsize; ++s)    {

    /* HHH */

    /* { double tmpf;

    tmpf= ACCESSW3D(Carray, (int)*truesize, (int)*truesize, r, c, s) =
	ACCESSW3D(HHH, (int)halfsize, (int)halfsize, r, c, s); */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize, r, c, s) =
        ACCESSW3D(HHH, (int)halfsize, (int)halfsize, r, c, s); 

    /*
    Rprintf("Carray[%d %d %d] = %lf (from HHH)\n",
        r,c,s,tmpf);
    */
    /* } */

    /* GHH */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c, s) =
        ACCESSW3D(GHH, (int)halfsize, (int)halfsize, r, c, s);
    
    /* HGH */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r, c+(int)halfsize, s) =
        ACCESSW3D(HGH, (int)halfsize, (int)halfsize, r, c, s);

    /*
    Rprintf("HGH[ %d, %d, %d] = %lf\n", r,c,s, ACCESSW3D(HGH, (int)halfsize,
            (int)halfsize, r,c,s));
    */

    /* GGH */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c+(int)halfsize, s) =
        ACCESSW3D(GGH, (int)halfsize, (int)halfsize, r, c, s);

    /*
    Rprintf("GGH[ %d, %d, %d] = %lf\n", r,c,s, ACCESSW3D(GGH, (int)halfsize,
            (int)halfsize, r,c,s));
    */

    /* HHG */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize, r, c, s+(int)halfsize)
        =
        ACCESSW3D(HHG, (int)halfsize, (int)halfsize, r, c, s);

    /* GHG */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c, s+(int)halfsize) =
        ACCESSW3D(GHG, (int)halfsize, (int)halfsize, r, c, s);
    
    /* HGG */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r, c+(int)halfsize, s+(int)halfsize) =
        ACCESSW3D(HGG, (int)halfsize, (int)halfsize, r, c, s);

    /* GGG */

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c+(int)halfsize, s+(int)halfsize) =
        ACCESSW3D(GGG, (int)halfsize, (int)halfsize, r, c, s);

    }
    

/* Free the third level memory */

free((void *)HHH);
free((void *)GHH);
free((void *)HGH);
free((void *)GGH);
free((void *)HHG);
free((void *)GHG);
free((void *)HGG);
free((void *)GGG);


/* Free c_in, c_out, d_out */

free((void *)c_in);
free((void *)c_out);
free((void *)d_out);

}


void wd3D(double *Carray, int *size, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *Carray::      Input and output coefficients        
int *size::           Dimension of this array          
double *H::           The wavelet coefficients         
int *LengthH::        Number of wavelet coefficients       
int *error::          Error code 0=o.k.                
 *---------------------*/
{
int insize;
void wd3Dstep(double *Carray, int *truesize, int *size, double *H,
	int *LengthH, int *error);

*error = 0;

insize = *size;

while(insize >= 2)  {

    /*
    Rprintf("Outsize is %ld\n", insize*2);
    */

    wd3Dstep(Carray, size, &insize, H, LengthH, error);

    if (*error != 0)
        return;

    insize /= 2;
    }

}

/*
 * Reconstruct 3D wavelet object in Carray
 */

void wr3D(double *Carray, int *truesize, double *H, int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *Carray::  Contains array of wavelet coefficients       
int *truesize::   Dimension of 3D array Carray             
double *H::       The wavelet filter coefficients          
int *LengthH::    Number of wavelet filter coefficients        
int *error::      Error code   (0=o.k.)                
         	  Memory errors from wr3Dstep              
         	  3035l the dimension of Carray is 1, therefore cannot 
		  do any further reconstruction           
 *---------------------*/
{
int sizeout;

void wr3Dstep(double *Carray, int *truesize, int *sizeout,
	double *H, int *LengthH, int *error);

*error = 0;

sizeout = 2;

if (*truesize < sizeout)    {
    *error = 3035;
    return;
    }

while(sizeout <= *truesize) {

    /*
    Rprintf("Outsize is %ld\n", sizeout);
    */

    wr3Dstep(Carray, truesize, &sizeout, H, LengthH, error);

    if (*error != 0)
        return;

    sizeout *= 2;
    }
}


/*
 * wr3Dstep:    Perform 3D wavelet reconstruction step
 */

void wr3Dstep(double *Carray, int *truesize, int *sizeout,
	double *H, int *LengthH, int *error) 
/*---------------------
 * Argument description
 *---------------------
double *Carray::  Array of wavelet coefficients and previous Cs to replace 
int *truesize::   True size of Carray                            
int *sizeout::    Size of answer array                     
double *H::       The wavelet coefficients                 
int *LengthH::    Number of wavelet coefficients               
int *error::      Error code. 0=o.k.
                  Memory errors 3018 to 3034
 *---------------------*/
{
register int r,c,s;
double *Ha,*Ga; /* Will be storage for third application of filters */

double *HH, *GH, *HG, *GG;  /* Storage for second application   */

double *HHH,*GHH,*HGH,*GGH,*HHG,*GHG,*HGG,*GGG; /* Third application    */
double *c_in, *d_in, *c_out;
int halfsize;
int type,bc;

void conbar(double *c_in, int LengthCin, int firstCin,
	double *d_in, int LengthDin, int firstDin,
	double *H, int LengthH,
	double *c_out, int LengthCout, int firstCout, int lastCout,
	int type, int bc);

*error = 0;
type = WAVELET;
bc = PERIODIC;

/*
 * Take the coefficients from the C array and store them in cubes
 * half the sizeout
 *
 * This is just the inverse of the last part of wd3Dstep
 */

halfsize = *sizeout/2;

/* Now create memory for first application */ 

if ((HHH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3018;
    return;
    }
if ((GHH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3019;
    return;
    }
if ((HGH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3020;
    return;
    }
if ((GGH = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3021;
    return;
    }
if ((HHG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3022;
    return;
    }
if ((GHG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3023;
    return;
    }
if ((HGG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3024;
    return;
    }
if ((GGG = CreateArray3D((int)halfsize, (int)halfsize, (int)halfsize, error))==NULL){
    *error = 3025;
    return;
    }

/* HHH */

for(r=0; r < (int)halfsize; ++r)
  for(c=0; c < (int)halfsize; ++c)
    for(s=0; s < (int)halfsize; ++s)    {

    /* HHH */

    ACCESSW3D(HHH, (int)halfsize, (int)halfsize, r, c, s) =
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize, r, c, s);

    /* GHH */

    ACCESSW3D(GHH, (int)halfsize, (int)halfsize, r, c, s) =
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
            r+(int)halfsize, c, s);
    
    /* HGH */

    ACCESSW3D(HGH, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
            r, c+(int)halfsize, s);
    /*
    Rprintf("HGH[ %d, %d, %d] = %lf\n", r,c,s, ACCESSW3D(HGH, (int)halfsize,
            (int)halfsize, r,c,s));
    */

    /* GGH */

    ACCESSW3D(GGH, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
            r+(int)halfsize, c+(int)halfsize, s);

    /*
    Rprintf("GGH[ %d, %d, %d] = %lf\n", r,c,s, ACCESSW3D(GGH, (int)halfsize,
            (int)halfsize, r,c,s));
    */

    /* HHG */

    ACCESSW3D(HHG, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize, r, c,
            s+(int)halfsize);

    /* GHG */

    ACCESSW3D(GHG, (int)halfsize, (int)halfsize, r, c, s) =
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
            r+(int)halfsize, c, s+(int)halfsize);
    
    /* HGG */

    ACCESSW3D(HGG, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
            r, c+(int)halfsize, s+(int)halfsize);

    /* GGG */

    ACCESSW3D(GGG, (int)halfsize, (int)halfsize, r, c, s) =
    ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c+(int)halfsize, s+(int)halfsize);

    }

/* Now create memory for HH, GH, HG and GG */

if ((HH = CreateArray3D((int)halfsize, (int)halfsize, (int)*sizeout, error))
        ==NULL){
    *error = 3026;
    return;
    }
if ((GH = CreateArray3D((int)halfsize, (int)halfsize, (int)*sizeout, error))
        ==NULL){
    *error = 3027;
    return;
    }
if ((HG = CreateArray3D((int)halfsize, (int)halfsize, (int)*sizeout, error))
        ==NULL){
    *error = 3028;
    return;
    }
if ((GG = CreateArray3D((int)halfsize, (int)halfsize, (int)*sizeout, error))
        ==NULL){
    *error = 3029;
    return;
    }

/* We now have to reconstruct HH, GH, HG and GG */

/* Create c_in, d_in and c_out */

if ((c_in = (double *)malloc((unsigned)(int)halfsize*sizeof(double)))==NULL){
    *error = 3030;
    return;
    }

if ((d_in = (double *)malloc((unsigned)(int)halfsize*sizeof(double)))==NULL){
    *error = 3031;
    return;
    }

if ((c_out = (double *)malloc((unsigned)(int)*sizeout*sizeof(double)))==NULL){
    *error = 3032;
    return;
    }


/* Fill up HH by wavelet reconstruction of HHH and GHH */

for(r=0; r < (int)halfsize; ++r)
    for(c=0; c < (int)halfsize; ++c)    {

        /* Fill up c_in and d_in */ 

        for(s=0; s < (int)halfsize; ++s)    {
            *(c_in+s) = ACCESSW3D(HHH, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+s) = ACCESSW3D(GHH, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up HH */

        for(s=0; s < (int)*sizeout; ++s)    {
            ACCESSW3D(HH, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
            /*
            Rprintf("HH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(HH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
            }


        /* Next row and column combination */
        }

/* Fill up GH by wavelet reconstruction of HGH and GGH */

for(r=0; r < (int)halfsize; ++r)
    for(c=0; c < (int)halfsize; ++c)    {

        /* Fill up c_in and d_in */ 

        for(s=0; s < (int)halfsize; ++s)    {
            *(c_in+s) = ACCESSW3D(HGH, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+s) = ACCESSW3D(GGH, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up GH */

        for(s=0; s < (int)*sizeout; ++s)    {
            ACCESSW3D(GH, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);
            /*
            Rprintf("GH[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(GH, (int)halfsize,
                    (int)halfsize, r,c,s));
            */
        }


        /* Next row and column combination */
        }

/* Fill up HG by wavelet reconstruction of HHG and GHG */

for(r=0; r < (int)halfsize; ++r)
    for(c=0; c < (int)halfsize; ++c)    {

        /* Fill up c_in and d_in */ 

        for(s=0; s < (int)halfsize; ++s)    {
            *(c_in+s) = ACCESSW3D(HHG, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+s) = ACCESSW3D(GHG, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up HG */

        for(s=0; s < (int)*sizeout; ++s)
            ACCESSW3D(HG, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);


        /* Next row and column combination */
        }

/* Fill up GG by wavelet reconstruction of HGG and GGG */

for(r=0; r < (int)halfsize; ++r)
    for(c=0; c < (int)halfsize; ++c)    {

        /* Fill up c_in and d_in */ 

        for(s=0; s < (int)halfsize; ++s)    {
            *(c_in+s) = ACCESSW3D(HGG, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+s) = ACCESSW3D(GGG, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up GG */

        for(s=0; s < (int)*sizeout; ++s)
            ACCESSW3D(GG, (int)halfsize, (int)halfsize, r, c, s) =
                *(c_out+s);


        /* Next row and column combination */
        }

/* Now we can rid ourselves of HHH, GHH, HGH, GGH, HHG, GHG, HGG, and GGG */

free((void *)HHH);
free((void *)GHH);
free((void *)HGH);
free((void *)GGH);
free((void *)HHG);
free((void *)GHG);
free((void *)HGG);
free((void *)GGG);

/* Now create memory for Ha and Ga */

if ((Ha = CreateArray3D((int)halfsize, (int)*sizeout, (int)*sizeout, error))
    ==NULL){
    *error = 3033;
        return;
        }
if ((Ga = CreateArray3D((int)halfsize, (int)*sizeout, (int)*sizeout, error))
    ==NULL){
    *error = 3034;
        return;
        }


/* Fill up Ha by wavelet reconstruction of HH and GH */

for(r=0; r < (int)halfsize; ++r)
    for(s=0; s < (int)*sizeout; ++s)    {

        /* Fill up c_in and d_in */ 

        for(c=0; c < (int)halfsize; ++c)    {
            *(c_in+c) = ACCESSW3D(HH, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+c) = ACCESSW3D(GH, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up Ha */

        for(c=0; c < (int)*sizeout; ++c)    {
            ACCESSW3D(Ha, (int)halfsize, (int)*sizeout, r, c, s) =
                *(c_out+c);
            /*
            Rprintf("Ha[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(Ha, (int)halfsize,
                    (int)*sizeout, r,c,s));
            */
        }


        /* Next row and side combination */
        }

/* Fill up Ga by wavelet reconstruction of HG and GG */

for(r=0; r < (int)halfsize; ++r)
    for(s=0; s < (int)*sizeout; ++s)    {

        /* Fill up c_in and d_in */ 

        for(c=0; c < (int)halfsize; ++c)    {
            *(c_in+c) = ACCESSW3D(HG, (int)halfsize, (int)halfsize,
                r, c, s);
            *(d_in+c) = ACCESSW3D(GG, (int)halfsize, (int)halfsize,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up Ga */

        for(c=0; c < (int)*sizeout; ++c)    {
            ACCESSW3D(Ga, (int)halfsize, (int)*sizeout, r, c, s) =
                *(c_out+c);
            /*
            Rprintf("Ga[ %d, %d, %d] = %lf\n",
                r,c,s, ACCESSW3D(Ga, (int)halfsize,
                    (int)*sizeout, r,c,s));
            */
            }


        /* Next row and side combination */
        }

/* Now rid outselves of the 2nd level memory */

free((void *)HH);
free((void *)GH);
free((void *)HG);
free((void *)GG);

/* Now store the result of combining Ha,Ga back in the Carray */

for(c=0; c < (int)*sizeout; ++c)
    for(s=0; s < (int)*sizeout; ++s)    {

        /* Fill up c_in and d_in */ 

        for(r=0; r < (int)halfsize; ++r)    {
            *(c_in+r) = ACCESSW3D(Ha, (int)halfsize, (int)*sizeout,
                r, c, s);
            *(d_in+r) = ACCESSW3D(Ga, (int)halfsize, (int)*sizeout,
                r, c, s);
                }

        /* Do the wavelet reconstruction step */

        conbar(c_in, (int)halfsize, 0,
               d_in, (int)halfsize, 0,
               H, (int)*LengthH,
               c_out, (int)*sizeout, 0, (int)*sizeout - 1,
            type, bc);


        /* Now fill up Carray */

        for(r=0; r < (int)*sizeout; ++r)
            ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
                r, c, s) =
                *(c_out+r);


        /* Next column and side combination */
        }

/* Free the first level memory */

free((void *)Ha);
free((void *)Ga);

/* Free c_in, c_out, d_out */

free((void *)c_in);
free((void *)c_out);
free((void *)d_in);
}


void getARRel(double *Carray, int *size, int *level, double *GHH, double *HGH, double *GGH, double *HHG, double *GHG, double *HGG, double *GGG)
{
register int r,c,s;
int halfsize;

halfsize = 1 << *level;

/*
Rprintf("Halfsize is %ld\n", halfsize);
*/

for(r=0; r < (int)halfsize; ++r)
  for(c=0; c < (int)halfsize; ++c)
    for(s=0; s < (int)halfsize; ++s)    {
    /* GHH */

    ACCESSW3D(GHH, (int)halfsize, (int)halfsize, r, c, s) =
        ACCESSW3D(Carray, (int)*size, (int)*size, r+(int)halfsize, c, s);
    
    /* HGH */

    ACCESSW3D(HGH, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*size, (int)*size, r, c+(int)halfsize, s);

    /* GGH */

    ACCESSW3D(GGH, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*size, (int)*size,
            r+(int)halfsize, c+(int)halfsize, s);

    /* HHG */

    ACCESSW3D(HHG, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*size, (int)*size, r, c,
            s+(int)halfsize);

    /* GHG */

    ACCESSW3D(GHG, (int)halfsize, (int)halfsize, r, c, s) =
        ACCESSW3D(Carray, (int)*size, (int)*size,
            r+(int)halfsize, c, s+(int)halfsize);
    
    /* HGG */

    ACCESSW3D(HGG, (int)halfsize, (int)halfsize, r, c, s) = 
        ACCESSW3D(Carray, (int)*size, (int)*size,
            r, c+(int)halfsize, s+(int)halfsize);

    /* GGG */

    ACCESSW3D(GGG, (int)halfsize, (int)halfsize, r, c, s) =
    ACCESSW3D(Carray, (int)*size, (int)*size,
        r+(int)halfsize, c+(int)halfsize, s+(int)halfsize);


    }


}

#define IX_HHH  0
#define IX_GHH  1
#define IX_HGH  2
#define IX_GGH  3
#define IX_HHG  4
#define IX_GHG  5
#define IX_HGG  6
#define IX_GGG  7

void putarr(double *Carray, int *truesize, int *level,
	int *Iarrayix, double *Iarray)
{
register int r,c,s;
int halfsize;

halfsize = 1 << *level;

switch(*Iarrayix)   {

  case IX_HHH:

    Rprintf("Inserting HHH\n");

    ACCESSW3D(Carray, (int)*truesize, (int)*truesize, 0, 0, 0) =
        ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, 0, 0, 0);

    break;
    

  case IX_GHH:

    Rprintf("Inserting GHH\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c, s)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;

  case IX_HGH:

    Rprintf("Inserting HGH\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r, c+(int)halfsize, s)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;

  case IX_GGH:

    Rprintf("Inserting GGH\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c+(int)halfsize, s)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;

  case IX_HHG:

    Rprintf("Inserting HHG\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r, c, s+(int)halfsize)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;

  case IX_GHG:

    Rprintf("Inserting GHG\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c, s+(int)halfsize)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;


  case IX_HGG:

    Rprintf("Inserting HGG\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r, c+(int)halfsize, s+(int)halfsize)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;

  case IX_GGG:

    Rprintf("Inserting GGG\n");

    for(r=0; r < (int)halfsize; ++r)
      for(c=0; c < (int)halfsize; ++c)
        for(s=0; s < (int)halfsize; ++s)    {
          ACCESSW3D(Carray, (int)*truesize, (int)*truesize,
        r+(int)halfsize, c+(int)halfsize, s+(int)halfsize)=
            ACCESSW3D(Iarray, (int)halfsize, (int)halfsize, r, c, s); 
        }

    break;
  default:

    Rprintf("Unknown insertion type\n");
    break;
  }
}

/*
 * WaveThresh3 - Beginning of TRD's multiwavelet code.
 */

/*Multiple wavelet decomposition */
/*TRD November 1994              */
/*last updated May 1995          */

void multiwd(double *C, int *lengthc, double *D, int *lengthd, int *nlevels,
	int *nphi, int *npsi, int *ndecim, double *H, double *G, int *NH,
	int *lowerc, int *upperc, int *offsetc,
	int *lowerd, int *upperd, int *offsetd, int *nbc)
/*---------------------
 * Argument description
 *---------------------
double *C::       C coefficients matrix 
int *lengthc::    number of coefficients in C 
double *D::       D coefficients matrix 
int  *lengthd::   number of coefficients in D 
int  *nlevels::   number of levels in decomposition 
int  *nphi::      number of scaling functions 
int  *npsi::      number of wavelet functions
int  *ndecim::    amount of decimation at each level
double  *H::      Band pass filter
double  *G::      High pass filter
int  *NH::        number of coeff matrices in the filter 
int  *lowerc::    for each level the lowest C coefficient 
int  *upperc::    for each level the highest C coefficient 
int  *offsetc::   amount to offset to access each level  
int  *upperd::    for each level the lowest C coefficient 
int  *lowerd::    for each level the highest C coefficient 
int  *offsetd::   amount to offset to access each level  
int  *nbc::       boundary conds 1=period 2=symm. 
 *---------------------*/
{
 int level,prevlvl,prevoffsetc,index,base,k,l,m,n;
 void TRDerror(char *s);

int trd_reflect(int a, int b);

/* ...  see L10209,102010 */
int trd_module(int a, int b);

 for(level=*nlevels-1;level >=0;level--)
  {
   /*some frequently used values computed here */

   prevlvl=level+1;                           /* previous level */
   prevoffsetc=*(offsetc+prevlvl);            /*offset of C for previous level*/
   /* prevoffsetd=*(offsetd+prevlvl); NOT USED            offset of D for previous level*/

   for(k=*(lowerc+level);k<=*(upperc+level);k++){ /*k index of new vector */
    for(l=0; l<*nphi;l++){                         /*l index of new elmt   */ 
     C[(*(offsetc+level)+k-*(lowerc+level))*(*nphi)+l]=0.0;
     for(m=*ndecim*k;m<*ndecim*k+*NH;m++){  /*index of already known vectors */
                  /*using periodic boundary conditions index = m mod #coeffs */
      index = m-*(lowerc+prevlvl);
      base=1+*(upperc+prevlvl)-*(lowerc+prevlvl);
      if(index >= base || index < 0){
       if(*nbc==1) index = trd_module(index,base);
       else if(*nbc==2) index = trd_reflect(index,base);
       else TRDerror("bad boundary conditions\n");
      }
      for(n=0;n<*nphi;n++){                   /* index of already known elemnt */
       C[(*(offsetc+level)+k-*(lowerc+level))*(*nphi)+l]+=
              H[((m-*ndecim*k)*(*nphi)+l)*(*nphi)+n]*C[*nphi*(prevoffsetc+index)+n];
      }
    }
    }
   }
   for(k=*(lowerd+level);k<=*(upperd+level);k++){ /* repeat for D */
    for(l=0;l<*npsi;l++){ 
     D[(*(offsetd+level)+k-*(lowerd+level))**npsi+l]=0.0;
     for(m=*ndecim*k;m<*ndecim*k+*NH;m++){
      index = m-*(lowerc+prevlvl);
      base = 1+*(upperc+prevlvl)-*(lowerc+prevlvl);
      if(index >= base || index < 0){
       if(*nbc==1) index = trd_module(index,base);
       else if(*nbc==2) index = trd_reflect(index,base);
       else TRDerror("bad boundary conditions\n");
      }
      for(n=0;n<*nphi;n++){
       D[(*(offsetd+level)+k-*(lowerd+level))* *npsi+l]+=
              G[((m-*ndecim*k)* *npsi+l)**nphi+n]*C[(prevoffsetc+index)* *nphi+n];
      }
     }
    }
   }
 }
}
 
/*Double wavelet reconstruction */
/*By T Downie November 1994 */
/*updated Jan 95            */

void multiwr(double *C, int *lengthc, double *D, int *lengthd, int *nlevels,
	int *nphi, int *npsi, int *ndecim, double *H, double *G, int *NH,
	int *lowerc, int *upperc, int *offsetc,
	int *lowerd, int *upperd, int *offsetd, int *nbc, int *startlevel)
/*---------------------
 * Argument description
 *---------------------
double *C::       C coefficients matrix 
int *lengthc::    number of coefficients in C 
double *D::       D coefficients matrix 
int  *lengthd::   number of coefficients in D 
int  *nlevels::   number of levels in decomposition 
int  *nphi::      number of scaling functions 
int  *npsi::      number of wavelet functions
int  *ndecim::    amount of decimation at each level
double  *H::      Band pass filter
double  *G::      High pass filter
int  *NH::        number of coeff matrices in the filter 
int  *lowerc::    for each level the lowest C coefficient 
int  *upperc::    for each level the highest C coefficient 
int  *offsetc::   amount to offset to access each level  
int  *upperd::    for each level the lowest C coefficient 
int  *lowerd::    for each level the highest C coefficient 
int  *offsetd::   amount to offset to access each level  
int  *nbc::       boundary conds 1=period 2=symm. 
int  *startlevel; level at which to start the wavelet reconstruction
 *---------------------*/
{
 int level,offslvlc,offslvld,index,base,newck,newcl,oldck,oldcl,olddl,lim;

int trd_module(int a, int b);
int trd_reflect(int a, int b);

  for(level=*startlevel; level<*nlevels; level++){     /*level=level of convolution*/
    offslvlc=*(offsetc+level); /*ammount to offset C for this level */
    offslvld=*(offsetd+level); /*ammount to offste D for this level */
    for(newck=*(lowerc+level+1); newck<=*(upperc+level+1);newck++){
                               /*newck=position of the new c coeff*/
      for(newcl=0; newcl< *nphi;newcl++){
                               /*newcl=element of the new c coeff vector*/
        lim= newck+1-*NH;
        while(lim % *ndecim != 0) lim++;
        for(oldck=lim/ *ndecim; oldck<= ((float) newck) / *ndecim;oldck++){ 
         
                               /*oldck=position of the c/d coeff in conv. */
          for(oldcl=0;oldcl< *nphi;oldcl++){
                               /*oldcl=element of the c coeff in conv.*/
            index=oldck- *(lowerc+level);
            base= 1+*(upperc+level)-*(lowerc+level);
            if(index < 0 || index >= base){
              if(*nbc == 1) index=trd_module(index,base);
              else index=trd_reflect(index,base);
            }
            C[(*(offsetc+level+1)+newck)* *nphi + newcl] +=
                     H[((newck-*ndecim*oldck) * *nphi+oldcl) * *nphi + newcl]*
                     C[(offslvlc+index)* *nphi + oldcl];
          }
          for(olddl=0;olddl< *npsi;olddl++){
                               /*olddl=element of the d coeff in conv.*/
            index=oldck- *(lowerd+level);
            base= 1+*(upperd+level)-*(lowerd+level);
            if(index < 0 || index >= base){
              if(*nbc == 1) index=trd_module(index,base);
              else index=trd_reflect(index,base);
            }
            C[(*(offsetc+level+1)+newck)* *nphi + newcl] +=
                     G[((newck-*ndecim*oldck) * *nphi+olddl) * *npsi + newcl]*
                     D[(offslvld+index)* *npsi + olddl];
          }
        }
      }
    }
  }
}

int trd_reflect(int a, int b)
{
int trd_module(int a, int b);

  if(b <= 0) return (-1);
  else {
   if (a < -b || a > 2*b)  a=trd_module(a,2*b);
   if (a < 0) a=-1*a-1;
   if (a > b) a=2*b-a-1;
  }
  return(a);
}


int trd_module(int a, int b)
{
 /* robust modulus function */
 /* returns a (mod b) for b >0 and any integer a */
 /* returns -1 if b <= 0 */

 if (b <= 0) return(-1);
 else if(a > 0) while(a >= b)  a -= b;
 else if(a < 0) while(a < 0)  a +=b ;
 return(a);
}

/*
 * WaveThresh3 - End of TRD's multiwavelet code
 */

/*
 * 
 * IsPowerOfTwo(n)
 *
 * Returns log to the base 2 of n
 *
 * e.g. if n = 2^J then IsPowerOfTwo(n) is J
 *
 * If n is not a power of two or is not positive then -1 is returned.
 *
 * Author: GPN
 *
 */

int IsPowerOfTwo(int n)
{
int cnt = 0;

if (n<=0)
    return((-1));

while (!(0x01 & n)) {

    ++cnt;
    n >>= 1;
    }

if (n > 1)
    return((-1));
else
    return(cnt);
}

void TRDerror(char *s)
{
REprintf("Module TRDerror in WaveThresh\n");
REprintf("%s", s);
error("This should not happen. Stopping.\n");
}

/* Following functions are to do Complex-valued non-decimated
 * wavelet transform PACKET version (i.e. though wst/AvBasis
 */


#define POINTDR(l,i) (DataR + (*LengthData*(l)) + (i))
#define POINTDI(l,i) (DataI + (*LengthData*(l)) + (i))
#define POINTCR(l,i) (CaR + (*LengthData*(l)) + (i))
#define POINTCI(l,i) (CaI + (*LengthData*(l)) + (i))



/*
 * COMWST:  Complex-valued packet-ordered non-decimated transform
 */

void comwst(double *CaR, double *CaI, double *DataR, double *DataI,
	int *LengthData, int *levels,
	double *HR, double *HI, double *GR, double *GI,
	int *LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *CaR::     Will contain bottom most Cs (real)          
double *CaI::     Will contain bottom most Cs (imaginary)     
double *DataR::   This is a 2D array. Zeroeth level contains data 
double *DataI::   This is a 2D array. Zeroeth level contains data 
int *LengthData:: Length of Data, this is power of 2              
int *levels::     The number of levels, 2^(*levels)=LengthData    
double *HR::      Smoothing filter (real)             
double *HI::      Smoothing filter (imag)             
double *GR::      Detail filter (real)                
double *GI::      Detail filter (imag)    
int *LengthH::    Length of filter                
int *error::      Error code, if non-zero then it's a mem error   
 *---------------------*/
{
int startin, outstart1, outstart2;
register int i;
double *bookR, *bookI;  /* Bookkeeping vectors, one for R and I        */

void comwvpkstr(double *CaR, double *CaI, double *DataR, double *DataI,
	int startin, int lengthin, int outstart1, int outstart2, int level,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	int *LengthData, double *bookR, double *bookI, int *error);

*error = 0;

/*
Rprintf("This routine is wavepackst\n");
Rprintf("Length of data is %ld\n", *LengthData);
Rprintf("Number of levels is %ld\n", *levels);
Rprintf("Data array is:\n");
for(i= (int)*levels; i>=0; --i)
    for(j=0; j< *LengthData; ++j)   {
        Rprintf("Level %d, Item %d is %lf\n", i,j, ACCESSD(i,j));
        }
*/

/* Create a bookeeping vector. That contains the C,C' level smooths
   thoughout the algorithm. One for imag as well */

if ((bookR = (double *)malloc((unsigned)*LengthData*sizeof(double)))==NULL){
    *error = 1;
    return;
    }
    
if ((bookI = (double *)malloc((unsigned)*LengthData*sizeof(double)))==NULL){
    *error = 2;
    return;
    }
    
/* Copy original data to book keeping vector */

for(i=0; i< *LengthData; ++i)   {
    *(bookR+i) = *POINTDR(*levels, i);
    *(bookI+i) = *POINTDI(*levels, i);
    }


startin = 0;
outstart1 = 0;
outstart2 = ((int)*LengthData)/2;

comwvpkstr(CaR, CaI, DataR, DataI, startin, (int)*LengthData,
    outstart1, outstart2,
    (int)*levels,
    HR, HI, GR, GI, (int)*LengthH, LengthData, bookR, bookI, error);

if (*error != 0)
    return;
else    {
    free((void *)bookR);
    free((void *)bookI);
    }
}

void comwvpkstr(double *CaR, double *CaI, double *DataR, double *DataI,
	int startin, int lengthin, int outstart1, int outstart2, int level,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	int *LengthData, double *bookR, double *bookI, int *error)
{
register int i;
int lengthout;
double *book1R, *book1I, *book2R, *book2I;

void comconC(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *HR, double *HI, int LengthH,
	double *c_outR, double *c_outI,
	int LengthCout, int firstCout, int lastCout,
	int type, int step_factor, int bc);
void comconD(double *c_inR, double *c_inI,
	int LengthCin, int firstCin,
	double *GR, double *GI, int LengthH,
	double *d_outR, double *d_outI,
	int LengthDout, int firstDout, int lastDout,
	int type, int step_factor, int bc);
void comrotater(double *bookR, double *bookI, int length);
void comwvpkstr(double *CaR, double *CaI, double *DataR, double *DataI,
	int startin, int lengthin, int outstart1, int outstart2, int level,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	int *LengthData, double *bookR, double *bookI, int *error);

/*
Rprintf("wvpkstr entry\n");
Rprintf("lengthout is %d\n", lengthout);
*/

lengthout = lengthin/2;

if ((book1R = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL) {
    *error = 3;
    return;
    }

else if ((book1I = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL) {
    *error = 4;
    return;
    }

else if ((book2R = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL){
    *error = 5;
    return;
    }
else if ((book2I = (double *)malloc((unsigned)lengthout*sizeof(double)))==NULL){
    *error = 6;
    return;
    }

comconC(bookR, bookI, lengthin, 0, HR, HI, LengthH,
    book1R, book1I, lengthout, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);


for(i=0; i < lengthout; ++i)    {
    *POINTCR(level-1, (outstart1+i)) = *(book1R+i);
    *POINTCI(level-1, (outstart1+i)) = *(book1I+i);
    }

/*
Rprintf("book1 coefficients \n");
for(i=0; i<lengthout; ++i)
    Rprintf("%lf ", *(book1+i));;
Rprintf("\n");
*/

comconD(bookR, bookI, lengthin, 0, GR, GI, LengthH,
    POINTDR(level-1, outstart1),
    POINTDI(level-1, outstart1),
    lengthout, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

/* Now cycle book around one, and do the convolutions again */

/* COMMENT OUT
 * tmp = *(book+lengthin-1);
 *
 * for(i=lengthin-1; i>0; --i)
 *  *(book+i) = *(book+i-1);
 * book = tmp;
 */

/* COMMENT OUT (replaced by rotater function) tmp = *book;
 * for(i=0; i<lengthin-1; ++i)
 *  *(book+i) = *(book+i+1);
 *
 * *(book+lengthin-1) = tmp;
 */

comrotater(bookR, bookI, lengthin);

comconC(bookR, bookI, lengthin, 0, HR, HI, LengthH,
    book2R, book2I, lengthout, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);

for(i=0; i < lengthout; ++i)    {
    *POINTCR(level-1, (outstart2+i)) = *(book2R+i);
    *POINTCI(level-1, (outstart2+i)) = *(book2I+i);
    }

/*
Rprintf("book2 coefficients \n");
for(i=0; i<lengthout; ++i)
    Rprintf("%lf ", *(book2+i));
Rprintf("\n");
*/


comconD( bookR, bookI, lengthin, 0, GR, GI, LengthH,
    POINTDR(level-1, outstart2),
    POINTDI(level-1, outstart2),
    lengthout, 0, lengthout-1, 
    WAVELET, 1, PERIODIC);


if (lengthout!=1)   {
    /*
     * Now pass book1 and book2 to the next level 
     */

    comwvpkstr(CaR, CaI, DataR, DataI, outstart1, lengthout,
        outstart1, outstart1+lengthout/2,
        level-1,
        HR, HI, GR, GI, LengthH, LengthData,
        book1R, book1I, error); 

    if (*error != 0)
        return;

    comwvpkstr(CaR, CaI, DataR, DataI, outstart2, lengthout,
        outstart2, outstart2+lengthout/2,
        level-1,
        HR, HI, GR, GI, LengthH, LengthData,
        book2R, book2I, error); 

    if (*error != 0)
        return;
    }


free((void *)book1R);
free((void *)book1I);
free((void *)book2R);
free((void *)book2I);
}

/*
 * COMROTATER: complex version of rotater
 */

void comrotater(double *bookR, double *bookI, int length)
{
register int i;
double tmpR,tmpI;

tmpR = *bookR;
tmpI = *bookI;

for(i=0; i<length-1; ++i)   {
        *(bookR+i) = *(bookR+i+1);
        *(bookI+i) = *(bookI+i+1);
    }
        
*(bookR+length-1) = tmpR;
*(bookI+length-1) = tmpI;
}

struct complex  {
    double *realval;
    double *imagval;
    };



void comAB_WRAP(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int *LengthData, int *level,
    double *HR, double *HI, double *GR, double *GI, int *LengthH,
    double *answerR, double *answerI, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wstR::        Wavelet coefficients - real          
double *wstI::        Wavelet coefficients - imag          
double *wstCR::       Father coeffs - real               
double *wstCI::       Father coeffs - imag             
int *LengthData::
int *level::
double *HR, *HI::     Smoothing filter, real and imag      
double *GR, *GI::     Detail filter, real and imag         
int *LengthH::
double *answerR, *answerI::   Real and imag of answer      
int *error::
 *---------------------*/
{
register int i;
int nlevels;
struct complex *acopy;
struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error);
void destroycomplex(struct complex *a);

nlevels = 2 + (int)*level;

acopy =  comAB(wstR, wstI, wstCR, wstCI, nlevels, (int)*level, 0, 1,
        HR, HI, GR, GI, (int)*LengthH, error);

for(i=0; i< (int)*LengthData; ++i)  {
    *(answerR+i) = *(acopy->realval+i);
    *(answerI+i) = *(acopy->imagval+i);
    }

destroycomplex(acopy);
}

void destroycomplex(struct complex *a)
{
free((void *)a->realval);
free((void *)a->imagval);
free((void *)a);
}

/* comAB Do the basis averaging for complex WST*/

/*
 * Error codes
 *
 * 1,2  -   Memory error in creating clR, clI
 * 3,4  -   Memory error in creating crR, crI
 * 3    -   Memory error in creating packet (getpacket)
 */


struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error)
/*---------------------
 * Argument description
 *---------------------
double *wstR::    Wavelet coefficients, non-dec, real            
double *wstI::    Wavelet coefficients, non-dec, imag            
double *wstCR::   Father wav. coeffs, non-dec, real            
double *wstCI::   Father wav. coeffs, non-dec, imag            
int nlevels::     The original length of the data         
int level::       The level to reconstruct             
int ix1::         The "left" packet index              
int ix2::         The "right" packet index             
double *HR,*HI::  Smoothing filter                 
double *GR,*GI::  Detail filter                    
int LengthH::     The length of the filter             
int *error::      Error code                       
 *---------------------*/
{
register int i;
double *clR, *clI;
double *crR, *crI;
struct complex *genericC;
struct complex *answer;
double *genCR, *genCI;  /* Generic Cs for when we need real and imag */
double *genDR, *genDI;  /* Generic Cs for when we need real and imag */
int LengthC;
int LengthCin;

void comcbr(double *c_inR, double *c_inI,
	int LengthCin, int firstCin, int lastCin,
	double *d_inR, double *d_inI,
	int LengthDin, int firstDin, int lastDin,
	double *HR, double *HI, double *GR, double *GI, int LengthH,
	double *c_outR, double *c_outI, int LengthCout, int firstCout,
	int lastCout, int type, int bc);
double *getpacket(double *wst, int nlevels, int level, int index, int *error);
struct complex *comAB(double *wstR, double *wstI, double *wstCR, double *wstCI,
    int nlevels, int level, int ix1, int ix2,
    double *HR, double *HI, double *GR, double *GI, int LengthH, int *error);
void rotateback(double *book, int length);
void destroycomplex(struct complex *a);

*error = 0;

/*
 * Now we must create cl and cr. These will contain the reconstructions
 * from the left and right packets respectively. The length of these
 * vectors depends upon the level we're at.
 */

LengthC = 1 << (level+1);
LengthCin = 1 << level;

/*
 * Create cl and cr: real and imaginary
 */

if ((clR = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 1;
    return(NULL);
    }

if ((clI = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 2;
    return(NULL);
    }

if ((crR = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 3;
    return(NULL);
    }

if ((crI = (double *)malloc((unsigned)LengthC*sizeof(double)))==NULL) {
    *error = 4;
    return(NULL);
    }

/*
 * What we do next depends on the level.
 *
 * If level is zero then we've recursed all the way down to the bottom of
 * the tree. And we can reconstruct the 2-vectors one-up-the-tree by using
 * good old conbar().
 *
 * If the level is not zero then we construct at that stage using conbar()
 * but to obtain the Cs we recurse. 
 */

if (level != 0) {

    /* Get C's at this level by asking the next level down. */

    genericC = comAB(wstR, wstI, wstCR, wstCI,
            nlevels, level-1, 2*ix1, 2*ix1+1,
            HR, HI, GR, GI, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genDR = getpacket(wstR, nlevels, level, ix1, error);
    genDI = getpacket(wstI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genericC->realval, genericC->imagval, LengthCin, 0, LengthCin-1, 
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           clR, clI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    destroycomplex(genericC);
    free((void *)genDR);
    free((void *)genDI);

    /* Now do the RHS */
    
    genericC = comAB(wstR, wstI, wstCR, wstCI, nlevels, level-1,
        2*ix2, 2*ix2+1,
        HR, HI, GR, GI, LengthH, error);

    if (*error != 0)
        return(NULL); 

    /* Get D's straight from the wst matrix */

    genDR = getpacket(wstR, nlevels, level, ix2, error);
    genDI = getpacket(wstI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genericC->realval, genericC->imagval, LengthCin, 0, LengthCin-1,
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           crR, crI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(crR, LengthC);
    rotateback(crI, LengthC);

    /* Can get rid of generics now */

    destroycomplex(genericC);
    free((void *)genDR);
    free((void *)genDI);
    }

else    {
    /* Have to really do it! */

    genCR = getpacket(wstCR, nlevels, level, ix1, error);
    genCI = getpacket(wstCI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    genDR = getpacket(wstR, nlevels, level, ix1, error);
    genDI = getpacket(wstI, nlevels, level, ix1, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genCR, genCI, LengthCin, 0, LengthCin-1, 
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           clR, clI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    free((void *)genCR);
    free((void *)genCI);
    free((void *)genDR);
    free((void *)genDI);

    genCR = getpacket(wstCR, nlevels, level, ix2, error);
    genCI = getpacket(wstCI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    genDR = getpacket(wstR, nlevels, level, ix2, error);
    genDI = getpacket(wstI, nlevels, level, ix2, error);

    if (*error != 0)
        return(NULL);

    /* Do the reconstruction */

    comcbr(genCR, genCI, LengthCin, 0, LengthCin-1,
           genDR, genDI, LengthCin, 0, LengthCin-1,
           HR, HI, GR, GI, LengthH,
           crR, crI, LengthC, 0, LengthC-1,
           WAVELET, PERIODIC);

    /* Rotate the RHS back */

    rotateback(crR, LengthC);
    rotateback(crI, LengthC);

    free((void *)genCR);
    free((void *)genCI);
    free((void *)genDR);
    free((void *)genDI);
    }

for(i=0; i<LengthC; ++i)    {
    *(clR+i) = ((double)0.5)*( *(clR+i) + *(crR+i) );
    *(clI+i) = ((double)0.5)*( *(clI+i) + *(crI+i) );
    }

if ((answer=(struct complex *)malloc((unsigned)sizeof(struct complex)))==NULL) {
    *error = 5l;
    return(NULL);
    }

answer->realval = clR;
answer->imagval = clI;

return(answer);
}
