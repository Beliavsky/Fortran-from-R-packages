

#define STRICT_R_HEADERS

/* This file contains the additional code needed to perform wavelet
 * density estimation in SPlus.  The WaveThresh package must be installed.
 * The locations of the functions are:
 * 
 * SFDE5                line 0053
 * SFDE6                line 0096
 * PLDE2                line 0148
 * phi                  line 0194
 * diad                 line 0246
 * T                    line 0263
 * StoDCDS              line 0281
 * DensityCovarianceDecomposeStep   line 0334
 * AXSDCV               line 0444
 * StoIDS_dh                line 0470
 * ImageDecomposeStep_dh        line 0535
 * wavedecomp_dh            line 0756
 * convolveC_dh             line 0874
 * convolveD_dh             line 0945
 * reflect_dh               line 1035
 * access0              line 1122
 * waverecons_dh            line 1137
 * conbar_dh                line 1238
 */

#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Error condition      */
#define OK              (0)

/* For boundary condition handling */
#define PERIODIC        1
#define SYMMETRIC       2
#define ZERO        3

/* For the type of wavelet decomposition */
#define WAVELET         1       /* The standard decomposition */
#define STATION         2       /* The stationary decomposition */

#define ACCESSC_DH(c, firstC, lengthC, ix, bc) *(c+reflect_dh(((ix)-(firstC)),(lengthC),(bc)))
#define AXSMAT(a, nrow, i, j) (a + (i) + (nrow)*(j))
#define ACCESS(image, size, i, j)       *(image + (i)*(size) + (j))

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

/* SFDE5 calculates empirical scaling function coefficients from data,
 * using the Daubechies-Lagarias algorithm */

void SFDE5(double *x, int *nx, double *p, double *filter, int *nf, int *prec,
	double *chat, int *kmin, int *kmax, double *philh, double *phirh,
	int *error)
/* double *x;   	The data                     			*/
/* int *nx;    		Number of data points                		*/
/* double *p;   	The primary resolution               		*/
/* double *filter;	Vector of filter coefficients            	*/
/* int *nf;     	Number of filter coefficients - 1		*/
/* int *prec;   	Precision used in evaluating phi		*/
/* double *chat; 	Vector to put coefficient estimates in		*/
/* int *kmin;  		minimum value of k                   		*/
/* int *kmax;  		maximum value of k                    		*/
/* double *philh;  	Left hand end of scaling function support    	*/
/* double *phirh;	Right hand end of scaling function support   	*/
/* int *error;		Error code - mostly out of memory		*/
{
    void phi(double y, double *filt, double *out, int *pre, int *n, int *error);
    register int i, j, k;
    register int min, max;
    register double z;
    double *phix;

    phix = (double *) calloc(*nf+1, sizeof(double));

    if (phix == NULL)	{
	    *error = 1;
	    return;
    }

/* calculate coefficient estimates */

    for(i=0; i < *nx; i++) {

        for(j=0; j < *nf; j++) {
            *(phix+j) = 0.0;
        }

        z = *p * *(x+i);
        min = ceil(z-*phirh);
        max = floor(z-*philh);
        phi(z, filter, phix, prec, nf, error);
	if (*error != 0)
		return;
        for (k=min; k <= max; k++)
            *(chat+(k-*kmin)) += sqrt(*p) * *(phix + k - min) / *nx;
    }
    free((void *)phix);
}


/* As SFDE5, but also calculates covariances of the coefficients */

void SFDE6(double *x, int *nx, double *p, double *filter, int *nf, int *prec,
	double *chat, double *covar, int *kmin, int *kmax, double *philh,
	double *phirh, int *error)
/* double *x;   	The data                     			*/
/* int *nx;     	Number of data points                		*/
/* double *p;   	The primary resolution               		*/
/* double *filter;  	Vector of filter coefficients            	*/
/* int *nf;     	Number of filter coefficients - 1            	*/
/* int *prec;   	Precision used in evaluating phi   		*/
/* double *chat;    	Vector to put coefficient estimates in          */
/* double *covar;   	Matrix to put covariance estimates in        	*/
/* int *kmin;   	minimum value of k                   		*/
/* int *kmax;      	maximum value of k                    		*/
/* double *philh;  	Left hand end of scaling function support       */
/* double *phirh;	Right hand end of scaling function support	*/
/* int *error;		Error code -- mostly out of memory		*/
{
    void phi(double y, double *filt, double *out, int *pre, int *n, int *error);
    register int i, j, k, l;
    register int min, max;
    register double z, phijk, phijl;
    double *phix;

    *error = 0;

    phix = (double *) calloc(*nf+1, sizeof(double));

    if (phix == NULL)	{
	    *error = 1;
	    return;
    }


/* calculate coefficient estimates */

    for(i=0; i < *nx; i++) {

        for(j=0; j < *nf; j++) {
            *(phix+j) = 0.0;
        }

        z = *p * *(x+i);
        min = ceil(z-*phirh);
        max = floor(z-*philh);
        phi(z, filter, phix, prec, nf, error);
	if (*error != 0)
		return;
        for (k=min; k <= max; k++) {
            phijk = sqrt(*p) * *(phix + k - min);
            *(chat+(k-*kmin)) += phijk / *nx;
            for (l=k; (l < (k+*nf)) && (l <= max); l++) {
                phijl = sqrt(*p) * *(phix + l - min);
                *AXSMAT(covar, (*kmax-*kmin+1), (k-*kmin), (l-k)) +=
                                                phijk * phijl / (*nx * *nx);
            }
        }
    }
    free((void *)phix);
}


/* Function to get plotting information for density estimate from high
 * level scaling function coefficients */

void PLDE2(double *C, double *p, double *filter, int *nf, int *prec,
	int *kmin, int *kmax, double *gx, double *gy, int *ng,
	double *philh, double *phirh, int *error)
/* double *C;      High resolution scaling function coefficients        */
/* double *p;      The primary resolution                               */
/* double *filter; Vector of filter coefficients                        */
/* int *nf;        Number of filter coefficients - 1                    */
/* int *prec;      Precision used in evaluating phi                     */
/* int *kmin;      minimum value of k                                   */
/* int *kmax;      maximum value of k                                   */
/* double *gx;     grid for drawing density estimate                    */
/* double *gy;     Vector to put density values in                      */
/* int *ng;        Length of above grids                                */
/* double *philh;  Left hand end of scaling function support            */
/* double *phirh;  Right hand end of scaling function support           */
/* int *error;     Error Code						*/
{
    void phi(double y, double *filt, double *out, int *pre, int *n, int *error);
    register int i, j, k;
    register int min, max;
    register double z;
    double *phix;

    *error = 0;

    phix = (double *) calloc(*nf+1, sizeof(double));

    if (phix == NULL)	{
	    *error = 1;
	    return;
	    }

/* Evaluate density estimate over the grid provided */

    for (i=0; i < *ng; i++) {

        for(j=0; j < *nf; j++) {
            *(phix+j) = 0.0;
        }

        z = *p * *(gx+i);
        min = ceil(z-*phirh);
        max = floor(z-*philh);
        if(min<*kmin)
            min = *kmin;
        phi(z, filter, phix, prec, nf, error);

	if (*error != 0)
		return;

	{
		double a, b;

        for (k=min; k<=max && k<=*kmax; k++)	{
	
	    a = *(C+(k-*kmin));

	    b = *(phix + k - min);

            *(gy+i) += a * sqrt(*p) * b;
	    }

	}
    }
    free((void *)phix);

    return;
}


/* Function to evaluate phi_Jk(x) for all k for which it is non-zero */

void phi(double y, double *filt, double *out, int *pre, int *n, int *error)
{
    double T(int index, double *filter, int *n, int j, int k);
    void diad(double x, int *prec, int *out);
    int i, j, k, l, *dix;
    double z, *ans, *tmp;

    dix = (int *) calloc(*pre, sizeof(int));

    if (dix == NULL)	{
	    *error = 2;
	    return;
    }

    ans = (double *) calloc((*n * *n), sizeof(double));

    if (ans == NULL)	{
	    free((void *)dix);
	    *error = 3;
	    return;
	    }

    tmp = (double *) calloc((*n * *n), sizeof(double));

    if (tmp == NULL)	{
	    free((void *)dix);
	    free((void *)ans);
	    *error = 4;
	    return;
	    }

    for(i=0; i < *n; i++) {
        for(j = 0; j < *n; j++) {
            if(i==j) *AXSMAT(ans, *n, i, j) = 1.0;
            else *AXSMAT(ans, *n, i, j) = 0.0;
        }
    }

    z = y - floor(y);
    diad(z, pre, dix);

    for(i=0; i<*pre; i++) {

        for(j=0; j < *n; j++) {
            for(k=0; k < *n; k++) {

                *AXSMAT(tmp, *n, j, k) = 0.0;
                for(l=0; l < *n; l++) {
                    *AXSMAT(tmp, *n, j, k) += *AXSMAT(ans, *n, j, l) * 
                                            T(*(dix+i), filt, n, (l+1), (k+1));
                }
            }
        }

        for(j=0; j < *n; j++) {
            for(k=0; k < *n; k++) {
                *AXSMAT(ans, *n, j, k) = *AXSMAT(tmp, *n, j, k);
            }
        }
    }

    for(i=0; i< *n; i++) {
        for(j=0; j < *n; j++) {
            *(out + *n - 1 - i) += *AXSMAT(ans, *n, i, j) / *n;
        }
    }
    free((void *)dix);
    free((void *)ans);
    free((void *)tmp);
}


/* Function to find diadic representation of a number in (0,1) */

void diad(double x, int *prec, int *out)

{
    double nu;
    int i;

    nu = x;
    for(i=0; i<*prec; i++) {
        nu = 2*nu;
        *(out + i) = floor(nu);
        nu = nu - floor(nu);
    }
}


/* Function to find T_i(j,k) from filter */

double T(int index, double *filter, int *n, int jj, int kk)
{

int ind;

ind=-1;     /* MAN: added as a initialization, but *should* be set to true value
            assuming that index can only be either zero or one */
    if(index==0) 
        ind = 2*jj-kk-1;
    else if(index==1)
        ind = 2*jj-kk;

    if(ind < 0 || ind > *n)
        return(0.0);
    else return sqrt(2) * *(filter + ind);
}


/* Function for calling DensityCovarianceDecomposeStep from Splus */

void StoDCDS(double *C, int *Csize, int *firstCin, double *H, int *LengthH,
	int *LengthCout, int *firstCout, int *lastCout,
        int *LengthDout, int *firstDout, int *lastDout,
	double *ImCC, double *ImDD,
	int *bc, int *type, int *error)
{
    register int i,j;
    double *cc_out, *dd_out;

    void DensityCovarianceDecomposeStep(double *C, int Crow, int firstCin,	
	double *H, int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
        double **cc_out, double **dd_out,
	int bc, int type, int *error);

    double AXSDCV(double *a, int nr, int nc, int i, int j);

    DensityCovarianceDecomposeStep(C, *Csize, *firstCin, H, *LengthH,
        *LengthCout, *firstCout, *lastCout,
        *LengthDout, *firstDout, *lastDout,
        &cc_out, &dd_out, *bc, *type,
        error);

/* Copy images */

    for(i=0; i<(int)*LengthDout; ++i) {
        for(j=0; j < (*LengthH-1); ++j)
            *AXSMAT(ImDD, (int)*LengthDout, i, j) =
                                        *AXSMAT(dd_out, *LengthDout, i, j);
    }

    for(i=0; i<(int)*LengthCout; ++i) {
        for(j=0; j < (*LengthH-1); ++j)
            *AXSMAT(ImCC, (int)*LengthCout, i, j) =
                                        *AXSMAT(cc_out, *LengthCout, i, j);
    }

    R_Free(cc_out);
    R_Free(dd_out);
}


/* Function for decomposing the covariance matrix of scaling function coefs */

void DensityCovarianceDecomposeStep(double *C, int Crow, int firstCin,	
	double *H, int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
        double **cc_out, double **dd_out,
	int bc, int type, int *error)

/* double *C;       Input data image                                    */
/* int Crow;        Number of rows of covariance matrix                 */
/* int firstCin;    Index number of first element in input "C" image    */
/* double *H;       Filter coefficients                                 */
/* int LengthH;     Length of filter                                    */
/* Details about output image						*/
/* int LengthCout;  Length of C part of output image                    */
/* int firstCout;   Index number of first element in output "C" image   */
/* int lastCout;    Index number of last element                        */
/* int LengthDout;  Length of D part of output image                    */
/* int firstDout;   Index number of first element in output "D" image   */
/* int lastDout;    Index number of last element                        */
/* double **cc_out; Smoothed output image                               */
/* double **dd_out; Diagonal detail                                     */
/* int bc;          Method of boundary correction                       */
/* int type;        Type of transform, wavelet or stationary            */
/* int *error;      Error code                                          */
{
    register int k,l,kmin,kmax,lmin,lmax,row,col;
    double *afterCC,*afterDD;     /* Results              */
    /* int step_factor; NOT USED This should always be 1 for the WAVELET trans*/

    double AXSDCV(double *a, int nr, int nc, int i, int j);

    *error = 0l;

    if ((afterCC =
        (double *)R_Calloc( LengthCout*(LengthH-1), double))==NULL)
        {
    *error = 6l;
        return;
    }

    for (row=0; row < LengthCout; row++) {
        for(col=0; col < (LengthH-1); col++) {
            *AXSMAT(afterCC, LengthCout, row, col) = 0.0;
        }
    }

    if ((afterDD =
    (double *)R_Calloc( LengthDout*(LengthH-1),double))==NULL){
        *error = 9l;
        return;
    }

    for (row=0; row < LengthCout; row++) {
        for(col=0; col < (LengthH-1); col++) {
            *AXSMAT(afterDD, LengthDout, row, col) = 0.0;
        }
    }

/* Link this memory to the returning pointers */

    *cc_out = afterCC;
    *dd_out = afterDD;


    for(row=firstCin; row < (firstCin+Crow); row++) {
        for(col=max(row-LengthH+2, firstCin);
                col < min(row+LengthH-1, firstCin+Crow); col++) {
            kmin = (int)ceil((0.5 * (double)(row+1-LengthH)));
            kmax = (int)floor((0.5 * (double)(row)));
            lmin = (int)ceil((0.5 * (double)(col+1-LengthH)));
            lmax = (int)floor((0.5 * (double)(col)));

            for(k=kmin; k <= kmax; k++) {
                for(l=max(lmin, k); l<=min(lmax, k+LengthH-1); l++) {
                    *AXSMAT(afterCC, LengthCout, (k-firstCout), (l-k)) += 
                        *(H+row-2*k) * *(H+col-2*l)  *
                        AXSDCV(C, Crow, LengthH-1, row-firstCin, col-firstCin);

                }
            }
        }
    }


    for(row=firstCin; row < (firstCin+Crow); row++) {
        for(col=max(row-LengthH+2, firstCin);
                col < min(row+LengthH-1, firstCin+Crow); col++) {
            kmin = (int)ceil((0.5 * (double)(row-1)));
            kmax = (int)floor((0.5 * (double)(row-2+LengthH)));
            lmin = (int)ceil((0.5 * (double)(col-1)));
            lmax = (int)floor((0.5 * (double)(col-2+LengthH)));

            for(k=kmin; k <= kmax; k++) {
                for(l=max(lmin, k); l<=min(lmax, k+LengthH-1); l++) {
                    *AXSMAT(afterDD, LengthDout, (k-firstDout), (l-k)) +=
                        (int)pow(-1, row+col) * *(H+2*k+1-row) * *(H+2*l+1-col)  *
                        AXSDCV(C, Crow, LengthH-1, row-firstCin, col-firstCin);

                }
            }
        }
    }

/* That should be it ! */

    return;
}


/* Function for accessing elements of the covariance matrix */

double AXSDCV(double *a, int nr, int nc, int i, int j)

/* double *a;      Pointer to covariance object                         */
/* int nr;         Number of rows of a                                  */
/* int nc;         Number of columns of a                               */
/* int i;          First index                                          */
/* int j;          Second index                                         */

{
    int ti;

    if(i > j) {
        ti = i;
        i = j;
        j = ti;
    }

    if((j-i) < nc)
        return *(a + i + nr*(j-i));
    else 
        return 0.0;
}


/* As WaveThresh StoIDS, but allows for zero boundary conditions */

void StoIDS_dh(double *C, int *Csize, int *firstCin, double *H, int *LengthH,
        int *LengthCout, int *firstCout, int *lastCout,
        int *LengthDout, int *firstDout, int *lastDout,
        double *ImCC, double *ImCD, double *ImDC, double *ImDD,
	int *bc, int *type, int *error)
{
    register int i,j;
    double *cc_out, *cd_out, *dc_out, *dd_out;

    void ImageDecomposeStep_dh(double *C, int Csize, int firstCin, double *H,
	int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
	double **cc_out, double **cd_out, double **dc_out, double **dd_out,
	int bc, int type, int *error);

    ImageDecomposeStep_dh(C, *Csize, *firstCin, H, *LengthH,
        *LengthCout, *firstCout, *lastCout,
        *LengthDout, *firstDout, *lastDout,
        &cc_out, &cd_out, &dc_out, &dd_out, *bc, *type,
        error);

/* Copy images */

    for(i=0; i<(int)*LengthDout; ++i) {
        for(j=0; j<(int)*LengthDout; ++j)
            *AXSMAT(ImDD, (int)*LengthDout, i, j) =
                                        *AXSMAT(dd_out, *LengthDout, i, j);

        for(j=0; j<(int)*LengthCout; ++j)
            *AXSMAT(ImDC, (int)*LengthDout, j, i) =
                                        *AXSMAT(dc_out, *LengthDout, j, i);
    }

    for(i=0; i<(int)*LengthCout; ++i) {
        for(j=0; j<(int)*LengthDout; ++j)
            *AXSMAT(ImCD, (int)*LengthCout, j, i) =
                                        *AXSMAT(cd_out, *LengthCout, j, i);

        for(j=0; j<(int)*LengthCout; ++j)
            *AXSMAT(ImCC, (int)*LengthCout, j, i) =
                                        *AXSMAT(cc_out, *LengthCout, j, i);
    }

    R_Free(cc_out);
    R_Free(cd_out);
    R_Free(dc_out);
    R_Free(dd_out);
}


/* As WaveThresh ImageDecomposeStep, but allows for zero boundary conditions */

void ImageDecomposeStep_dh(double *C, int Csize, int firstCin, double *H,
	int LengthH,
	int LengthCout, int firstCout, int lastCout,
	int LengthDout, int firstDout, int lastDout,
	double **cc_out, double **cd_out, double **dc_out, double **dd_out,
	int bc, int type, int *error)
/* double *C;       Input data image                                    */
/* int Csize;       Size of image (side length)                         */
/* int firstCin;    Index number of first element in input "C" image    */
/* double *H;       Filter coefficients                                 */
/* int LengthH;     Length of filter                                    */
/* Details about output image						*/
/* int LengthCout;  Length of C part of output image                    */
/* int firstCout;   Index number of first element in output "C" image   */
/* int lastCout;    Index number of last element                        */
/* int LengthDout;  Length of D part of output image                    */
/* int firstDout;   Index number of first element in output "D" image   */
/* int lastDout;    Index number of last element                        */
/* double **cc_out; Smoothed output image                               */
/* double **cd_out; Horizontal detail                                   */
/* double **dc_out; Vertical detail                                     */
/* double **dd_out; Diagonal detail                                     */
/* int bc;          Method of boundary correction                       */
/* int type;        Type of transform, wavelet or stationary            */
/* int *error;      Error code                                          */
{
    register int j,row,col;
    double *ccopy;  /* Used to copy input data to convolution routines      */
    double *ccopy_out;/* Used to copy output data to afterC after conv.     */
    double *dcopy_out;/* Used to copy output data to afterD after conv.     */
    double *afterC; /* Temporary store for image data after C convolution   */
    double *afterD; /* Temporary store for image data after D convolution   */
    double *afterCC,*afterCD,*afterDC,*afterDD;     /* Results              */
    int step_factor;        /* This should always be 1 for the WAVELET trans*/

    void convolveC_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *c_out,
	int firstCout, int lastCout, int type, int step_factor, int bc);

    void convolveD_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

    *error = 0l;

    step_factor = 1;

/* Get memory for afterC */

    if ((afterC =
    (double *)R_Calloc(Csize*LengthCout,double))==NULL){
        *error = 1l;
        return;
    }

/* Get memory for afterD */

    if ((afterD =
    (double *)R_Calloc(Csize*LengthDout,double))==NULL){
        *error = 2l;
        return;
    }

/* Get memory for row of image to pass to convolution routines */

    if ((ccopy =
    (double *)R_Calloc(Csize,double)) == NULL) {
        *error = 3l;
        return;
    }

/* Get memory for output row after C convolution */

    if ((ccopy_out =
    (double *)R_Calloc(LengthCout,double))==NULL) {
        *error = 4l;
        return;
    }

/* Get memory for output row after D convolution */

    if ((dcopy_out = (double *)R_Calloc(LengthDout,double))==NULL) {
        *error = 5l;
        return;
    }

/* Do convolutions on rows of C */

    for(row=0; row < (int)Csize; ++row)     {

/* Copy row of C into ccopy */

        for(j=0; j<Csize; ++j)
            *(ccopy+j) = *AXSMAT(C, Csize, row, j);

/* Now convolve this row with C filter */

        convolveC_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                ccopy_out, (int)firstCout, (int)lastCout, (int)type,
                step_factor, (int)bc);

/* Now convolve this row with D filter */

        convolveD_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                dcopy_out, (int)firstDout, (int)lastDout, (int)type,
                step_factor, (int)bc);

/* Copy answer back to arrays */

        for(j=0; j<(int)LengthCout; ++j)
            *AXSMAT(afterC, (int)Csize, row, j) = *(ccopy_out + j);

        for(j=0; j<(int)LengthDout; ++j)
            *AXSMAT(afterD, (int)Csize, row, j) = *(dcopy_out + j);

    }

/* Now we have to apply both the C and D filters to afterC and afterD.
 * We get four answers. First we get the necessary memory
 */

    if ((afterCC = (double *)R_Calloc(LengthCout*LengthCout,double))==NULL){
        *error = 6l;
        return;
    }

    if ((afterCD = (double *)R_Calloc(LengthDout*LengthCout,double))==NULL){
        *error = 7l;
        return;
    }

    if ((afterDC = (double *)R_Calloc(LengthCout*LengthDout,double))==NULL){
        *error = 8l;
        return;
    }

    if ((afterDD = (double *)R_Calloc(LengthDout*LengthDout,double))==NULL){
        *error = 9l;
        return;
    }

/* Link this memory to the returning pointers */

    *cc_out = afterCC;
    *cd_out = afterCD;
    *dc_out = afterDC;
    *dd_out = afterDD;


/* Apply the filters, first to afterC to get afterCC and afterCD */

    for(col=0; col < (int)LengthCout; ++col) {

/* Copy column to ccopy */

        for(j=0; j<(int)Csize; ++j)
            *(ccopy + j) = *AXSMAT(afterC, (int)Csize, j, col);

/* Apply C filter */

        convolveC_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                ccopy_out, (int)firstCout, (int)lastCout, (int)type,
                step_factor, (int)bc);

/* Apply D filter */

        convolveD_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                dcopy_out, (int)firstDout, (int)lastDout, (int)type,
                step_factor, (int)bc);

/* Copy answer back */

        for(j=0; j<(int)LengthCout; ++j)
            *AXSMAT(afterCC, (int)LengthCout, j, col) = *(ccopy_out+j);

        for(j=0; j<(int)LengthDout; ++j)
            *AXSMAT(afterCD, (int)LengthCout, j, col) = *(dcopy_out+j);
    }

/* Apply the filters, now to afterD to get afterDC and afterDD */

    for(col=0; col < (int)LengthDout; ++col)        {

/* Copy column to ccopy */

        for(j=0; j<(int)Csize; ++j)
            *(ccopy + j) = *AXSMAT(afterD, (int)Csize, j, col);

/* Apply C filter */

        convolveC_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                ccopy_out, (int)firstCout, (int)lastCout, (int)type,
                step_factor, (int)bc);

/* Apply D filter */

        convolveD_dh(ccopy, (int)Csize, (int)firstCin, H, (int)LengthH,
                dcopy_out, (int)firstDout, (int)lastDout, (int)type,
                step_factor, (int)bc);

/* Copy answer back */

        for(j=0; j<(int)LengthCout; ++j)
            *AXSMAT(afterDC, (int)LengthDout, j, col) = *(ccopy_out+j);

        for(j=0; j<(int)LengthDout; ++j)
            *AXSMAT(afterDD, (int)LengthDout, j, col) = *(dcopy_out+j);
    }

/* That should be it ! */

    R_Free(afterD);
    R_Free(afterC);
    R_Free(dcopy_out);
    R_Free(ccopy_out);
    R_Free(ccopy);

    return;
}


/* As WaveThresh wavedecomp, but allows for zero boundary conditions */

void wavedecomp_dh(double *C, double *D, double *H, int *LengthH, int *levels,
	int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/* double *C;              Input data, and the subsequent smoothed data */
/* double *D;              The wavelet coefficients                     */
/* double *H;              The smoothing filter H                       */
/* int *LengthH;           Length of smoothing filter                   */
/* int *levels;            The number of levels in this decomposition   */
/* int *firstC;            The first possible C coef at a given level   */
/* int *lastC;             The last possible C coef at a given level    */
/* int *offsetC;           Offset from C[0] for certain level's coeffs  */
/* int *firstD;            The first possible D coef at a given level   */
/* int *lastD;             The last possible D coef at a given level    */
/* int *offsetD;           Offset from D[0] for certain level's coeffs  */
/* int *type;              The type of wavelet decomposition            */
/* int *bc;                Method of boundary correction                */
/* int *error;             Error code                                   */
{
register int next_level,at_level;
register int step_factor;       /* Controls width of filter for station */
register int verbose;   /* Controls message printing, passed in error var*/

void convolveC_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *c_out,
	int firstCout, int lastCout, int type, int step_factor, int bc);

void convolveD_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc);

if (*error == 1l)       /* Error switches on verbosity */
        verbose = 1;
else
        verbose = 0;

switch(*bc)     {

        case PERIODIC:  /* Periodic boundary conditions */
                if (verbose) Rprintf("Periodic boundary method\n");
                break;

        case SYMMETRIC: /* Symmetric boundary conditions */
                if (verbose) Rprintf("Symmetric boundary method\n");
                break;

        case ZERO:  /* Zero boundary conditions */
        if (verbose) Rprintf("Zero boundary method\n");
        break;

        default:        /* The bc must be one of the above */
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

        default:        /* The type must be of one the above */
                if (verbose) Rprintf("Unknown decomposition type\n");
                *error = 2;
                return;
        }

if (verbose) Rprintf("Decomposing into level: ");

*error = 0l;

step_factor = 1;        /* This variable should *always* be 1 for standard
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

        convolveC_dh( (C+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                H,
                (int)*LengthH,
                (C+*(offsetC+next_level)),
                (int)(*(firstC+next_level)),
                (int)(*(lastC+next_level)) , (int)*type,
                step_factor, (int)*bc);

        convolveD_dh( (C+*(offsetC+at_level)),
                (int)(*(lastC+ at_level) - *(firstC+at_level)+1),
                (int)(*(firstC+at_level)),
                H,
                (int)*LengthH,
                (D+*(offsetD+next_level)),
                (int)(*(firstD+next_level)),
                (int)(*(lastD+next_level)), (int)*type,
                step_factor, (int)*bc );

        if (*type == STATION)
                step_factor *= 2;       /* Any half decent compiler should
                                         * know what to do here ! */
        }
if (verbose)
        Rprintf("\n");
return;
}


/* As WaveThresh convolveC, but allows for zero boundary conditions */

void convolveC_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *c_out,
	int firstCout, int lastCout, int type, int step_factor, int bc)
/* double *c_in;    Input data                                           */
/* int LengthCin;   Length of this array                                 */
/* int firstCin;    The first C value                                    */
/* double *H;       Filter                                               */
/* int LengthH;     Length of filter                                     */
/* double *c_out;   Output data                                          */
/* int firstCout;   First index of C array                               */
/* int lastCout;    Last index of C array                                */
/* int type;        Type of wavelet decomposition                        */
/* int step_factor; For stationary wavelets only                         */
/* int bc;          Method of boundary correction PERIODIC, SYMMETRIC    */
{
double sum;
register int k;
register int count_out;
register int m;
register int cfactor;   /* This determines what sort of dilation we do  */
                        /* and depends on the type argument             */
int reflect_dh(int n, int lengthC, int bc);
double access0(double *c, int lengthC, int n);

count_out = 0;

switch(type)    {

        case WAVELET:   /*  Ordinary wavelets   */
                        cfactor = 2;    /* Pick every other coefficient */
                        break;

        case STATION:   /* Stationary wavelets  */
                        cfactor = 1;    /* Pick every coefficient       */
                        break;


        default:        /* This is an error, one of the above must have */
                        /* been picked */
                        /* However, this must be tested in a previous   */
                        /* routine.                                     */
                        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
                        break;
                }
if (bc==ZERO) {
    for(k=firstCout; k<=lastCout; ++k)      {
        sum = 0.0;
        for(m=0; m<LengthH; ++m)        {
            sum += *(H+m) * access0(c_in, LengthCin,
                ((step_factor*m)+(cfactor*k)-firstCin));
        }
        *(c_out + count_out) = sum;
        ++count_out;
    }
}

else {
    for(k=firstCout; k<=lastCout; ++k)      {
        sum = 0.0;
        for(m=0; m<LengthH; ++m)        {
            sum += *(H+m) * ACCESSC_DH(c_in, firstCin, LengthCin,
                ((step_factor*m)+(cfactor*k)),bc);
        }

        *(c_out + count_out) = sum;
        ++count_out;
    }
}
}


/* As WaveThresh convolveD, but allows for zero boundary conditions */

void convolveD_dh(double *c_in, int LengthCin, int firstCin, double *H,
	int LengthH, double *d_out, int firstDout, int lastDout,
	int type, int step_factor, int bc)
/* double *c_in;    Input data                                           */
/* int LengthCin;   Length of this array                                 */
/* int firstCin;    The first index of the C input array                 */
/* double *H;       Filter                                               */
/* int LengthH;     Length of filter                                     */
/* double *d_out;   Output data                                          */
/* int firstDout;   First index of D array                               */
/* int lastDout;    Last index of D array                                */
/* int type;        Type of wavelet decomposition                        */
/* int step_factor; For stationary wavelets only                         */
/* int bc;          Method of boundary correction PERIODIC or SYMMETRIC  */
{
double sum;
double tmp;
register int k;
register int count_out;
register int m;
register int cfactor;

int reflect_dh(int n, int lengthC, int bc);
double access0(double *c, int lengthC, int n);

count_out = 0;

switch(type)    {

        case WAVELET:   /*  Ordinary wavelets   */
                        cfactor = 2;    /* Pick every other coefficient */
                        break;

        case STATION:   /* Stationary wavelets  */
                        cfactor = 1;    /* Pick every coefficient       */
                        break;


        default:        /* This is an error, one of the above must have */
                        /* been picked */
                        /* However, this must be tested in a previous   */
                        /* routine.                                     */
                        cfactor=0;       /* MAN: added for total cover: shouldn't happen */
                        break;
                }

if (bc==ZERO) {
    for(k=firstDout; k<=lastDout; ++k)      {
        sum = 0.0;

        for(m=0; m<LengthH; ++m)        {

            tmp = access0(c_in, LengthCin, 
                (cfactor*k+(step_factor*(1-m))-firstCin));

            if (m&1)        /* odd */
                sum += *(H+m) *  tmp;
            else
                sum -= *(H+m) *  tmp;

            }

        *(d_out + count_out) = sum;
        ++count_out;
    }
}

else {
    for(k=firstDout; k<=lastDout; ++k)      {
        sum = 0.0;

        for(m=0; m<LengthH; ++m)        {

            tmp = ACCESSC_DH(c_in, firstCin, LengthCin,
                        (cfactor*k+(step_factor*(1-m))),bc);

            if (m&1)        /* odd */
                sum += *(H+m) *  tmp;
            else
                sum -= *(H+m) *  tmp;

            }

        *(d_out + count_out) = sum;
        ++count_out;
    }
}
}


/* As WaveThresh reflect, but allows for zero boundary conditions */

int reflect_dh(int n, int lengthC, int bc)
{

if ((n >= 0) && (n < lengthC))
        return(n);
else if (n<0)   {
        if (bc==PERIODIC)       {
                /*
                n = lengthC+n;
                */
                n = n%lengthC + lengthC*((n%lengthC)!=0);
                if (n < 0)      {
                        REprintf("reflect: access error (%d,%d)\n",
                                n,lengthC);
                        REprintf("reflect: left info from right\n");
                        
			error("This should not happen: stopping.\n");
                        }
                else
                        return(n);
                }

        else if (bc==SYMMETRIC) {
                n = -1-n;
                if (n >= lengthC)       {
                        REprintf("reflect: access error (%d,%d)\n",
                                n,lengthC);
			error("This should not happen: stopping.\n");
                        }
                else
                        return(n);
                }

        else    {
                REprintf("reflect: Unknown boundary correction");
                REprintf(" value of %d\n", bc);
		error("This should not happen: stopping.\n");
                }

        }
else    {
        if (bc==PERIODIC)       {
                /*
                printf("periodic extension, was %d (%d) now ",n,lengthC);
                n = n - lengthC;
                */
                n %= lengthC;
                /*
                printf("%d\n", n);
                */
                if (n >= lengthC)       {
                        REprintf("reflect: access error (%d,%d)\n",
                                n,lengthC);
                        REprintf("reflect: right info from left\n");
			error("This should not happen: stopping.\n");
                        }
                else
                        return(n);
                }
        else if (bc==SYMMETRIC) {
                n = 2*lengthC - n - 1;
                if (n<0)        {
                        REprintf("reflect: access error (%d,%d)\n",
                                n,lengthC);
			error("This should not happen: stopping.\n");
                        }
                else
                        return(n);
                }
        else    {
                REprintf("reflect: Unknown boundary correction\n");
		error("This should not happen: stopping.\n");
                }


        }
/* Safety */
REprintf("reflect: SHOULD NOT HAVE REACHED THIS POINT\n");
error("This should not happen: stopping.\n");
return(0); /* for lint only */
}


/* Returns c(n), or 0 if n outside given range */

double access0(double *c, int lengthC, int n)

/* double *c;   data vector      */
/* int lengthC;	length of vector */
/* int n;       index wanted     */

{
    if ((n>=0) && (n<lengthC))
        return(*(c+n));
    else return(0.0);
}


/* As WaveThresh waverecons, but allows for zero boundary conditions */

void waverecons_dh(double *C, double *D, double *H, int *LengthH, int *levels,
        int *firstC, int *lastC, int *offsetC,
	int *firstD, int *lastD, int *offsetD,
	int *type, int *bc, int *error)
/* double *C;              Input data, and the subsequent smoothed data */
/* double *D;              The wavelet coefficients                     */
/* double *H;              The smoothing filter H                       */
/* int *LengthH;           Length of smoothing filter                   */
/* int *levels;            The number of levels in this decomposition   */
/* int *firstC;            The first possible C coef at a given level   */
/* int *lastC;             The last possible C coef at a given level    */
/* int *offsetC;           Offset from C[0] for certain level's coeffs  */
/* int *firstD;            The first possible D coef at a given level   */
/* int *lastD;             The last possible D coef at a given level    */
/* int *offsetD;           Offset from D[0] for certain level's coeffs  */
/* int *type;              The type of wavelet decomposition            */
/* int *bc;                Which boundary handling are we doing         */
/* int *error;             Error code                                   */
{
register int next_level, at_level;
register int verbose;   /* Printing messages, passed in error           */

void conbar_dh(double *c_in, int LengthCin, int firstCin,
           double *d_in, int LengthDin, int firstDin,
           double *H, int LengthH,
           double *c_out, int LengthCout, int firstCout, int lastCout,
	   int type, int bc);

if (*error == 1l)
        verbose = 1;
else
        verbose = 0;

switch(*bc)     {

        case PERIODIC:  /* Periodic boundary conditions */
                if (verbose) Rprintf("Periodic boundary method\n");
                break;

        case SYMMETRIC: /* Symmetric boundary conditions */
                if (verbose) Rprintf("Symmetric boundary method\n");
                break;

        case ZERO: /*Zero boundary conditions */
                if (verbose) Rprintf("Zero boundary method\n");
                break;

        default:        /* The bc must be one of the above */
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

        default:        /* The type must be of one the above */
                if (verbose) Rprintf("Unknown decomposition type\n");
                *error = 2;
                return;
        }

if (verbose) Rprintf("Building level: ");

*error = 0l;

for(next_level = 1; next_level <= *levels; ++next_level)        {


        if (verbose)
                Rprintf("%d ", next_level);

        at_level = next_level - 1;

        conbar_dh( (C+*(offsetC+at_level)),
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


/* As WaveThresh conbar, but allows for zero boundary conditions */

#define CEIL(i) ( ((i)>0) ? ( ((i)+1)/2):((i)/2) )

void conbar_dh(double *c_in, int LengthCin, int firstCin,
           double *d_in, int LengthDin, int firstDin,
           double *H, int LengthH,
           double *c_out, int LengthCout, int firstCout, int lastCout,
	   int type, int bc)
/* int firstCout;           This determines summation over n             */
/* int lastCout;            and this does too                            */
/* int type;                The type of wavelet reconstruction           */
{
    register int n,k;
    register int cfactor;
    double sumC, sumD;

    int reflect_dh(int n, int lengthC, int bc);
    double access0(double *c, int lengthC, int n);

    switch(type)    {

        case WAVELET:   /* Standard wavelets */
            cfactor = 2;
            break;

        case STATION:   /* Stationary wavelets */
            cfactor = 1;
            break;

        default:        /* This should never happen */
            cfactor=0;       /* MAN: added for total cover: shouldn't happen */
            break;
    }


    /* Compute each of the output C for ZERO bcs */

    if(bc==ZERO) {

        for(n=firstCout; n<=lastCout; ++n)      {

        /* We want  n+1-LengthH <= 2*k to start off (n-2k<=LengthH-1) */

            k = CEIL(n+1-LengthH);
            sumC = 0.0;

            while( cfactor*k <= n ) {
                sumC += *(H+n-cfactor*k)*access0(c_in, LengthCin, k-firstCin);
                ++k;
            }

        /* Now do D part */

            k = CEIL(n-1);
            sumD = 0.0;
            while( cfactor*k <= (LengthH +n -2) )   {
                sumD += *(H+1+cfactor*k-n) * access0(d_in, LengthDin,
                                                     k-firstDin);
                ++k;
            }

            if (n & 1)              /* n odd */
                sumC -= sumD;
            else
                sumC += sumD;

            *(c_out+(n-firstCout)) = sumC;
        }
    }

    /* Now for other bcs */
 
    else {

        for(n=firstCout; n<=lastCout; ++n)      {

        /* We want  n+1-LengthH <= 2*k to start off (n-2k<=LengthH-1) */

            k = CEIL(n+1-LengthH);
            sumC = 0.0;

            while( cfactor*k <= n ) {
                sumC += *(H+n-cfactor*k)*ACCESSC_DH(c_in, firstCin, LengthCin, 
                                                    k, bc);
                ++k;
            }

        /* Now do D part */

            k = CEIL(n-1);
            sumD = 0.0;
            while( cfactor*k <= (LengthH +n -2) )   {
                sumD += *(H+1+cfactor*k-n) * ACCESSC_DH(d_in, firstDin,
                                                        LengthDin, k, bc);
                ++k;
            }

            if (n & 1)              /* n odd */
                sumC -= sumD;
            else
                sumC += sumD;

            ACCESSC_DH(c_out, firstCout, LengthCout, n, bc) = sumC;
        }
    } 
}

