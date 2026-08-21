#include <math.h>

void c_contfrac(const double *a, const double *b, const int *n, double *f, double *tol)
{
	double TINY = 1e-30;
	double EPS = 2.22044604925031e-16;
	double C,D,Delta;
	int j;

	*f = TINY;
	C = *f;
	D = 0.0;

	for(j=0 ; j < *n ; j++){
		D = b[j] + a[j]*D;
		if(D == 0.0){
			D = TINY;
		}
		C = b[j] + a[j]/C;
		if(C == 0.0){
			C = TINY;
		}
		D = 1.0 / D;
		Delta = C*D;
		*f = (*f) * Delta;
		if( ( (Delta - 1.0) <= EPS) && ((1.0 - Delta) <= EPS)){ 
			*tol = Delta -1.0;
			for(j++ ; j < *n ; j++){
			}
			return;
		}
	}
	*tol = Delta -1.0;
	return;
}

void c_contfrac_complex(const double *ar, const double *ai, const double *br, const double *bi, const int *n, double *fr, double *fi, double *tol)
{
	double TINY = 1e-30;
	double EPS = 2.22044604925031e-16;
	double Cr, Ci, Dr, Di, Deltar, Deltai;
	double jj, tempr, tempi, Cinvr, Cinvi;
	int j;

	*fr = TINY;
	*fi = 0.0;
	Cr = *fr;
	Ci = *fi;
	Dr = 0.0;
	Di = 0.0;

	for(j=0 ; j < *n ; j++){

		tempr = Dr;
		tempi = Di;

		Dr = br[j] + ar[j]*tempr - ai[j]*tempi;
		Di = bi[j] + ar[j]*tempi + ai[j]*tempr;  /* D = b[j] + a*D */

		if( (Dr == 0.0) && (Di == 0.0)){
			Dr = TINY;
		}
		jj = Cr*Cr + Ci*Ci;
		Cinvr = Cr/jj;
		Cinvi = -Ci/jj;   /* Cinv = 1/C */
		Cr = br[j] + ar[j]*Cinvr - ai[j]*Cinvi;
		Ci = bi[j] + ar[j]*Cinvi + ai[j]*Cinvr;  /* C = b[j] + a[j]/C [or C = b + a*Cinv ] */

		if((Cr == 0.0) && (Ci == 0.0)){
			Cr = TINY;
		}

		jj = Dr*Dr + Di*Di;
		Dr = Dr/jj;
		Di = -Di/jj;   /* D=1/D */
		
		Deltar = Cr*Dr - Ci*Di;
		Deltai = Cr*Di + Ci*Dr;  /* Delta = C*D */
		
		tempr = *fr;
		tempi = *fi;

		*fr = tempr*Deltar - tempi*Deltai;
		*fi = tempr*Deltai + tempi*Deltar;   /* f = f*D */

		/* diff = (Deltar-1.0)*(Deltar-1.0) + Deltai*Deltai; 
		   diff = mod(Delta-1)^2 */

		if(  ( (Deltar-1.0) <= EPS) && 
		     ( (1.0-Deltar) <= EPS) &&
		     ( (0.0-Deltai) <= EPS) &&
		     ( (Deltai-0.0) <= EPS)  ){ 
			*tol = sqrt((Deltar-1.0)*(Deltar-1.0) + Deltai*Deltai);
			return;
		}
	}
	
	*tol = sqrt((Deltar-1.0)*(Deltar-1.0) + Deltai*Deltai);
	return;
}

void c_convergents(const double *a, const double *b, const double *b0, const int *n, double *A, double *B){
	
	A[0] = *b0; 
	B[0] = 1.0;
	
	A[1] = b[0]*A[0] + a[0]*1.0;  /* because A_{-1} = 1 */
	B[1] = b[0]*B[0] + a[0]*0.0;  /* because B_{-1} = 0 */
	
	for(int j=2; j < (*n)+1 ; j++){
		A[j] = b[j-1]*A[j-1] + a[j-1]*A[j-2];
		B[j] = b[j-1]*B[j-1] + a[j-1]*B[j-2];
	}
	
}

void c_convergents_complex(const double *ar, const double *ai, const double *br, const double *bi, const double *b0r, const double *b0i, const int *n, double *Ar, double *Ai, double *Br, double *Bi){

	Ar[0] = *b0r;
	Ai[0] = *b0i;

	Br[0] = 1.0;
	Bi[0] = 0.0;   /* B[0] = 1 */

	Ar[1] = br[0]*Ar[0] - bi[0]*Ai[0] + ar[0];
	Ai[1] = br[0]*Ai[0] + bi[0]*Ar[0] + ai[0];  

	Br[1] = br[0]*Br[0] - bi[0]*Bi[0];
	Bi[1] = br[0]*Bi[0] + bi[0]*Br[0];

	for(int j=2 ; j < (*n)+1 ; j++){
		Ar[j] = br[j-1]*Ar[j-1] - bi[j-1]*Ai[j-1] + ar[j-1]*Ar[j-2] - ai[j-1]*Ai[j-2];
		Ai[j] = br[j-1]*Ai[j-1] + bi[j-1]*Ar[j-1] + ar[j-1]*Ai[j-2] + ai[j-1]*Ar[j-2]; 

		Br[j] = br[j-1]*Br[j-1] - bi[j-1]*Bi[j-1] + ar[j-1]*Br[j-2] - ai[j-1]*Bi[j-2];
		Bi[j] = br[j-1]*Bi[j-1] + bi[j-1]*Br[j-1] + ar[j-1]*Bi[j-2] + ai[j-1]*Br[j-2]; 
	}
       
}
