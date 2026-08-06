#include <stdlib.h>
#include <math.h>

void hosking(double *Xt, int *N, double *vin)
{
  int i, j, t;
  int nrl = 1, nrh = *N-1, ncl = 1, nch = *N-1;
  int nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double *vt, *mt, *Nt, *Dt, *rhot;
  double **phi; /* = dmatrix(1, *N-1, 1, *N-1); */


  vt = (double *) malloc((size_t) ((*N + 2) * sizeof(double)));
  mt = (double *) malloc((size_t) ((*N + 2) * sizeof(double)));
  Nt = (double *) malloc((size_t) ((*N + 2) * sizeof(double)));
  Dt = (double *) malloc((size_t) ((*N + 2) * sizeof(double)));
  rhot = (double *) malloc((size_t) ((*N + 2) * sizeof(double)));

  /*** Begin dmatrix code ***/

  /* allocate pointers to rows */
  phi=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
  /* if (!phi) nrerror("allocation failure 1 in matrix()"); */
  phi += 1;
  phi -= nrl;
  
  /* allocate rows and set pointers to them */
  phi[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
  /* if (!phi[nrl]) nrerror("allocation failure 2 in matrix()"); */
  phi[nrl] += 1;
  phi[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) phi[i]=phi[i-1]+ncol;

  /*** End dmatrix code ***/

  for(i = 1; i <= *N-1; i++)
    for(j = 1; j <= *N-1; j++)
      phi[i][j] = 0.0;

  vt[0] = vin[0];
  Nt[0] = 0.0; Dt[0] = 1.0;
  Xt[0] *= sqrt(vt[0]);

  rhot[0] = 1.0;
  /* phi[1][1] = d / (1.0 - d); */

  for(t = 1; t <= *N-1; t++) {
    rhot[t] = vin[t] / vin[0];
    Nt[t] = rhot[t];
    if(t > 1)
      for(j = 1; j <= t-1; j++)
	Nt[t] -= phi[t-1][j] * rhot[t-j];
    Dt[t] = Dt[t-1] - (Nt[t-1] * Nt[t-1]) / Dt[t-1];
    phi[t][t] = Nt[t] / Dt[t];
    for(j = 1; j <= t-1; j++)
      phi[t][j] = phi[t-1][j] - phi[t][t] * phi[t-1][t-j];
  }

  for(t = 1; t <= *N-1; t++) {
    mt[t] = 0.0;
    for(j = 1; j <= t; j++)
      mt[t] += phi[t][j] * Xt[t-j];
    vt[t] = (1.0 - phi[t][t] * phi[t][t]) * vt[t-1];
    Xt[t] = Xt[t] * sqrt(vt[t]) + mt[t];
  }
  free((char*) (vt));
  free((char*) (mt));
  free((char*) (Nt));
  free((char*) (Dt));
  free((char*) (rhot)); 
  free((char*) (phi[1]));
  free((char*) (phi));
}
