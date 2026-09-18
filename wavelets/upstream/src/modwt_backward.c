# include <math.h>

void modwt_backward (double *W, double *V, int *j,
		     int *N, double *h, double *g,
		     int *L, double *Vj
		     )
{
  int t, k, n;
  double k_div;

  for(t=0; t < *N; t++){
    k = t;
    Vj[t] = h[0]*W[k] + g[0]*V[k];
    for(n=1; n < *L; n++){
      k += (int) pow(2,(*j-1));
      k_div = (double) k/(double) *N;
      if(k >= *N) k -= (int) floor(k_div) * *N;
      Vj[t] += h[n]*W[k] + g[n]*V[k];
    }
  }
}
