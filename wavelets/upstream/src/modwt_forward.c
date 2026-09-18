# include <math.h>

void modwt_forward (double *V, int *N, int *j,
		    double *h, double *g, int *L,
		    double *Wj, double *Vj
		    )
{
  int t, k, n;
  double k_div;

  for(t=0; t < *N; t++){
    k = t;
    Wj[t] = h[0]*V[k];
    Vj[t] = g[0]*V[k];
    for(n=1; n < *L; n++){
      k -= (int) pow(2,(*j-1));
      k_div = -k/(double) *N;
      if(k < 0) k += (int) ceil(k_div) * *N;
      Wj[t] += h[n]*V[k];
      Vj[t] += g[n]*V[k];
    }
  }
}
