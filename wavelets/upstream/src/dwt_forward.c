# include <math.h>

void dwt_forward (double *V, int *M, double *h,
                  double *g, int *L, double *Wj,
                  double *Vj
		  )
{
  int t, u, n;

  for(t=0; t < *M/2; t++){
    u = 2*t + 1;
    Wj[t] = h[0]*V[u];
    Vj[t] = g[0]*V[u];
    for(n=1; n < *L; n++){
      u--;
      if(u < 0) u = (*M-1);
      Wj[t] += h[n]*V[u];
      Vj[t] += g[n]*V[u];
    }
  }
}
