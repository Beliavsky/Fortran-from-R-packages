# include <math.h>

void dwt_backward (double *W, double *V, int *M,
		   double *h, double *g, int *L,
		   double *Vj
		   )
{
  int l, m, t, u, i, k, n;

  l = -2;
  m = -1;

  for(t=0; t < *M; t++){
    l += 2;
    m += 2;
    u = t;
    i = 1;
    k = 0;
    Vj[l] = h[i]*W[u] + g[i]*V[u];
    Vj[m] = h[k]*W[u] + g[k]*V[u];
    if(*L > 2){
      for(n=1; n < (*L/2); n++){
	u++;
	if(u >= *M) u = 0;
	i += 2;
	k += 2;
	Vj[l] += h[i]*W[u] + g[i]*V[u];
	Vj[m] += h[k]*W[u] + g[k]*V[u];
      }
    }
  }
}
