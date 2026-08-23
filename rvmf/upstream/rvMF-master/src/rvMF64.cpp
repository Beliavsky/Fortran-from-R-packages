#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

#define dg(m,k) ((m>>(30-6*k))&63)            /* gets k-th digit of m (base 64), k=1,...,5    */

static int *P;                                /* Probabilities as an array of 30-bit integers */
static int size,t1,t2,t3,t4,t5,offset;        /* size of P[], limits for table lookups        */
static unsigned short *AA,*BB,*CC,*DD,*EE;    /* Tables for condensed table-lookup            */

/* Creates the 5 tables after array P[size] has been created */
void get5tbls(){
  int i,j,m,na=0,nb=0,nc=0,nd=0,ne=0;//k
  /* get table sizes, calloc */
  for(i=0;i<size;i++){m=P[i];na+=dg(m,1);nb+=dg(m,2);nc+=dg(m,3);nd+=dg(m,4);ne+=dg(m,5);}
  AA=(unsigned short*)calloc(na,sizeof(short));BB=(unsigned short*)calloc(nb,sizeof(short));
  CC=(unsigned short*)calloc(nc,sizeof(short));DD=(unsigned short*)calloc(nd,sizeof(short));
  EE=(unsigned short*)calloc(ne,sizeof(short));
  t1=na<<24; t2=t1+(nb<<18); t3=t2+(nc<<12); t4=t3+(nd<<6); t5=t4+ne;
  na=nb=nc=nd=ne=0;
  /* Fill tables AA,BB,CC,DD,EE */
  for(i=0;i<size;i++){
    m=P[i]; //k=i+offset;
    for(j=0;j<dg(m,1);j++) AA[na+j]=i; na+=dg(m,1);
    for(j=0;j<dg(m,2);j++) BB[nb+j]=i; nb+=dg(m,2);
    for(j=0;j<dg(m,3);j++) CC[nc+j]=i; nc+=dg(m,3);
    for(j=0;j<dg(m,4);j++) DD[nd+j]=i; nd+=dg(m,4);
    for(j=0;j<dg(m,5);j++) EE[ne+j]=i; ne+=dg(m,5);
  }
}

/* Creates Pi Probabilities
 * P's are 30-bit integers, assumed denominator 2^30
 * Fill the array P, and set size, offset variable
 * Assume p=0 if p < 2^(-31)
 */
void PiP(double kappa, int m, double log_conf){ /* S^(m) */
int i,j=-1,imax,last=0;
  double p=1.,t=1.;
  /* If pi(0)>=1/2^30 generate P from 0 */
  if(log_conf<=20.79442){ /* 20.79442 = log(2^(30)) */
  p=t=exp(-log_conf);
    for(i=1;t*2147483648.>1;i++) t*=(2*kappa*(m/2.0+i-1)/(m+i-1)/i);  /* 2147483648 = 2^31 */
  size=i-1; offset=0;// last=i-2;
  /* Given size, calloc and fill P array, (30-bit integers) */
  P=(int*)calloc(size,sizeof(int)); P[0]=exp(-log_conf)*1073741824+.5; /* add .5 for rounding */
  for(i=1;i<size;i++) {p*=(2*kappa*(m/2.0+i-1)/(m+i-1)/i); P[i]=p*1073741824+.5;} /* 1073741824 = 2^30 */
  }

  /* If pi(0)<1/2^30, generate from largest P up, then largest down */
  else{
    imax=ceil((-m-1+2.*kappa+sqrt((m+1-2*kappa)*(m+1-2*kappa)-4*m*(1-kappa)))/2.); /* pi(l) attains max at l=imax */
  p=t=exp(-log_conf+lgamma(m)-lgamma(m/2.0)+imax*log(2*kappa)+lgamma(m/2.0+imax)-lgamma(imax+1)-lgamma(m+imax));
  /* imax to up */
  for(i=imax+1;t*2147483648.>1;i++) t*=(2*kappa*(m/2.0+i-1)/(m+i-1)/i);
  last=i-2; t=p; //j=-1;
  /* imax to down */
  for(i=imax-1;i>=0;i--){t*=((i+1)*(m+i)/(2*kappa)/(m/2.0+i)); if(t*2147483648.<1){j=i;break;} }
  offset=j+1;  size=last-offset+1;
  /* Given size, calloc and fill P array, (30-bit integers) */
  P=(int*)calloc(size,sizeof(int));
  t=p; P[imax-offset]=p*1073741824+.5;
  for(i=imax+1;i<=last;i++){t*=(2*kappa*(m/2.0+i-1)/(m+i-1)/i);P[i-offset]=t*1073741824+.5;}
  t=p;
  for(i=imax-1;i>=offset;i--){t*=((i+1)*(m+i)/(2*kappa)/(m/2.0+i));P[i-offset]=t*1073741824+.5;}
  }
}

// [[Rcpp::export(.rvMF64)]]
NumericVector rvMF64(int n, int p, double kappa, double log_const) {
  PiP(kappa, p-1, log_const); get5tbls();
  NumericVector l(n);
  for(int i = 0; i < n; i++){
    int j = floor(R::runif(0,1)*t5);
    int u = 0;
    if(j<t1) u = AA[j>>24];
    else if(j<t2) u= BB[(j-t1)>>18];
    else if(j<t3) u= CC[(j-t2)>>12];
    else if(j<t4) u= DD[(j-t3)>>6];
    else u= EE[j-t4];
    l[i] = 2.*R::rbeta((p-1)/2.0 + u + offset, (p-1)/2.0)-1.;
  }
  free(AA); free(BB); free(CC); free(DD); free(EE); free(P);
  return l;
}
