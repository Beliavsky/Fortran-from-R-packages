#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

/****************************************************************************************/
/* This file contains functions to calculate tests of independence and randomness for   */
/*  arbitrary data                                                                      */
/*                                                                                      */
/*                                                                                      */
/****************************************************************************************/




int nchoosek(int n, int k)
{
   if(k==0)
     return 1;
    else
        return (int)( n * nchoosek(n-1,k-1)/k);

    }

int tot_trunc(int d, int trunc)
    {
        int k,s=0;
        for(k=2;k<= trunc;k++)
        {
            s += nchoosek(d,k);
        }
     return s;
    }

int tot_trunc_serial(int p, int trunc)
    {
        int k,s=0;
        for(k=2;k<= trunc;k++)
        {
            s += nchoosek(p-1,k-1);
        }
     return s;
    }

   double maxi(double u, double v)

   {
      if( u> v)
         return u;
      else
         return v;

   }
   int maxint(int u, int v)

   {
      if( u> v)
         return u;
      else
         return v;

   }

   double mini(double u, double v)

   {
      if( u> v)
         return v;
      else
         return u;

   }




void quick_sort(double *t, int lo, int hi)
{
    int i,j;
    double mid;
    double tmp;


    i=lo; j=hi;
    mid = t[  (lo + hi) >> 1 ];


    do
    {
        while (t[i] < mid) i++;
        while (mid < t[j]) j--;

        if (i <= j)
        {

            tmp = t[i];
            t[i++] = t[j];
            t[j--] = tmp;
        }
    } while (i <=j);

    if (lo < j) quick_sort(t, lo, j);
    if (i < hi) quick_sort(t, i, hi);
}





void unique(double *x, int *n, double *values, int *m)
{

int i,n1,k;
double *y = calloc(n[0],sizeof(double));
for(i=0;i<n[0];i++)
y[i]=x[i];
n1 = n[0]-1;

quick_sort(y,0,n1);
k=0;
values[0]=y[0];
for(i=1;i<n[0];i++)
  {
     if(y[i]>values[k])
       {
         k++;
         values[k]=y[i];
       }
  }
m[0]=k+1;
free(y);
}
void prepare_data(double *x, int *n, double *values, int *m, double *Fn, double *fn)
   {
        int i,k, somme;
        double v;
        for(k=0;k<m[0];k++){
            somme=0;
            v = values[k];
            for(i=0;i<n[0];i++){
                    somme += (x[i]<= v);
            }
            Fn[k] = ((double) somme)/((double)n[0]);
        }
        fn[0]=Fn[0];
        for(k=1;k<m[0];k++){
            fn[k] = Fn[k]-Fn[k-1];
        }

   }




   void rank(double *x, double *r, int n)

   {



      int i, j;
      int count;

      for(i=0;i<n;i++)
      {
         count=0;
         for(j=0;j<n;j++)
         {
            if(x[j] <= x[i])
               ++count;
         }
         r[i] = (double)count;

      }
   }



   double mean(double *x, int n)

   {
      int i;
      double sum = 0.0;

      for(i=0;i<n;i++)
         sum += x[i];

      return sum/((double) n);
   }



   double sum(double *x, int n)

   {
      int i;
      double s = 0.0;

      for(i=0;i<n;i++)
         s += x[i];


      return s;
   }

double stdev(double *x, int n)

{
  int i;
  double m, sum = 0.0;

  m = mean(x,n);

  for(i=0;i<n;i++)
    sum += (x[i]-m)*(x[i]-m);

  return sqrt(sum/((double) n));
}




   double maxvec(double *x, int n)

   {
      int i;
      double y, s;

      s = 0.0;
      for(i=0;i<n;i++)
      {
         y = fabs(x[i]);

         if(s <y)
            s = y;
        /* printf("s = %f\n",s); */


      }

      return s;
   }


   void multvec(double *x, double *y, double *xy, int n)

   {
      int i;


      for(i=0;i<n;i++)
         xy[i] = x[i]*y[i];


   }



/* Kendall's tau and Spearman's rho for arbitrary data*/
   void estdep(double *x, double *y, int *n, double *tau, double *rho, double *std)
   {
       int i,j;
       double x0,y0, sum , sum1, sum2, sum3, a1, a2, c1, c2, s1, s2;


       sum  = 0.0;
       sum3 = 0.0;
       s1 = 0.0;
       s2 = 0.0;
       for(i=0;i<n[0];i++)
       {
           sum1 = 0.0;
           sum2 = 0.0;
           x0 = x[i]; y0=y[i];

           for(j=0;j<n[0];j++)
           {
               a1 = ( x[j] <= x0) + ( x[j] < x0);
               a2 = ( y[j] <= y0) + ( y[j] < y0);
               sum  +=  a1*a2;
               sum1 += a1;
               sum2 += a2;

           }
          c1 = ((double)sum1)/((double)n[0]) -1.0;
          c2 = ((double)sum2)/((double)n[0]) -1.0;
           sum3 += c1*c2;
           s1 +=  c1*c1;
           s2 +=  c2*c2;
       }
        s1 = ((double) s1)/((double)n[0]);
        s2 = ((double) s2)/((double)n[0]);
     std[0] = sqrt(s1*s2);
     tau[0] = -1.0+ ((double)sum)/((double)n[0]*n[0]);
     rho[0] = ((double)sum3)/((double)n[0]) /std[0];
     }



   /* Dn function */
    double Dn(double x, double y, double *z, double *f, double *F, int *m)
    {
      int i;
      double zz, A, B, sum = 1.0/3.0;
      for(i=0;i<m[0];i++)
      {
          zz = z[i];
          A = (x<=zz)+(x<zz);
          B = (y<=zz)+(y<zz);
          sum += f[i]*(-0.5*F[i]*(A+B) +(1.0/6.0)*( A*B +(x<zz)*(y<zz)+(x<=zz)*(y<=zz)+f[i]*(A+B+(x<zz)+(y<zz))));
      }
          return(sum);
    }


    void Ifun(double *x, int *n, double *values, int *m, double *I1, double *I1point, double *I4)
    {
      int i, j, k, l;
      double v, xx,yy,zz, somme,somme2, a,am,b,bm, cte=1.0/3.0;
      double *fn = calloc(m[0],sizeof(double));
      double *Fn = calloc(m[0],sizeof(double));


      for(k=0;k<m[0];k++)
         {
            somme=0.0;
            v = values[k];
            for(i=0;i<n[0];i++){
                    somme += (x[i]<= v);

            }
            Fn[k] =  somme/((double)n[0]);
         }

        fn[0]=Fn[0];
        for(k=1;k<m[0];k++){
            fn[k] = Fn[k]-Fn[k-1];
        }


      l=0;
      for(j=0;j<n[0];j++)
      {
          somme2=0.0;
          yy = x[j];
          for(i=0;i<n[0];i++)
          {
              xx=x[i];
              somme=0.0;
              for(k=0;k<m[0];k++)
                 {
                     zz = values[k];
                     a = (xx<=zz);
                     am = (xx<zz);

                     b = (yy<=zz);
                     bm = (yy<zz);

                     somme += fn[k]*( (a+am)*(b+bm)+ am*bm +a*b )/6.0;
                   }
                   I1[l]=somme;

                   l++;
            somme2 += somme;
          }
          I1point[j] = somme2/((double)n[0]);
      }

        l=0;
      for(j=0;j<n[0];j++)
      {
         zz = I1point[j];
          for(i=0;i<n[0];i++)
          {
              I4[l] = I1[l]-zz - I1point[i] + cte;
              l++;
          }
       }
            free(fn);
            free(Fn);
    }


void Ifun0(double *x, int *n, double *values, int *m, double *I1, double *I1point)
    {
      int i, j, k, l;
      double v, xx,yy,zz, somme,somme2, a,am,b,bm;
      double *fn = calloc(m[0],sizeof(double));
      double *Fn = calloc(m[0],sizeof(double));


      for(k=0;k<m[0];k++){
            somme=0.0;
            v = values[k];
            for(i=0;i<n[0];i++){
                    somme += (x[i]<= v);

            }
            Fn[k] =  somme/((double)n[0]);
         }
        fn[0]=Fn[0];
        for(k=1;k<m[0];k++){
            fn[k] = Fn[k]-Fn[k-1];
        }

      l=0;
      for(j=0;j<n[0];j++)
      {
          somme2=0.0;
          yy = x[j];
          for(i=0;i<n[0];i++)
          {
              xx=x[i];
              somme=0.0;
              for(k=0;k<m[0];k++)
                 {
                     zz = values[k];
                     a = (xx<=zz);
                     am = (xx<zz);

                     b = (yy<=zz);
                     bm = (yy<zz);

                     somme += fn[k]*( (a+am)*(b+bm)+ am*bm +a*b )/6.0;
                   }
                   I1[l]=somme;
                   l++;
            somme2 += somme;
          }
          I1point[j] = somme2/((double)n[0]);
      }


            free(fn);
            free(Fn);
    }

     void Amat(int **A, double *cardA, int d, int *trunc)
   /*** Fills the matrix A which identifies the subsets A used to compute the statistics Tn_A  ***/
   /*** A is defined using the binary decomposition. Note that for a given m there are         ***/
   /*** 2^d-1-d subsets. 1 is included in all subsets and the other elements varies            ***/
   /*** so we have sets like {1,3}, {2,3} ,{1,2}, {1,2,4}... They all can be written in binary ***/
   /*** format for example {1,3} will correspond to [1,0,1,0,0,..,0]=                          ***/
   /*** {1,2,4} =[1,1,0,1,0,0,..,0] and so on.                                                 ***/
   /*** Ghoudi and Remillard 2010                                                              ***/
   {
      int i,j,ii,MA,l;
      int **AA, *carda;

      MA=pow(2,d)-1;

      AA = (int **)calloc(MA,sizeof(int *));
       for(j=0; j<MA; j++)
          AA[j] = (int *)calloc(d,sizeof(int));

      carda = (int *)calloc(MA,sizeof(int));


      for(i=0;i<MA;i++)
      {
         AA[i][0]=1;
         carda[i]=0.0;
      /* printf("\n %d ",A[i][0]); */
         ii=i+1;
         for(j=0; j< d; j++)
         {
            AA[i][j]=ii % 2;
            ii=ii/2;
            carda[i] +=AA[i][j];

         }

      }

      l=0;
      for(i=0;i<MA;i++){
            if((carda[i]<=trunc[0])&&(carda[i]>1)){
                cardA[l]=carda[i];
                 for(j=0;j<d;j++)
                    {
                        A[l][j]=AA[i][j];
                    }
                l++;
            }
      }

      free(carda);
      for(j=0;j<MA;j++){free(AA[j]);}
      free(AA);
    }

   void Sn_A(double *IV, int *n, int *d, int *trunc, double *stat, double *cardA, double *M, double *Asets)
   /* This procedure computes the statistics $Sn_A$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of serial dependence so sets A all
   start with 1. i.e A={1,...}. It also fills the matrix M(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, ll, lll, n2, cA,**A;
      double ss;

      cA= tot_trunc(d[0],trunc[0]);

       n2= n[0]*n[0];

      A = (int **)calloc(cA,sizeof(int *));
       for(j=0; j<cA; j++)
          A[j] = (int *)calloc(d[0],sizeof(int));

           Amat(A,cardA,d[0],trunc);


      l=0;
      for(j=0;j<d[0];j++)
         {
           for(i=0;i<cA;i++)
            {
                Asets[l] = A[i][j];
                l++;
            }
      }




     for(k=0;k< cA;k++)
        {
           stat[k]=0.0;
        }

      ll=0;
      lll=0;
      for(j=0; j<n[0]; j++)
      {
         for(i=0; i<n[0]; i++)
         {
            for(k=0;k< cA;k++)
            {
               ss=1.0;

               for(l=0;l<d[0];l++)
                  {
                     if(A[k][l]) ss *= IV[n2*l+lll];
                   }
               stat[k] += ss;
               M[ll]=ss;
               ll++;
            }
            lll++;

         }
      }

      for(k=0;k< cA;k++)
        stat[k]=stat[k]/((double)n[0]);

      for(k=0;k< cA;k++){free(A[k]);}
        free(A);



   }

void Sn(double *I1, double *I1point, int *n, int *d, double *stat, double *J)
   /* This procedure computes the statistic $Sn$ of the multilinear copula.
   It also fills the matrix J(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, n2;
      double ss,ss1,s1,s2,s12,somme,zz,z1,c1,c2,c0;
      double *prodI1point = calloc(n[0],sizeof(double));
      double *I1pointsomme = calloc(n[0],sizeof(double));
      c0 = 1.0/pow(3.0,(double)d[0]);  /*  1/3^(d) */
      c1 = 3.0*c0;  /*  1/3^(d-1) */
      c2 = 3.0*c1; /*  1/3^(d-2) */

      n2= n[0]*n[0];


      for(j=0; j<n[0]; j++)
        {
         ss1 = 1.0;
         ss  = 0.0;
         for(k=0;k< d[0];k++)
            {
               zz = I1point[n[0]*k+j];
               ss1 *= zz;
               ss  += zz;
            }
         prodI1point[j] = ss1;
         I1pointsomme[j] = ss;
        }


      l=0;

      somme = 0.0;

      for(j=0; j<n[0]; j++)
      {
         zz = prodI1point[j];
         for(i=0; i<n[0]; i++)
         {
            ss=1.0;
            s1 = 0.0;
            s2 = 0.0;
            s12 = 0.0;
            for(k=0;k< d[0];k++)
            {
               z1 = I1[n2*k+l];
               ss *= z1;
               s12 += I1point[n[0]*k+i]*I1point[n[0]*k+j];
               s1 += zz*z1/I1point[k*n[0]+j]+ z1*prodI1point[i]/I1point[n[0]*k+i];
               s2 += z1;
            }
            somme += ss - prodI1point[i]-zz + c0 ;
            J[l]= ss -s1 +s2*c1  +(I1pointsomme[i]*I1pointsomme[j]-s12)*c2;
            l++;

         }
      }

      stat[0]=somme/((double)n[0]);  /* Sn  */



   }

void Sn0(double *I1, double *I1point, int *n, int *d, double *Sn)
   /* This procedure computes the statistic $Sn$ of the multilinear copula.
       J is not computed.
   */
   {
      int i,j,k,l, n2;
      double ss,ss1,somme,zz,c0;
      double *prodI1point = calloc(n[0],sizeof(double));
      c0 = 1.0/pow(3.0,(double)d[0]);  /*  1/3^(d) */


      n2= n[0]*n[0];


      for(j=0; j<n[0]; j++)
        {
         ss1 = 1.0;
         ss  = 0.0;
         for(k=0;k< d[0];k++)
            {
               zz = I1point[n[0]*k+j];
               ss1 *= zz;
               ss  += zz;
            }
         prodI1point[j] = ss1;
       }




      somme = 0.0;
      l = 0.0;
      for(j=0; j<n[0]; j++)
      {
         zz = prodI1point[j];
         for(i=0; i<n[0]; i++)
         {
            ss=1.0;
            for(k=0;k< d[0];k++)
            {
               ss *= I1[n2*k+l];

            }
            somme += ss - prodI1point[i]-zz + c0 ;
             l++;
          }
      }

      Sn[0]=somme/((double)n[0]);  /* Sn  */

    free(prodI1point);

   }

void stats_nonserial(double *x, int *n, int *d, int *trunc,double *stat, double *cardA, double *M, double *Asets, double *sn, double *J)
{

int i,k,n2, *m;
double *y, *values, *I4temp, *I1temp, *I1pointtemp, *I4, *I1, *I1point;

n2 = n[0]*n[0];

     y      =  calloc(n[0],sizeof(double));
     values =  calloc(n[0],sizeof(double));
     m      =  calloc(1,sizeof(int));


        I4temp       = calloc(n2,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*d[0],sizeof(double));
    I1point      = calloc(n[0]*d[0],sizeof(double));
    I1           = calloc(n2*d[0],sizeof(double));



     for(k=0;k<d[0];k++)
     {
         for(i=0;i<n[0];i++)
         {
             y[i] = x[k*n[0]+i];
         }
         unique(y,n,values,m);

         Ifun(y, n, values, m, I1temp, I1pointtemp, I4temp);


        for(i=0;i<n2;i++)
         {
             I4[k*n2+i] = I4temp[i];
             I1[k*n2+i] = I1temp[i];
         }
         for(i=0;i<n[0];i++)
         {
             I1point[k*n[0]+i] = I1pointtemp[i];
         }

     }




       Sn(I1,I1point,n,d,sn,J);

       Sn_A(I4, n, d, trunc, stat, cardA, M, Asets);

       free(m);
       free(I1); free(I4);  free(I1point);
       free(I1temp); free(I4temp); free(I1pointtemp);
       free(y); free(values);


}

   void statsim(double *M, double *J, double *xi, int *n, int *cA, double *stat, double *sn)
	{
	     int i,j,k, l, lsn;
	     double a,b,somme, m;

	     m = mean(xi,n[0]);

	     for(i=0;i<n[0];i++)
            xi[i] = xi[i]-m;

         for(k=0;k< cA[0];k++)
           stat[k]=0.0;


        l=0;
           lsn=0;
           somme=0.0;
	     for(j=0;j<n[0];j++)
         {
             a = xi[j];
             for(i=0;i<n[0];i++)
             {
                 b=xi[i];
                 for(k=0;k<cA[0];k++)
                 {
                     stat[k] +=a*b*M[l];
                     l++;
                 }
                 somme += a*b*J[lsn];
                 lsn++;
             }
         }

      for(k=0;k< cA[0];k++)
        {
           stat[k]=stat[k]/((double)n[0]);
        }
        sn[0] = somme/((double)n[0]);
	}


/***********************  functions for serial case **************************/
/* Serial versions of Kendall's tau and Spearman's rho  for arbitrary data*/
    void estdep_serial(double *x, int *n, int *lag, double *tau, double *rho, double *s2)
   {
       int i,j;

       double x0, y0, *y, sum , sum1,sum2, sum3, a1, a2, c1, c2, ss2;

       y = calloc(n[0], sizeof(double));
       /* circular x with lag */
       for(i=0;i<n[0]-lag[0];i++)
           y[i]=x[i+lag[0]];

       for(i=0;i<lag[0];i++)
           y[n[0]-lag[0]+i] = x[i];




       sum  = 0.0;
       sum3 = 0.0;

       ss2 = 0.0;
       for(i=0;i<n[0];i++)
       {
           sum1 = 0.0;
           sum2 = 0.0;
           x0 = x[i]; y0=y[i];

           for(j=0;j<n[0];j++)
           {
               a1 = ( x[j] <= x0) + ( x[j] < x0);
               a2 = ( y[j] <= y0) + ( y[j] < y0);
               sum  +=  a1*a2;
               sum1 += a1;
               sum2 += a2;

           }
          c1 = ((double)sum1)/((double)n[0]) -1.0;
          c2 = ((double)sum2)/((double)n[0]) -1.0;
           sum3 += c1*c2;
           ss2 +=  c2*c2;
       }

        s2[0] = ((double) ss2)/((double)n[0]);

     tau[0] = -1.0+ ((double)sum)/((double)n[0]*n[0]);
     rho[0] = ((double)sum3)/((double)n[0]) /s2[0];

     free(y);
   }

void DDn(double *x, int *n, double *values, int *m, double *IV)
    {
      int i, j, k, l=0;
      double xx,yy,zz, A, B, somme, a,am,b,bm;
      double *fn = calloc(m[0],sizeof(double));
      double *Fn = calloc(m[0],sizeof(double));

      double v;
        for(k=0;k<m[0];k++){
            somme=0.0;
            v = values[k];
            for(i=0;i<n[0];i++){
                    somme += (x[i]<= v);
            }
            Fn[k] = ((double) somme)/((double)n[0]);
        }
        fn[0]=Fn[0];
        for(k=1;k<m[0];k++){
            fn[k] = Fn[k]-Fn[k-1];
        }

      for(j=0;j<n[0];j++)
      {
          yy = x[j];
          for(i=0;i<n[0];i++)
          {
              xx=x[i];
              somme=1.0/3.0;
              for(k=0;k<m[0];k++)
                 {
                     zz = values[k];
                     a = (xx<=zz);
                     am = (xx<zz);
                     A = a+am;
                     b = (yy<=zz);
                     bm = (yy<zz);
                     B = b+bm;
                     somme += fn[k]*(-0.5*Fn[k]*(A+B) +(1.0/6.0)*( A*B +am*bm+a*b+fn[k]*(A+B+am+bm)));
                   }
                   IV[l]=somme;
                   l++;
          }
      }

            free(fn);
            free(Fn);
    }

   void Amatserial(int **A, double *cardA, int p, int *trunc)
   /*** Fills the matrix A which identifies the subsets A used to compute the statistics Tn_A ***/
   /*** A is defined using the binary decomposition. Note that for a given m there are        ***/
   /*** 2^(m-1)-1 subsets. 1 is included in all subsets and the other elements varies         ***/
   /*** so we have sets like {1,3}, {1,2}, {1,2,4}... They all can be written in binary       ***/
   /*** format for example {1,3} will correspond to [1,0,1,0,0,..,0]=                         ***/
   /*** {1,2,4} =[1,1,0,1,0,0,..,0] and so on.                                                ***/
   /*** Ghoudi and Remillard 2010                                                             ***/
   {
      int i,j,ii,MA,l=0;
      int **AA, *carda;

      MA=pow(2,p-1)-1;

      AA = (int **)calloc(MA,sizeof(int *));
       for(j=0; j<MA; j++)
          AA[j] = (int *)calloc(p,sizeof(int));

      carda = (int *)calloc(MA,sizeof(int));


      for(i=0;i<MA;i++)
      {
         AA[i][0]=1;
         carda[i]=1.0;
      /* printf("\n %d ",A[i][0]); */
         ii=i+1;
         for(j=1;j< p;j++)
         {
            AA[i][j]=ii % 2;
            ii=ii/2;
            carda[i] +=AA[i][j];

         }
      /* printf("% d ", cardA[i]); */
      }

      for(i=0;i<MA;i++){
            if(carda[i]<=trunc[0]){
                cardA[l]=carda[i];
                 for(j=0;j<p;j++)
                    {
                        A[l][j]=AA[i][j];
                 }
                l++;
            }
      }

      free(carda);
      for(j=0;j<MA;j++){free(AA[j]);}
      free(AA);
    }

  	void statsim_serial(double *M, double *J, double *xi, int *n, int *p, int *trunc, double *stat, double *sn)
	{
	     int i,j,k,cA, l,lsn;
	     double a,b,somme, m=0.0;

	     m = mean(xi,n[0]);

	     for(i=0;i<n[0];i++)
            xi[i] = xi[i]-m;

         cA= tot_trunc_serial(p[0],trunc[0]);

	     for(k=0;k< cA;k++)
           stat[k]=0.0;


           l=0;
           lsn=0;
           somme=0.0;
	     for(j=0;j<n[0];j++)
         {
             a = xi[j];
             for(i=0;i<n[0];i++)
             {
                 b=xi[i];
                 for(k=0;k<cA;k++)
                 {
                     stat[k] +=a*b*M[l];
                     l++;
                 }
                 somme += a*b*J[lsn];
                 lsn++;
             }
         }
      for(k=0;k< cA;k++)
        {
           stat[k]=stat[k]/((double)n[0]);
        }
        sn[0] = somme/((double)n[0]);
	}

	void Sn_A_serial(double *x, int *n, int *p, int *trunc, double *values, int *m, double *stat, double *cardA, double *M, double *Asets)
   /* This procedure computes the statistics $Sn_A$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of serial dependence so sets A all
   start with 1. i.e A={1,...}. It also fills the matrix M(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, ll, lll, cA,**A;
      double ss, *xlag, **M1;

      cA= tot_trunc_serial(p[0],trunc[0]);



      A = (int **)calloc(cA,sizeof(int *));
       for(j=0; j<cA; j++)
          A[j] = (int *)calloc(p[0],sizeof(int));

      xlag = (double *)calloc(n[0],sizeof(double));

      M1 = (double **)calloc(p[0],sizeof(double *));
      for(k=0; k<p[0]; k++)
        M1[k] = (double *)calloc(n[0]*n[0],sizeof(double));





      Amatserial(A,cardA,p[0],trunc);


      l=0;
      for(j=0;j<p[0];j++)
         {
           for(i=0;i<cA;i++)
            {
                Asets[l] = A[i][j];
                l++;
            }
      }


for(k=0;k<p[0];k++)
{
   for(i=0;i<n[0]-k;i++)
       xlag[i] = x[i+k];


   for(i=0;i<k;i++)
       xlag[n[0]-k+i] = x[i];


  DDn(xlag,n,values,m,M1[k]);


}



     for(k=0;k< cA;k++)
        stat[k]=0.0;

      ll=0;
      lll=0;
      for(j=0; j<n[0]; j++)
      {
         for(i=0; i<n[0]; i++)
         {
            for(k=0;k< cA;k++)
            {
               ss=1.0;
               for(l=0;l<p[0];l++)
                  {
                      if(A[k][l]) ss *= M1[l][lll];

                  }
               stat[k] += ss;
               M[ll]=ss;
               ll++;
            }
            lll++;

         }
      }

      for(k=0;k< cA;k++)
        stat[k]=stat[k]/((double)n[0]);

      for(k=0;k< cA;k++){free(A[k]);}
        free(A);

    free(xlag);
    for(k=0;k<p[0];k++) free(M1[k]) ;
    free(M1);

   }

void Sn_A_serialvec(double *I4, int *n, int *p, int *trunc, double *stat, double *cardA, double *M, double *Asets)
   /* This procedure computes the statistics $Sn_A$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of serial dependence so sets A all
   start with 1. i.e A={1,...}. It also fills the matrix M(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, ll, lll, n2, cA,**A;
      double ss;

      cA= tot_trunc_serial(p[0],trunc[0]);

       n2= n[0]*n[0];

      A = (int **)calloc(cA,sizeof(int *));
       for(j=0; j<cA; j++)
          A[j] = (int *)calloc(p[0],sizeof(int));

           Amatserial(A,cardA,p[0],trunc);


      l=0;
      for(j=0;j<p[0];j++)
         {
           for(i=0;i<cA;i++)
            {
                Asets[l] = A[i][j];
                l++;
            }
      }




     for(k=0;k< cA;k++)
        {
           stat[k]=0.0;
        }

      ll=0;
      lll=0;
      for(j=0; j<n[0]; j++)
      {
         for(i=0; i<n[0]; i++)
         {
            for(k=0;k< cA;k++)
            {
               ss=1.0;

               for(l=0;l<p[0];l++)
                  {
                     if(A[k][l])
                     {
                        ss *= I4[n2*l+lll];
                     }
                   }
               stat[k] += ss;
               M[ll]=ss;
               ll++;
            }
            lll++;

         }
      }

      for(k=0;k< cA;k++)
        stat[k]=stat[k]/((double)n[0]);

      for(k=0;k< cA;k++){free(A[k]);}
        free(A);



   }




void Sn_A_serialvec0(double *IV, int *n, int *p, int *trunc, double *stat, double *cardA, double *Asets)
   /* This procedure computes the statistics $Sn_A$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of serial dependence so sets A all
   start with 1. i.e A={1,...}. It also fills the matrix M(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l,  lll, n2, cA,**A;
      double ss;

      cA= tot_trunc_serial(p[0],trunc[0]);

       n2= n[0]*n[0];

      A = (int **)calloc(cA,sizeof(int *));
       for(j=0; j<cA; j++)
          A[j] = (int *)calloc(p[0],sizeof(int));

           Amatserial(A,cardA,p[0],trunc);


      l=0;
      for(j=0;j<p[0];j++)
         {
           for(i=0;i<cA;i++)
            {
                Asets[l] = A[i][j];
                l++;
            }
      }




     for(k=0;k< cA;k++)
        {
           stat[k]=0.0;
        }


      lll=0;
      for(j=0; j<n[0]; j++)
      {
         for(i=0; i<n[0]; i++)
         {
            for(k=0;k< cA;k++)
            {
               ss=1.0;

               for(l=0;l<p[0];l++)
                  {
                     if(A[k][l])
                     {
                        ss *= IV[n2*l+lll];
                     }
                   }
               stat[k] += ss;

            }
            lll++;

         }
      }

      for(k=0;k< cA;k++)
        stat[k]=stat[k]/((double)n[0]);

      for(k=0;k< cA;k++){free(A[k]);}
        free(A);



   }



  void stats_serial0(double *X, int *n, int *p, int *trunc, double *stat, double *cardA, double *Asets, double *sn)
{


     int n2, d, *m;
     double *x2, *y, *values, *I4temp, *I1temp, *I1pointtemp, *I4, *I1, *I1point;
     int i,k;

      n2 = n[0]*n[0];


     x2     =  calloc(n[0]*2,sizeof(double));
     y      =  calloc(n[0],sizeof(double));
     values =  calloc(n[0],sizeof(double));
     m      =  calloc(1,sizeof(int));

      d = p[0];

     unique(X,n,values,m);


    I4temp       = calloc(n2,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*d,sizeof(double));
    I1point      = calloc(n[0]*d,sizeof(double));
    I1           = calloc(n2*d,sizeof(double));



      for(i=0;i<n[0];i++)
      {
          x2[i]   = X[i];
          x2[i+n[0]] = X[i];


      }

     for(k=0;k<d;k++)
     {
         for(i=0;i<n[0];i++)
         {
             y[i] = x2[n[0]+i-k];
         }


         Ifun(y, n, values, m, I1temp, I1pointtemp, I4temp);


        for(i=0;i<n2;i++)
         {
             I4[k*n2+i] = I4temp[i];
             I1[k*n2+i] = I1temp[i];
         }
         for(i=0;i<n[0];i++)
         {
             I1point[k*n[0]+i] = I1pointtemp[i];
         }

     }

      Sn0(I1,I1point,n,p,sn);


       Sn_A_serialvec0(I4, n, p, trunc, stat, cardA, Asets);




       free(m);
       free(I1); free(I4);  free(I1point);
       free(I1temp); free(I4temp); free(I1pointtemp);
       free(y); free(x2); free(values);

}


void Sn_serial0(double *X, int *n, int *p, double *sn)
{


     int n2, d, *m;
     double *x2, *y, *values, *I1temp, *I1pointtemp, *I1, *I1point;
     int i,k;

      n2 = n[0]*n[0];


     x2     =  calloc(n[0]*2,sizeof(double));
     y      =  calloc(n[0],sizeof(double));
     values =  calloc(n[0],sizeof(double));
     m      =  calloc(1,sizeof(int));

      d = p[0];

     unique(X,n,values,m);


    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I1point      = calloc(n[0]*d,sizeof(double));
    I1           = calloc(n2*d,sizeof(double));



      for(i=0;i<n[0];i++)
      {
          x2[i]   = X[i];
          x2[i+n[0]] = X[i];


      }

     for(k=0;k<d;k++)
     {
         for(i=0;i<n[0];i++)
         {
             y[i] = x2[n[0]+i-k];
         }


         Ifun0(y, n, values, m, I1temp, I1pointtemp);


        for(i=0;i<n2;i++)
         {
            I1[k*n2+i] = I1temp[i];
         }
         for(i=0;i<n[0];i++)
         {
             I1point[k*n[0]+i] = I1pointtemp[i];
         }

     }

      Sn0(I1,I1point,n,p,sn);







       free(m);
       free(I1); free(I1point);
       free(I1temp);  free(I1pointtemp);
       free(y); free(x2); free(values);

}


void snsim_serial(double *I1, double *I1point,  double *xi, int *n, int *p,  double *sn)
	{
	     int i,j,c1, j1, i1,l1,l2, k1,k2, r, ind1, ind2,**A;
	     double a,b,somme,cte = 1.0/3.0;
        double  m,*cardA;

	     m = mean(xi,n[0]);

	     for(i=0;i<n[0];i++)
            xi[i] = xi[i]-m;


        c1 = (int) pow((double)2,(double)p[0])-p[0]-1;

         A = calloc(c1,sizeof(int *));
             for(j=0; j<c1; j++)
             {
                A[j] = calloc(p[0],sizeof(int));
             }

        cardA =  calloc(c1,sizeof(int));

        Amat(A, cardA, p[0], p);






      somme=0.0;
       for(k1=0;k1<c1;k1++)
       {
          l1 = A[k1][0];
          j  = 1;
          while(l1)
          {
               l1+=A[k1][j];
               j++;

          }
          /*  first nonzero component of A[k1] */

          for(k2=0;k2<c1;k2++)
           {
            l2 = A[k2][0];
            j  = 1;
            while(l2)
             {
               l2+=A[k2][j];
               j++;

             }
          /*  first nonzero component A[k2] */
                for(j=0;j<n[0];j++)
                  {
                     a = xi[j];
                     for(i=0;i<n[0];i++)
                       {
                         b=xi[i];
                       for(r=0;r<p[0];r++)
                          {
                             ind1 = A[k1][r];
                             ind2 = A[k2][r];
                             i1 = i+l1-r;
                             if(i1<0)
                               {
                                 i1 = n[0]+i1;
                               }
                             j1 = j+l2-r;
                             if(j1<0)
                               {
                                j1= n[0]+j1;
                               }
                              somme += a*b*(  ind1*ind2*I1[n[0]*j1+i1]+(1-ind1)*ind2*(I1point[j1]-cte)+
                              (1-ind2)*ind1*(I1point[11]-cte)+(1-ind1)*(1-ind2)*cte);
                          }

                        }
                  }
            }
       }



        sn[0] = somme/((double)n[0]);
        free(cardA);
        for(j=0;j<c1;j++)
        {
           free(A[j]);
        }
        free(A);
	}



   void Sn_serial(double *I4, double *I1, double *I1point, int *n, int *p, double *sn, double *J)
   /* This procedure computes the statistic $Sn$ of the multilinear copula.
   It also fills the matrix J(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, n2, r, i1, j1, k1, k2, ind1, ind2, m1,*lA,  **A;
      double ss,ss1, somme,somme1,zz,*cardA, c0, cte=1.0/3.0;
      double *prodI1point = calloc(n[0],sizeof(double));
      double *I1pointsomme = calloc(n[0],sizeof(double));

      c0 = 1.0/pow(3.0,(double)p[0]);  /*  1/3^(d) */



      n2 = n[0]*n[0];


      m1 = tot_trunc(p[0],p[0]);



      lA = calloc(m1,sizeof(int));

         A = calloc(m1,sizeof(int *));
             for(j=0; j<m1; j++)
             {
                A[j] = calloc(p[0],sizeof(int));
             }

        cardA =  calloc(m1,sizeof(double));

        Amat(A, cardA, p[0], p);






      /* computation of the first nonzero element of A[k] */
       for(k=0;k<m1;k++)
        {
          l = A[k][0];
          j  = 1;
          while(l==0)
            {
               l += A[k][j];
               j++;
             }
            lA[k] = j-1;

        }

       for(j=0; j<n[0]; j++)
        {
         ss1 = 1.0;
         ss  = 0.0;
         for(k=0;k< p[0];k++)
            {
               zz = I1point[n[0]*k+j];
               ss1 *= zz;
               ss  += zz;
            }
         prodI1point[j] = ss1;
         I1pointsomme[j] = ss;
        }


      l=0;

      somme = 0.0;

      for(j=0; j<n[0]; j++)
      {
         zz = prodI1point[j];
         for(i=0; i<n[0]; i++)
         {
            ss=1.0;

            for(k=0;k< p[0];k++)
            {
              ss *= I1[n2*k+l];
            }
            somme += ss - prodI1point[i]-zz + c0 ;

            somme1 = 0.0;
            for(k1=0;k1<m1;k1++)
            {
               for(k2=0;k2<m1;k2++)
                {
                   ss1 = 1.0;
                   for(r=0;r<p[0];r++)
                    {
                      ind1 = A[k1][r];
                      ind2 = A[k2][r];
                      i1 = i+lA[k1]-r;
                       if(i1<0){ i1 = n[0]+i1; }
                       else if(i1>=n[0]){
                          i1 = i1-n[0];
                       }
                      j1 = j+lA[k2]-r;
                       if(j1<0){ j1 = n[0]+j1; }
                       else if(j1>=n[0]){
                          j1 = j1-n[0];
                        /*printf("out of bound\n"); */
                       }


                      if( (ind1==1) && (ind2==1))
                       {
                        ss1 *= I4[n[0]*j1+i1];
                       } else if( (ind1==1) && (ind2==0) )
                               {
                                  ss1*= (I1point[i1]-cte);
                               }
                        else if( (ind1==0) && (ind2==1) )
                              {
                                ss1 *= (I1point[j1]-cte);
                              }
                              else
                                 {
                                 ss1*= cte;
                                 }
                         }
                      somme1 += ss1;
                }
            }
           J[l]= somme1;

           l++;

         }
      }

     sn[0]=somme/((double)n[0]);  /* Sn  */

         free(lA);
         for(k=0;k<m1;k++)
          {
             free(A[k]);
             }
         free(A);
         free(cardA);

   }


 void stats_serial(double *X, int *n, int *p, int *trunc, double *stat, double *cardA, double *M, double *Asets, double *sn, double *J)
{


     int n2, d, *m;
     double *x2, *y, *values, *I4temp, *I1temp, *I1pointtemp, *I4, *I1, *I1point;
     int i,k;

      n2 = n[0]*n[0];


     x2     =  calloc(n[0]*2,sizeof(double));
     y      =  calloc(n[0],sizeof(double));
     values =  calloc(n[0],sizeof(double));
     m      =  calloc(1,sizeof(int));

      d = p[0];

     unique(X,n,values,m);


    I4temp       = calloc(n2,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*d,sizeof(double));
    I1point      = calloc(n[0]*d,sizeof(double));
    I1           = calloc(n2*d,sizeof(double));



      for(i=0;i<n[0];i++)
      {
          x2[i]   = X[i];
          x2[i+n[0]] = X[i];


      }

     for(k=0;k<d;k++)
     {
         for(i=0;i<n[0];i++)
         {
             y[i] = x2[n[0]+i-k];
         }


         Ifun(y, n, values, m, I1temp, I1pointtemp, I4temp);


        for(i=0;i<n2;i++)
         {
             I4[k*n2+i] = I4temp[i];
             I1[k*n2+i] = I1temp[i];
         }
         for(i=0;i<n[0];i++)
         {
             I1point[k*n[0]+i] = I1pointtemp[i];
         }

     }

      /*Sn(I1,I1point,n,p,sn,J);*/
      Sn_serial(I4,I1,I1point,n,p,sn,J);

       Sn_A_serialvec(I4, n, p, trunc, stat, cardA, M, Asets);




       free(m);
       free(I1); free(I4);  free(I1point);
       free(I1temp); free(I4temp); free(I1pointtemp);
       free(y); free(x2); free(values);

}

void stats_serial_cvm(double *X, int *n, int *p, int *trunc, double *stat, double *cardA, double *M, double *Asets)
{


     int n2, d, *m;
     double *x2, *y, *values, *I4temp, *I4;
     double *I1temp, *I1pointtemp;
     int i,k;

      n2 = n[0]*n[0];


     x2     =  calloc(n[0]*2,sizeof(double));
     y      =  calloc(n[0],sizeof(double));
     values =  calloc(n[0],sizeof(double));
     m      =  calloc(1,sizeof(int));

      d = p[0];

     unique(X,n,values,m);


    I4temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*d,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));




      for(i=0;i<n[0];i++)
      {
          x2[i]   = X[i];
          x2[i+n[0]] = X[i];


      }

     for(k=0;k<d;k++)
     {
         for(i=0;i<n[0];i++)
         {
             y[i] = x2[n[0]+i-k];
         }


         Ifun(y, n, values, m, I1temp, I1pointtemp, I4temp);


        for(i=0;i<n2;i++)
         {
             I4[k*n2+i] = I4temp[i];

         }


     }

      /*Sn(I1,I1point,n,p,sn,J);*/
     /* Sn_serial(I4,I1,I1point,n,p,sn,J);*/

       Sn_A_serialvec(I4, n, p, trunc, stat, cardA, M, Asets);




       free(m);
       free(I4); free(I1temp); free(I4temp); free(I1pointtemp);

       free(y); free(x2); free(values);

}

 void Sn_serial_binmat(int *p, double *p0, double *J)
   /* This procedure computes the matrix for $Sn$ of the multilinear copula when the margin is Bernoulli.
     */
   {
      int k,l, r, b, j, k1, k2, ind1, ind2, m0, m1,*lA,  **A, **A0,*num0, *num1,*puiss2, *modA;
      double s, ss1, somme,*cardA, *cardA0,  cte=1.0/3.0;
      double cte2;

      m0 = tot_trunc_serial(p[0],p[0]);
      m1 = tot_trunc(p[0],p[0]);

     s = sqrt( p0[0]*(1-p0[0]));
     cte2 = (1.0+p0[0])/6.0;

      lA = calloc(m1,sizeof(int));

      puiss2 = calloc(p[0],sizeof(int));

      num0 = calloc(m0,sizeof(int));
      num1 = calloc(m1,sizeof(int));
      modA = calloc(m1,sizeof(int));

      puiss2[0]=1;
      for(r=1;r<p[0];r++)
       {
          puiss2[r]=2*puiss2[r-1];
       }


         A = calloc(m1,sizeof(int *));
             for(k=0; k<m1; k++)
             {
                A[k] = calloc(p[0],sizeof(int));
             }

       A0 = calloc(m0,sizeof(int *));
             for(k=0; k<m0; k++)
             {
                A0[k] = calloc(p[0],sizeof(int));
             }

        cardA =  calloc(m1,sizeof(double));
        cardA0 =  calloc(m0,sizeof(double));

         Amat(A, cardA, p[0], p);

         for(k=0;k<m1;k++){
            somme=A[k][0];
            for(r=1;r<p[0];r++){
               somme += puiss2[r]*A[k][r];
            }
            num1[k] = somme;
         }
         Amatserial(A0, cardA0, p[0], p);

          for(k=0;k<m0;k++){
            somme=A0[k][0];
            for(r=1;r<p[0];r++){
               somme += puiss2[r]*A0[k][r];
            }
            num0[k] = somme;
         }

     l=0;
      for(k1=0;k1<m0;k1++)
      {
        for(k1=0;k1<m0;k1++){
            J[l]=0.0;
            l++;
        }
       }

      /* computation of the first nonzero element of A[k] */


        for(k=0;k<m1;k++)
        {
          l = A[k][0];
          j  = 1;
          while(l==0)
            {
               l += A[k][j];
               j++;
             }
            lA[k] = j-1;

        }

/* computation of the class of A[k] */
       for(k=0;k<m1;k++)
        {
           b = num1[k]/puiss2[lA[k]];
           j=0;
           while(num0[j]!=b)
           {j++;}
           modA[k]=j;

        }


         for(k2=0;k2<m1;k2++)
            {

               for(k1=0;k1<m1;k1++)
                {
                   ss1 = pow(s,cardA[k1]+cardA[k2]);
                   for(r=0;r<p[0];r++)
                    {
                      ind1 = A[k1][r];
                      ind2 = A[k2][r];


                      if( (ind1 + ind2 ==1) )  /* r in only one of A and B */
                          {
                                  ss1*= cte2;

                          } else
                                 {
                                 ss1*= cte;

                                 }
                     }

                    J[modA[k1] + p[0]*modA[k2]] += ss1;


                }

            }




         free(lA);
         for(k=0;k<m1;k++)
          {
             free(A[k]);
             }
         free(A);
         free(cardA);

         for(k=0;k<m0;k++)
          {
             free(A0[k]);
             }
         free(A0);
         free(cardA0);

   }

   void ind(int j, int *n, int *indices)
   {
      int i;
      for(i=0;i<j;i++)
         {
            indices[i] = n[0]+i-j;

          }
     for(i=j;i<n[0];i++)
          {
             indices[i] = i-j;

          }

   }

   void stats_serial_bin(double *X, int *n, int *p, double *p0, double *ZnA, double *J, double *Asets, double *cardA)
   {

    double *values, s, somme, prod, sn;
    int i, j, k, l, **A,  *m, cA, *indices;

     values = calloc(2,sizeof(double));
      m     = calloc(1,sizeof(int));

     indices     = calloc(n[0],sizeof(int));
      sn = sqrt((double)n[0]);
     unique(X,n,values,m);

     somme = 0.0;
     for(i=0;i<n[0];i++)
     {
        somme += (X[i]==values[0]);
     }
     p0[0] = somme /((double) n[0]);


     cA = tot_trunc_serial(p[0],p[0]);


     s = p0[0]*(1-p0[0]);





         A = calloc(cA,sizeof(int *));
             for(k=0; k<cA; k++)
             {
                A[k] = calloc(p[0],sizeof(int));
             }


Amatserial(A, cardA, p[0], p);

        for(k=0;k<cA;k++)
         {

           somme = 0.0;
            for(i=0;i<n[0];i++)
             {
              prod =1.0;
             for(j=0;j<p[0];j++)
              {
                ind(j,n,indices);
                if(A[k][j]){
                   prod *= ( (X[indices[i]]==values[0])-p0[0]);
                }
              }
             somme += prod;
            }
            ZnA[k] = somme/pow(s,0.5*cardA[k])/sn;
           }



        Sn_serial_binmat(p, p0, J);



      l=0;
      for(j=0;j<p[0];j++)
         {
           for(i=0;i<cA;i++)
            {

                Asets[l] = A[i][j];
                l++;
            }
      }



    for(k=0;k<cA;k++)
          {
             free(A[k]);
             }
         free(A);

   }

   void IfunVectors(double **x, int *n, int *d, double **values, int *m, double *I1, double *I1point, double *I4, double *D00)
    {
      int i, j, k, l,t;
      double v, xx,yy,zz, somme,somme2, a,am,b,bm, prod, somme3;
      double **fn = calloc(d[0],sizeof(double*));
      double **Fn = calloc(d[0],sizeof(double*));

    for(k=0;k<d[0];k++)
    {
       fn[k] = calloc(n[0],sizeof(double));
       Fn[k] = calloc(n[0],sizeof(double));
    }


      for(k=0;k<d[0];k++)
      {
        for(j=0;j<m[k];j++)
         {
            somme=0.0;
            v = values[k][j];
            for(i=0;i<n[0];i++){
                    somme += (x[k][i]<= v);

            }
            Fn[k][j] =  somme/((double)n[0]);
         }

        fn[k][0]=Fn[k][0];
        for(j=1;j<m[k];j++)
          {
            fn[k][j] = Fn[k][j]-Fn[k][j-1];
          }
      }



      l=0;
      somme3=0.0;
      for(t=0;t<n[0];t++)
      {
          somme2=0.0;

          for(i=0;i<n[0];i++)
          {

             prod  = 1.0;
             for(k=0;k<d[0];k++)
             {
              xx = x[k][i];
              yy = x[k][t];
              somme=0.0;
               for(j=0;j<m[k];j++)
                 {
                     zz = values[k][j];
                     a = (xx<=zz);
                     am = (xx<zz);

                     b = (yy<=zz);
                     bm = (yy<zz);

                     somme += fn[k][j]*( (a+am)*(b+bm)+ am*bm +a*b )/6.0;
                  }
                   prod*= somme;
             }
                I1[l] = prod;
                l++;
            somme2 += prod;
          }
          I1point[t] = somme2/((double)n[0]);
          somme3 += I1point[t] ;
      }

      somme3 = somme3/((double)n[0]);
      D00[0] = somme3;

        l=0;
      for(t=0;t<n[0];t++)
      {
         zz = I1point[t];
          for(i=0;i<n[0];i++)
          {
              I4[l] = I1[l]-zz - I1point[i] + somme3;
              l++;
          }
       }


      for(k=0;k<d[0];k++)
      {
         free(fn[k]);
         free(Fn[k]);
      }
            free(fn);
            free(Fn);

    }

    void Sn_serialVectors(double *I4, double *I1, double *I1point, double *D00, int *n, int *p, double *sn, double *J)
   /* This procedure computes the statistic $Sn$ of the multilinear copula.
   It also fills the matrix J(k,i,j) that can be used for the multiplier method.
   */
   {
      int i,j,k,l, n2, r, i1, j1, k1, k2, ind1, ind2, m1,*lA,  **A;
      double ss,ss1, somme,somme1,zz,*cardA, c0, cte;
      double *prodI1point = calloc(n[0],sizeof(double));
      double *I1pointsomme = calloc(n[0],sizeof(double));


      cte = D00[0];
      c0=1.0;
      for(j=0;j<p[0];j++)
      {
         c0 *= cte;
      }

      n2 = n[0]*n[0];


      m1 = tot_trunc(p[0],p[0]);



      lA = calloc(m1,sizeof(int));

         A = calloc(m1,sizeof(int *));
             for(j=0; j<m1; j++)
             {
                A[j] = calloc(p[0],sizeof(int));
             }

        cardA =  calloc(m1,sizeof(double));

        Amat(A, cardA, p[0], p);






      /* computation of the first nonzero element of A[k] */
       for(k=0;k<m1;k++)
        {
          l = A[k][0];
          j  = 1;
          while(l==0)
            {
               l += A[k][j];
               j++;
             }
            lA[k] = j-1;

        }

       for(j=0; j<n[0]; j++)
        {
         ss1 = 1.0;
         ss  = 0.0;
         for(k=0;k< p[0];k++)
            {
               zz = I1point[n[0]*k+j];
               ss1 *= zz;
               ss  += zz;
            }
         prodI1point[j] = ss1;
         I1pointsomme[j] = ss;
        }


      l=0;

      somme = 0.0;

      for(j=0; j<n[0]; j++)
      {
         zz = prodI1point[j];
         for(i=0; i<n[0]; i++)
         {
            ss=1.0;

            for(k=0;k< p[0];k++)
            {
              ss *= I1[n2*k+l];
            }
            somme += ss - prodI1point[i]-zz + c0 ;

            somme1 = 0.0;
            for(k1=0;k1<m1;k1++)
            {
               for(k2=0;k2<m1;k2++)
                {
                   ss1 = 1.0;
                   for(r=0;r<p[0];r++)
                    {
                      ind1 = A[k1][r];
                      ind2 = A[k2][r];
                      i1 = i+lA[k1]-r;
                       if(i1<0){ i1 = n[0]+i1; }
                       else if(i1>=n[0]){
                          i1 = i1-n[0];
                       }
                      j1 = j+lA[k2]-r;
                       if(j1<0){ j1 = n[0]+j1; }
                       else if(j1>=n[0]){
                          j1 = j1-n[0];
                        /*printf("out of bound\n"); */
                       }


                      if( (ind1==1) && (ind2==1))
                       {
                        ss1 *= I4[n[0]*j1+i1];
                       } else if( (ind1==1) && (ind2==0) )
                               {
                                  ss1*= (I1point[i1]-cte);
                               }
                        else if( (ind1==0) && (ind2==1) )
                              {
                                ss1 *= (I1point[j1]-cte);
                              }
                              else
                                 {
                                 ss1*= cte;
                                 }
                         }
                      somme1 += ss1;
                }
            }
           J[l]= somme1;

           l++;

         }
      }

     sn[0]=somme/((double)n[0]);  /* Sn  */

         free(lA);
         for(k=0;k<m1;k++)
          {
             free(A[k]);
             }
         free(A);
         free(cardA);

   }

        void stats_serialVectors(double *x, int *n, int *d, int *p, int *trunc, double *stat, double *cardA, double *M, double *Asets, double *sn, double *J)
{


     int n2, *m,*mm;
     double **X, **x2, **y, **values, *I4temp, *I1temp, *I1pointtemp, *I4, *I1, *I1point, *D00;
     int i,j,k,l;

      n2 = n[0]*n[0];

     X      = calloc(d[0],sizeof(double*));
     x2     = calloc(d[0],sizeof(double*));
     values = calloc(d[0],sizeof(double*));
     y      = calloc(d[0],sizeof(double*));
     D00    = calloc(1,sizeof(double));

     for(k=0;k<d[0];k++)
      {
          X[k]      = calloc(n[0],sizeof(double));
          values[k] = calloc(n[0],sizeof(double));
          x2[k]     = calloc(n[0]*2,sizeof(double));
          y[k]      = calloc(n[0],sizeof(double));
      }



      l=0;
      for(k=0;k<d[0];k++)
      {
         for(i=0;i<n[0];i++)
         {
            X[k][i]= x[l];
            l++;
         }
      }



     m      =  calloc(d[0],sizeof(int));
     mm     =  calloc(1,sizeof(int));


     for(k=0;k<d[0];k++)
     {
         unique(X[k],n,values[k],mm);
         m[k] = mm[0];

     }



    I4temp       = calloc(n2,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*p[0],sizeof(double));
    I1point      = calloc(n[0]*p[0],sizeof(double));
    I1           = calloc(n2*p[0], sizeof(double));


    for(k=0;k<d[0];k++)
     {
      for(i=0;i<n[0];i++)
        {
          x2[k][i]      = X[k][i];
          x2[k][i+n[0]] = X[k][i];
        }
     }



 for(j=0;j<p[0];j++)
     {
       for(i=0;i<n[0];i++)
         {
           for(k=0;k<d[0];k++)
            {
                y[k][i] = x2[k][n[0]+i-j];
            }
         }


         IfunVectors(y, n, d, values, m, I1temp, I1pointtemp, I4temp,D00);



        for(i=0;i<n2;i++)
         {
             I4[j*n2+i] = I4temp[i];
             I1[j*n2+i] = I1temp[i];
        }

         for(i=0;i<n[0];i++)
         {
             I1point[j*n[0]+i] = I1pointtemp[i];
         }

     }


     Sn_serialVectors(I4,I1,I1point,D00,n,p,sn,J);

     Sn_A_serialvec(I4, n, p, trunc, stat, cardA, M, Asets);




       free(m); free(mm);
       free(I1); free(I4);  free(I1point);
       free(I1temp); free(I4temp); free(I1pointtemp); free(D00);

       for(k=0;k<d[0];k++)
         {
         free(X[k]); free(values[k]); free(x2[k]); free(y[k]);
         }
       free(values); free(X); free(x2); free(y);

}


void stats_serialVectors_cvm(double *x, int *n, int *d, int *p, int *trunc, double *stat, double *cardA, double *M, double *Asets)
{


     int n2, *m,*mm;
     double **X, **x2, **y, **values, *I4temp, *I1temp, *I1pointtemp, *I4, *D00;
     int i,j,k,l;

      n2 = n[0]*n[0];

     X      = calloc(d[0],sizeof(double*));
     x2     = calloc(d[0],sizeof(double*));
     values = calloc(d[0],sizeof(double*));
     y      = calloc(d[0],sizeof(double*));
     D00    = calloc(1,sizeof(double));

     for(k=0;k<d[0];k++)
      {
          X[k]      = calloc(n[0],sizeof(double));
          values[k] = calloc(n[0],sizeof(double));
          x2[k]     = calloc(n[0]*2,sizeof(double));
          y[k]      = calloc(n[0],sizeof(double));
      }



      l=0;
      for(k=0;k<d[0];k++)
      {
         for(i=0;i<n[0];i++)
         {
            X[k][i]= x[l];
            l++;
         }
      }



     m      =  calloc(d[0],sizeof(int));
     mm     =  calloc(1,sizeof(int));


     for(k=0;k<d[0];k++)
     {
         unique(X[k],n,values[k],mm);
         m[k] = mm[0];

     }



    I4temp       = calloc(n2,sizeof(double));
    I1pointtemp  = calloc(n[0],sizeof(double));
    I1temp       = calloc(n2,sizeof(double));
    I4           = calloc(n2*p[0],sizeof(double));




    for(k=0;k<d[0];k++)
     {
      for(i=0;i<n[0];i++)
        {
          x2[k][i]      = X[k][i];
          x2[k][i+n[0]] = X[k][i];
        }
     }


 for(j=0;j<p[0];j++)
     {
       for(i=0;i<n[0];i++)
         {
           for(k=0;k<d[0];k++)
            {
                y[k][i] = x2[k][n[0]+i-j];
            }
         }


         IfunVectors(y, n, d, values, m, I1temp, I1pointtemp, I4temp,D00);



        for(i=0;i<n2;i++)
         {
             I4[j*n2+i] = I4temp[i];
         }



     }




     Sn_A_serialvec(I4, n, p, trunc, stat, cardA, M, Asets);




       free(m); free(mm);
       free(I4); free(I1temp); free(I4temp); free(I1pointtemp); free(D00);

       for(k=0;k<d[0];k++)
         {
         free(X[k]); free(values[k]); free(x2[k]); free(y[k]);
         }
       free(values); free(X); free(x2); free(y);

}



/****************************************************************************************/
/* This file contains functions to compute dependence measures used in                  */
/* tests of independence and randomness for  arbitrary data                             */
/*                                                                                      */
/*  By Bruno Remillard, October 1, 2022                                                 */
/****************************************************************************************/


/*  Inverse distribution function of the standard Gaussian r.v.  from Algorithm AS 241 */


double Ninv2(double u)
{
  double plow, phigh, q, r, *a, *b, *c, *d;

  a =(double *)calloc(7,sizeof(double));
  b =(double *)calloc(6,sizeof(double));
  c =(double *)calloc(7,sizeof(double));
  d =(double *)calloc(5,sizeof(double));

  /* Coefficients in rational approximations */
  a[1] = -3.969683028665376e+01;
  a[2] = 2.209460984245205e+02;
  a[3] = -2.759285104469687e+02;
  a[4] = 1.383577518672690e+02;
  a[5] = -3.066479806614716e+01;
  a[6] =  2.506628277459239;

  b[1] = -5.447609879822406e+01;
  b[2] =  1.615858368580409e+02;
  b[3] = -1.556989798598866e+02;
  b[4] =  6.680131188771972e+01;
  b[5] = -1.328068155288572e+01;

  c[1] = -7.784894002430293e-03;
  c[2] = -3.223964580411365e-01;
  c[3] = -2.400758277161838;
  c[4] = -2.549732539343734;
  c[5] = 4.374664141464968;
  c[6] =  2.938163982698783;

  d[1] = 7.784695709041462e-03;
  d[2] = 3.224671290700398e-01;
  d[3] = 2.445134137142996;
  d[4] =  3.754408661907416;

  /* Break-points */

  plow  = 0.02425;
  phigh = 1 - plow;


  if(u < plow)
  {  q = sqrt(-2*log(u));
    return  (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);
  }
  else{
    if(u <= phigh)
    { q = u - 0.5;
      r = q*q;
      return  (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q /(((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
    }
    else
    {
      q = sqrt(-2*log(1-u));
      return -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);
    }
  }
  free(a);
  free(b);
  free(c);
  free(d);
}




double phi(double x)   /* density of standard Gaussian    */

{
  double cte;

  cte = 1.0/sqrt(2.0*3.14159265359);
  return exp(-0.5*x*x)*cte;


}


double LE(double u)   /*Integral of inverse exponential   */

{
  if(u>0.0 )
    return u-u*log(u);
  else
    return 0.0;
}

double LG(double u)   /*Integral of inverse Gaussian    */

{
  if(u>0.0  && u < 1.0 )
    return phi(Ninv2(u));
  else
    return 0.0;
}

void cdfn(double *x, int *n, double *Fn, double *Fn0, double *fn)
{
  int i,k, s0,s1;
  double v;
  for(i=0;i<n[0];i++)
  {
    v = x[i];
    s0 = 0;
    s1 = 0;
    for(k=0;k<n[0];k++)
    {
      s1 += (x[k]<= v);
      s0 += (x[k]< v);
    }

    Fn[i]  = ((double) s1)/((double)n[0]);
    Fn0[i] = ((double) s0)/((double)n[0]);
    fn[i]  = Fn[i]-Fn0[i];
  }
}

void Stat_A(double *x, int *n, int *d, int *trunc, double *statS, double *statG, double *statE,  double *cardA, double *Asets)
  /* This procedure computes the statistics $rho_{A,n}$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of non-serial dependence.
   */
{
  int i,j,k,l, cA,**A;
  double ssG, ssS, ssE, *sgammaS,*sgammaG, *sgammaE, *snG, *snS, *snE;
  double **matx, **matG, **matS, **matE, **matFn, **matfn, **matFn0, **IG, **IS,**IE;

  snE =(double *)calloc(d[0],sizeof(double));

  snS =(double *)calloc(d[0],sizeof(double));

  snG =(double *)calloc(d[0],sizeof(double));

  matS = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matS[j] =(double *)calloc(n[0],sizeof(double));

  matE = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matE[j] =(double *)calloc(n[0],sizeof(double));

  matG = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matG[j] =(double *)calloc(n[0],sizeof(double));





  matx = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matx[j] = (double *)calloc(n[0],sizeof(double));



  matFn = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matFn[j] = (double *)calloc(n[0],sizeof(double));

  matFn0 = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matFn0[j] = (double *)calloc(n[0],sizeof(double));

  matfn = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matfn[j] = (double *)calloc(n[0],sizeof(double));


  IG = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IG[j] = (double *)calloc(n[0],sizeof(double));

  IS = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IS[j] = (double *)calloc(n[0],sizeof(double));

  IE = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IE[j] = (double *)calloc(n[0],sizeof(double));




  l=0;

  for(j=0;j<d[0];j++)
  {
    for(i=0;i<n[0];i++)
    {
      matx[j][i]= x[l];
      l++;
    }
  }



  for(j=0;j<d[0];j++)
  {
    cdfn(matx[j],n,matFn[j],matFn0[j],matfn[j]);
    for(i=0;i<n[0];i++)
    {
      IS[j][i]= 0.5*(matFn[j][i]+matFn0[j][i]-1.0);
      IG[j][i]= (LG(matFn[j][i])-LG(matFn0[j][i]))/matfn[j][i];
      IE[j][i]= (LE(matFn[j][i])-LE(matFn0[j][i]))/matfn[j][i]-1.0;

    }
    snS[j]=stdev(IS[j],n[0]);
    snG[j]=stdev(IG[j],n[0]);
    snE[j]=stdev(IE[j],n[0]);

  }


  cA= tot_trunc(d[0],trunc[0]);



  A = (int **)calloc(cA,sizeof(int *));
  for(j=0; j<cA; j++)
    A[j] = (int *)calloc(d[0],sizeof(int));

  Amat(A,cardA,d[0],trunc);


  sgammaE = (double*)calloc(cA,sizeof(double));
  sgammaG = (double*)calloc(cA,sizeof(double));
  sgammaS = (double*)calloc(cA,sizeof(double));


  l=0;
  for(j=0;j<d[0];j++)
  {
    for(i=0;i<cA;i++)
    {
      Asets[l] = A[i][j];
      l++;
    }
  }




  for(k=0;k< cA;k++)
  {
    statS[k]=0.0;
    statG[k]=0.0;
    statE[k]=0.0;
  }

  for(k=0;k< cA;k++)
  {
    ssS = 1.0;
    ssE = 1.0;
    ssG = 1.0;

    for(j=0;j<d[0];j++)
    {
      if(A[k][j])
      {
        ssS *= snS[j];
        ssG *= snG[j];
        ssE *= snE[j];
      }
    }
    sgammaS[k] += ssS;
    sgammaG[k] += ssG;
    sgammaE[k] += ssE;

  }


  for(i=0; i<n[0]; i++)
  {
    for(k=0;k< cA;k++)
    {
      ssS = 1.0;
      ssE = 1.0;
      ssG = 1.0;

      for(j=0;j<d[0];j++)
      {
        if(A[k][j])
        {
          ssS *= IS[j][i];
          ssG *= IG[j][i];
          ssE *= IE[j][i];
        }
      }
      statS[k] += ssS;
      statG[k] += ssG;
      statE[k] += ssE;

    }


  }


  for(k=0;k< cA;k++)
  {
    statS[k]=statS[k]/((double)n[0])/sgammaS[k];
    statG[k]=statG[k]/((double)n[0])/sgammaG[k];
    statE[k]=statE[k]/((double)n[0])/sgammaE[k];
  }



  free(sgammaE);
  free(sgammaG);
  free(sgammaS);
  free(snE);
  free(snG);
  free(snS);





  for(k=0;k< cA;k++)
    free(A[k]);

  free(A);

  for(j=0;j< d[0];j++)
  {
    free(IE[j]);
    free(IG[j]);
    free(IS[j]);
    free(matFn[j]);
    free(matFn0[j]);
    free(matfn[j]);
    free(matx[j]);
    free(matG[j]);
    free(matE[j]);
    free(matS[j]);

  }

  free(matx);
  free(matG);
  free(matE);
  free(matS);
  free(matFn);
  free(matFn0);
  free(matfn);
  free(IE);
  free(IG);
  free(IS);


}


void Stat_A_serial(double *y, int *n, int *d, int *trunc, double *statS, double *statG, double *statE,  double *cardA, double *Asets)
  /* This procedure computes the statistics $rho_{A,n}$ of the multilinear copula for all subsets A with a max lag of mm.
   This is done in the context of non-serial dependence.
   */
{
  int i,j,k,l, cA, **A;
  double ssG, ssS, ssE, *sgammaS,*sgammaG, *sgammaE, snG, snS, snE, *y2;
  double **matx, **matG, **matS, **matE, **matFn, **matfn, **matFn0, **IG, **IS,**IE;


  matS = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matS[j] =(double *)calloc(n[0],sizeof(double));

  matE = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matE[j] =(double *)calloc(n[0],sizeof(double));

  matG = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matG[j] =(double *)calloc(n[0],sizeof(double));


  IS = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IS[j] =(double *)calloc(n[0],sizeof(double));

  IE = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IE[j] =(double *)calloc(n[0],sizeof(double));

  IG = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    IG[j] =(double *)calloc(n[0],sizeof(double));



  matx = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matx[j] = (double *)calloc(n[0],sizeof(double));



  matFn = (double  **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matFn[j] = (double *)calloc(n[0],sizeof(double));

  matFn0 = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matFn0[j] = (double *)calloc(n[0],sizeof(double));

  matfn = (double **)calloc(d[0],sizeof(double *));
  for(j=0; j<d[0]; j++)
    matfn[j] = (double *)calloc(n[0],sizeof(double));


  y2     =(double *)  calloc(n[0]*2,sizeof(double));


  for(i=0;i<n[0];i++)
  {
    y2[i]   = y[i];
    y2[i+n[0]] = y[i];


  }

  l=0;

  for(j=0;j<d[0];j++)
  {
    for(i=0;i<n[0];i++)
    {
      matx[j][i]=  y2[n[0]+i-j];
      l++;
    }
  }



  for(j=0;j<d[0];j++)
  {
    cdfn(matx[j],n,matFn[j],matFn0[j],matfn[j]);
    for(i=0;i<n[0];i++)
    {
      IS[j][i]= 0.5*(matFn[j][i]+matFn0[j][i]-1.0);
      IG[j][i]= (LG(matFn[j][i])-LG(matFn0[j][i]))/matfn[j][i];
      IE[j][i]= (LE(matFn[j][i])-LE(matFn0[j][i]))/matfn[j][i]-1.0;

    }
    snS=stdev(IS[0],n[0]);
    snG=stdev(IG[0],n[0]);
    snE=stdev(IE[0],n[0]);

  }


  cA= tot_trunc_serial(d[0],trunc[0]);



  A = (int **)calloc(cA,sizeof(int *));
  for(j=0; j<cA; j++)
    A[j] = (int *)calloc(d[0],sizeof(int));

  Amatserial(A,cardA,d[0],trunc);





  sgammaE = (double*)calloc(cA,sizeof(double));
  sgammaG = (double*)calloc(cA,sizeof(double));
  sgammaS = (double*)calloc(cA,sizeof(double));


  l=0;
  for(j=0;j<d[0];j++)
  {
    for(i=0;i<cA;i++)
    {
      Asets[l] = A[i][j];
      l++;
    }
  }




  for(k=0;k< cA;k++)
  {
    statS[k]=0.0;
    statG[k]=0.0;
    statE[k]=0.0;
  }

  for(k=0;k< cA;k++)
  {
    ssS = 1.0;
    ssE = 1.0;
    ssG = 1.0;

    for(j=0;j<d[0];j++)
    {
      if(A[k][j])
      {
        ssS *= snS;
        ssG *= snG;
        ssE *= snE;
      }
    }
    sgammaS[k] += ssS;
    sgammaG[k] += ssG;
    sgammaE[k] += ssE;

  }


  for(i=0; i<n[0]; i++)
  {
    for(k=0;k< cA;k++)
    {
      ssS = 1.0;
      ssE = 1.0;
      ssG = 1.0;

      for(j=0;j<d[0];j++)
      {
        if(A[k][j])
        {
          ssS *= IS[j][i];
          ssG *= IG[j][i];
          ssE *= IE[j][i];
        }
      }
      statS[k] += ssS;
      statG[k] += ssG;
      statE[k] += ssE;

    }


  }


  for(k=0;k< cA;k++)
  {
    statS[k]=statS[k]/((double)n[0])/sgammaS[k];
    statG[k]=statG[k]/((double)n[0])/sgammaG[k];
    statE[k]=statE[k]/((double)n[0])/sgammaE[k];
  }


  free(y2);
  free(sgammaE);
  free(sgammaG);
  free(sgammaS);






  for(k=0;k< cA;k++)
    free(A[k]);

  free(A);

  for(j=0;j< d[0];j++)
  {
    free(IE[j]);
    free(IG[j]);
    free(IS[j]);
    free(matFn[j]);
    free(matFn0[j]);
    free(matfn[j]);
    free(matx[j]);
    free(matG[j]);
    free(matE[j]);
    free(matS[j]);

  }

  free(matx);
  free(matG);
  free(matE);
  free(matS);
  free(matFn);
  free(matFn0);
  free(matfn);
  free(IE);
  free(IG);
  free(IS);


}

