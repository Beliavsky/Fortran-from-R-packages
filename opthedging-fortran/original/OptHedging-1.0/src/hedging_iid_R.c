#include <stdio.h> 
#include <stdlib.h>
#include <math.h> 
#include <R.h>
#include <Rmath.h>


double Cn( double s, double K)
   { 
      if(s>K)
          return s-K;
      else 
          return 0.0;
     
   }

 
   double xbar( double *x, int n)
   {
      int i;
      double sum = 0.0;
   
      for(i=0;i<n;i++)
         sum += x[i];
   
      return sum/((double)n);
   
   }


   double x2bar( double *x, int n)
   {
      int i;
      double sum = 0.0;
   
      for(i=0;i<n;i++)
         sum += x[i]*x[i];
   
      return sum/((double)n);
   
   }


void interpolation1d( double *interpol, double *s, double *F, int *m, double *maxS, double *minS)
           
   {
   
      int j;
      double z, dt;
   
      
      dt = (maxS[0]-minS[0])/(m[0]-1.0);
   
      z = (s[0]-minS[0])/dt;
      j = floor(z);
   
    
   
      if(j < 0)
         interpol[0] = (1.0-z)*F[0]+ z* F[1];
      else
      {
         if(j>= m[0]-2)
             interpol[0] = F[m[0]-1]+ (s[0]-maxS[0])*(F[m[0]-1]-F[m[0]-2])/dt;
         else
             interpol[0] = F[j]+ (z-j)*(F[j+1]-F[j]);
      }
   
   }


/**************************************************************************************/
/* Hedging for i.i.d returns  to be used in R!                                        */
/* Input: minS, maxS, R (Levy returns), m (#pts of the grid), N (sample size)         */
/**************************************************************************************/

   void HedgingIID(double *R, double *T, double *K, double *r, int *put, int *n, int *m,
                   double *maxS, double *minS, int *N, double *S, double *Cvec,
                   double *avec, double *c1)

   {
      int i,j,k;
      double M1,M2,c2,Kp,suma,sumC, *xi, *Ck, *s0, z,dt, **C, **a;



             
      C = (double **)malloc((sizeof(double))*n[0]);
      a = (double **)malloc((sizeof(double))*n[0]);
   
      for(k=0;k<n[0];k++)
      { 
         C[k] = calloc(m[0], sizeof(double));
         a[k] = calloc(m[0], sizeof(double));
      }
   

    /* T: time (years) to maturity ; 1 year = 252 days ; 1 month  = 21 days */
    /* s: actual price                  */
    /* K: strike price                  */
    /* sigma: annual volatility         */
    /* r: annual rate   interest rate   */




      /* periodic values */
            

      Kp = K[0]*exp(-r[0]*T[0]);
      dt = (maxS[0]-minS[0])/(m[0]-1.0);



      for(i=0;i<m[0];i++)
         S[i] = minS[0]+ dt*i ;  /** Grid **/


     
        /*** calculation of  period n-1  **/

      xi =   (double *)malloc(N[0]*sizeof(double));
      Ck =   (double *)malloc(1*sizeof(double));
      s0  =  (double *)malloc(1*sizeof(double));
      for(i=0;i<N[0];i++)
         xi[i] = exp(R[i])-1.0;



         /** Less biased estimates !! **/
      M1  = xbar(xi,N[0]);
      M2  = x2bar(xi,N[0]);

      c1[0] = M1/M2;
      c2 = 1.0-M1*c1[0];

    /*  printf(" Period  %d \n", n[0]-1);  */

      for(i=0;i<m[0];i++)
      {
         suma = 0.0;
         sumC = 0.0;
         for(j=0;j<N[0];j++)
         {
            if(put[0])
               Ck[0] = Cn(Kp,S[i]*(1.0+xi[j]));
            else
               Ck[0] = Cn(S[i]*(1.0+xi[j]),Kp);


            z = (1.0-c1[0]*xi[j])/c2 ;
            suma += xi[j]*Ck[0];
            sumC += z*Ck[0];
         }
         a[n[0]-1][i] = suma/((double)N[0])/M2;
         C[n[0]-1][i] = sumC/((double)N[0]);

      }

      for(k=n[0]-2; k>=0;k--)
      {
       /* printf(" Period  %d \n", k); */


         for(i=0;i<m[0];i++)
         {


            suma = 0.0;
            sumC  = 0.0;
            for(j=0;j<N[0];j++)
            {
               s0[0] = S[i]*(1.0+xi[j]);
               interpolation1d(Ck,s0, C[(k+1)], m,maxS,minS);
               z = (1.0-c1[0]*xi[j])/c2 ;
               suma += xi[j]*Ck[0];
               sumC += z*Ck[0];
            }
            a[k][i]  = suma/((double)N[0])/M2;
            C[k][i]  = sumC/((double)N[0]);


         }
      }

      free(xi);
      

     i = 0;
      for(j=0;j<m[0];j++)
      { 
         for(k=0;k<n[0];k++)
         {
            avec[i] = a[k][j];
            Cvec[i] = C[k][j];
            i ++;
         }
      }



for(k=0;k<n[0];k++)
      { 
        free(C[k]); free(a[k]);
      }
free(a); free(C); free(Ck);
}

