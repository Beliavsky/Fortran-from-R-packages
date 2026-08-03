#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "indgenerror.h"


 /*********************************************************************/
 /* program: calculs_cvm_3d.c                                         */
 /*                                                                   */
 /* Computes Cramer-von Mises Mobius statistics                       */
 /*          for test of independence between 2 or 3 series.          */
 /*          Compute also cross-correlations.                         */
 /* By: Bruno Remillard, Kilani Ghoudi, February 28,  2011            */
 /*                                                                   */
 /*                                                                   */
 /*********************************************************************/



  /* Auxiliary functions*/

   double maxi(double u, double v)

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



   void multvec(double *x, double *y, double *z, int n)

   {
      int i;


      for(i=0;i<n;i++)
         z[i] = x[i]*y[i];


   }



   void kappa_xi(double *kappa)
   {
      kappa[1] = 1.0/36.0;
      kappa[2] = 1.0/4050.0;
      kappa[3] = 8.0/893025.0;
      kappa[4] = 4.0/7441875.0;
      kappa[5] = 128.0/2917512675.0;
      kappa[6] = 122235136.0/27179912769384375.0;

   }

   double kappa_xik(int k, int j)

   {
      double zeta[7], a[7];


      a[1]  = 1.0; a[2] = 2.0; a[3] = 8.0; a[4] = 48; a[5] = 384; a[6] = 3840;  /*  2^(n-1)*(n-1)! */
      zeta[1] = 1.0/6.0; zeta[2] = 1.0/90.0; zeta[3] = 1.0/945.0;
      zeta[4] = 1.0/9450.0; zeta[5] = 1.0/93555.0; zeta[6] = 691.0/638512875.0;    /* zeta(2n)/(Pi0^(2n)) */



      if(k==0)
         return a[j];    /*  cumulants of a chi-square with 1df */
      else
         return kappa_xik(k-1,j)*zeta[j];   /*  cumulants of a xik */

   }

/*============================================*/
   int fact(int n)
      {
      int z;
      if (n == 0)  z=1;
      else
      {
         if (n > 1)  z=n*fact(n-1);
         else z=1;
      }
      return z;
   }
/*============================================*/
   double H(int n, double x)
   {
      int m, in, sigm ;
      double z ;
      in = n/2;
      z=0.0;
      for ( m=0; m<=in; m++)
      {
         sigm=2*(2* (int)(m/2) -m )+1;
         z=z+pow(x,(double) (n-2.0*m))*sigm/
            pow(2.0, (double) m)/fact(m)/fact(n-2*m);
      }
      z=z*fact(n);
      return z;
   }
/*============================================*/
   double zeta(double x)
       {
      int m,mm;
      double z, eps;
      eps=0.000000001;
      mm=(int) pow((x-1.0)*eps, 1.0/(1.0-x))+1;
   /* printf("\n mm= %d \n",mm); */
      z=0.0;
      for (m=1; m<=mm; m++)
      {z=z+1.0/pow((double)  m,x); }
      return z;
   }
/*============================================*/
   double K(int n,int d)
      {
      double z ;
      if (n==1) z=1.0/pow(6.0,d);
      else z=pow(2.0,n-1.0)*fact(n-1)*pow(zeta(2.0*n), d)/pow(Pi0,2.0*d*n);
      return z;
   }

double Phi(double x)    /*Distribution function of the standard Gaussian r.v.      */
   {
      /* A&S formula 7.1.26 */

      double p  =  0.3275911;
      double a1 =  0.254829592;
      double a2 = -0.284496736;
      double a3 =  1.421413741;
      double a4 = -1.453152027;
      double a5 =  1.061405429;



      int sign = 1;
      if (x < 0)
         sign = -1;
      x = fabs(x)/sqrt(2.0);


      double t = 1.0/(1.0 + p*x);
      double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

      return 0.5*(1.0 + sign*y);
   }




/*============================================*/

   double edgeworth_cdf(double *y, double *kappa)

   {
      int j;
      double mu, sigma, z, gam[5], p1, p11, p2, p3, p4;


      mu    = kappa[1];
      sigma =  sqrt(kappa[2]);
      z     = (y[0]-mu)/sigma;

      for(j=1; j<=4;j++)
         gam[j] = kappa[j+2]/pow(sigma, j+2.0);

      p1  = (gam[1]/6.0)* H(2,z) +  (gam[2]/24.0)* H(3,z) + (gam[1]*gam[1]/72.0)*H(5,z) ;
      p11 = (gam[3]/120.0)*H(4,z) + (gam[1]*gam[2]/144.0)* H(6,z);
      p2  = (gam[1]*gam[1]*gam[1]/1296.0)* H(8,z) + (gam[4]/720.0)* H(5,z);
      p3  = (gam[2]*gam[2]/1152.0)* H(7,z) +(gam[1]*gam[3]/720.0)* H(7,z);
      p4  = (gam[1]*gam[1]*gam[2]/1728.0)*H(9,z) + (gam[1]*gam[1]*gam[1]*gam[1]/31104.0)* H(11,z) ;
      return Phi(z)-(p1+p11+p2+p3+p4)*exp(-0.5*z*z)/sqrt(2*Pi0);
   }





/*   Main functions*/

   double F(int n, double t, int p)
   {
     double out=0.0;
     /* if( n == 50)
         return maxi(0.0,mini(F50(t,p),0.9999));
      if( n == 100)
         return maxi(0.0,mini(F100(t,p),0.9999));
      if( n == 250)
         return maxi(0.0,mini(F250(t,p),0.9999));
      if( n == 500)
         return maxi(0.0,mini(F500(t,p),0.9999));
      if( n == 1000)
         return maxi(0.0,mini(F1000(t,p),0.9999));
      else
         return maxi(0.0, mini(F1000(t,p),0.9999));
   		*/

      if( n < 61)
        out=  F50(t,p);
      else
         out= F100(t,p);

      return out;
   }


   double Tln(double *r1, double *r2, int n, int lag)
   {

      int i,j;
      double de, cte, t12;
      double  mat1, mat2;


      de = 2.0*n*(n+1.0);
      cte = (2.0*n+1.0)/(6.0*n);


      t12=0.0;


      for(i=0;i<n;i++)
      {
         for(j=0;j<n;j++)
         {
            mat1 = r1[i]*(r1[i]-1)/de + r1[j]*(r1[j]-1)/de  + cte - maxi(r1[i],r1[j])/(n+1.0);
            mat2 = r2[i+lag]*(r2[i+lag]-1)/de + r2[j+lag]*(r2[j+lag]-1)/de  + cte - maxi(r2[i+lag],r2[j+lag])/(n+1.0);

            t12 += (mat1*mat2);
         }
      }

      return t12/((double)n);
   }





   double T3ln(double *r1, double *r2, double *r3, int n, int lag2, int lag3)

   {

      int i,j;
      double de, cte, t123;
      double  mat1, mat2, mat3;


      de = 2.0*n*(n+1.0);
      cte = (2.0*n+1.0)/(6.0*n);


      t123=0.0;


      for(i=0;i<n;i++)
      {
         for(j=0;j<n;j++)
         {
            mat1 = r1[i]*(r1[i]-1)/de + r1[j]*(r1[j]-1)/de  + cte - maxi(r1[i],r1[j])/(n+1.0);
            mat2 = r2[i+lag2]*(r2[i+lag2]-1)/de + r2[j+lag2]*(r2[j+lag2]-1)/de  + cte - maxi(r2[i+lag2],r2[j+lag2])/(n+1.0);
            mat3 = r3[i+lag3]*(r3[i+lag3]-1)/de + r3[j+lag3]*(r3[j+lag3]-1)/de  + cte - maxi(r3[i+lag3],r3[j+lag3])/(n+1.0);

            t123 += (mat1*mat2*mat3);
         }
      }

      return t123/((double)n);
   }


   void Bias_Tdn(int d,int n, double *bias)

   {
      bias[0]=pow(n-1.0,d)/pow(6.0*n,d) +(n-1.0)/pow(-6.0*n,d)-1.0/pow(6.0,d);
   }


/*** W Corrected for bias of TA!! ***/

   void cvm2d(double *x, double *y, int *pn, int *pl_max, double *T, double *PV, double *W, double *PVW, double *chi, double *PVchi)
   {
      int i,lag, sum3,n,l_max, m1,j;
      double  sum1, sum2, t11,t12, p;
      double *r1, *r2, *bias;
      double cum_W[7], cum_F[7];

      l_max = pl_max[0];
      n = pn[0];
      r1=calloc(n+l_max,sizeof(double));
      r2=calloc(n+l_max,sizeof(double));
      bias=calloc(1,sizeof(double));
      rank(x,r1,n);
      rank(y,r2,n);



      for(i=0;i<l_max;i++){r1[n+i]=r1[i]; r2[n+i]=r2[i];}

      sum1 = 0.0; sum2 = 0.0; sum3 = 0;
      for(lag = 0; lag <= l_max;lag++)   /*positive lags */
      {
         t11 = Tln(r1, r2, n, lag);
         t12 = 90.0*(t11-1.0/36.0);  /* normalized */
         T[sum3] = t12 ;
         p =  1.0-F(n,t12,2);
         PV[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;
      }

      for(lag = 1;lag <= l_max;lag++)   /*negative lags */
      {
         t11 = Tln(r2, r1, n, lag);
         t12 = 90.0*(t11-1.0/36.0); /* normalized */
         T[sum3] = t12;
         p =  1.0-F(n,t12,2);
         PV[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;
      }
      Bias_Tdn(2,n,bias);
      W[0]   = sum1 - ((double) sum3)*bias[0];
      chi[0] = sum2;
      free(r1); free(r2);

      m1 = 2*l_max+1;

      for(j=1;j<=6;j++)
      {
         cum_W[j]  = m1*kappa_xik(2,j);
         cum_F[j]  = 2*m1*kappa_xik(0,j);

      }

      PVchi[0]    = 1.0-edgeworth_cdf(chi,cum_F);


      PVW[0]    = 1.0-edgeworth_cdf(W,cum_W);
   }



/*** W Corrected for bias of TA!! ***/
   void cvm3d(double *x, double *y, double *z, int *pn, int *plag2, int *plag3, double *T12, double *T13,
             double *T23, double *T123, double *PV12, double *PV13, double *PV23, double *PV123, double *W12,
             double *W13, double *W23, double *W123,double *W,
             double *PVW12, double *PVW13, double *PVW23, double *PVW123,double *PVW,
             double *chi12, double *chi13, double *chi23,double *chi123,double *chi,
             double *PVF12, double *PVF13, double *PVF23,double *PVF123,double *PVF)
   {
      int i, l, j, lag, l_max, sum3,n,lag2,lag3,m1,m2;
      double sum1, sum2, t11,t12, p;
      double *r1, *r2, *r3, *bias;
      double cum_F12[7], cum_F123[7],cum_F[7], cum_W12[7], cum_W123[7], cum_W[7];
      double pi2, pii[7];

      n = pn[0];
      lag2 = plag2[0];
      lag3 = plag3[0];
      if(lag2 > 2*lag3)
         l_max = lag2;
      else l_max = 2*lag3;

      r1=calloc(n+l_max,sizeof(double));
      r2=calloc(n+l_max,sizeof(double));
      r3=calloc(n+l_max,sizeof(double));
      bias=calloc(1,sizeof(double));
      rank(x,r1,n);
      rank(y,r2,n);
      rank(z,r3,n);

      for(i=0;i<l_max;i++){r1[n+i]=r1[i]; r2[n+i]=r2[i]; r3[n+i]=r3[i];}




      pi2 = Pi0*Pi0;

      pii[1] = pi2;
      for(j=2;j<=6;j++)
         pii[j] = pii[j-1]*pi2;
      Bias_Tdn(2,n,bias);

      m1 = 2*lag2+1;
      m2 = (2*lag3+1)*(2*lag3+1);

      for(j=1;j<=6;j++)
      {
         cum_W12[j]  = m1*kappa_xik(2,j);
         cum_W123[j] = m2*kappa_xik(3,j)*pii[j];
         cum_W[j]    =  3*cum_W12[j]+ cum_W123[j];
         cum_F12[j]  = 2*m1*kappa_xik(0,j);
         cum_F123[j] = 2*m2*kappa_xik(0,j);
         cum_F[j] =  3*cum_F12[j]+ cum_F123[j];

      }



      /*  A = {1,2}  */
      sum1 = 0.0; sum2 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = Tln(r1, r2, n, lag);
         t12 = 90.0*(t11-1.0/36.0);  /* normalized */
         T12[sum3] = t12 ;
         p =  1.0-F(n,t12,2);
         PV12[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;

      }

      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = Tln(r2, r1, n, lag);
         t12 = 90.0*(t11-1.0/36.0); /* normalized */
         T12[sum3] = t12;
         p =  1.0-F(n,t12,2);
         PV12[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;

      }

      W12[0]   = sum1 -((double) sum3)*bias[0];
      chi12[0] = sum2;


    /*  A = {1,3}  */
      sum1 = 0.0; sum2 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = Tln(r1, r3, n, lag);
         t12 = 90.0*(t11-1.0/36.0);  /* normalized */
         T13[sum3] = t12 ;
         p =  1.0-F(n,t12,2);
         PV13[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;

      }


      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = Tln(r3, r1, n, lag);
         t12 = 90.0*(t11-1.0/36.0); /* normalized */
         T13[sum3] = t12;
         p =  1.0-F(n,t12,2);
         PV13[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;

      }

      W13[0]   = sum1-((double) sum3)*bias[0];
      chi13[0] = sum2;

     /*  A = {2,3}  */

      sum1 = 0.0; sum2 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = Tln(r2, r3, n, lag);
         t12 = 90.0*(t11-1.0/36.0);  /* normalized */
         T23[lag] = t12 ;
         p =  1.0-F(n,t12,2);
         PV23[lag] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;

      }

      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = Tln(r3, r2, n, lag);
         t12 = 90.0*(t11-1.0/36.0); /* normalized */
         T23[sum3] = t12;
         p =  1.0-F(n,t12,2);
         PV23[sum3] = p;
         sum1 += t11;
         sum2 += -2.0*log(p);
         sum3++;
      }

      W23[0]   = sum1-((double) sum3)*bias[0];
      chi23[0] = sum2;



   /*  A = {1,2,3}  */

      Bias_Tdn(3,n,bias);

      sum1 = 0.0; sum2 = 0.0; sum3 = 0;
      for(lag = 0;lag <= lag3; lag++)   /*positive lags */
      {
         for(l=0;l <=lag3;l++)           /*positive lags */
         {
            t11 = T3ln(r1, r2, r3, n, lag,l);
            t12 = 90.0*sqrt(90.0)*(t11-1.0/216.0);  /* normalized */
            T123[sum3] = t12 ;
            p =  1.0-F(n,t12,3);
            PV123[sum3] = p;
            sum1 += t11;
            sum2 += -2.0*log(p);
            sum3++;

         }
      }

      for(lag = 0; lag <=lag3;lag++)   /*positive lags */
      {
         for(l=1;l <=lag3;l++)   /*negative lags */
         {
            t11 = T3ln(r3, r1, r2, n, l, (l+lag));
            t12 = 90.0*sqrt(90.0)*(t11-1.0/216.0);   /* normalized */
            T123[sum3] = t12 ;
            p =  1.0-F(n,t12,3);
            PV123[sum3] = p;
            sum1 += t11;
            sum2 += -2.0*log(p);
            sum3++;

         }
      }

      for(lag=1;lag <=lag3;lag++)   /*negative lags */
      {
         for(l=0;l <=lag3;l++)         /*positive lags */
         {
            t11 = T3ln(r2, r1, r3, n, lag,(l+lag));
            t12 = 90.0*sqrt(90.0)*(t11-1.0/216.0);   /* normalized */
            T123[sum3] = t12 ;
            p =  1.0-F(n,t12,3);
            PV123[sum3] = p;
            sum1 += t11;
            sum2 += -2.0*log(p);
            sum3++;

         }
      }

      for(lag=1;lag <=lag3;lag++)   /*negative lags */
      {
         for(l=1;l <=lag;l++)   /*negative lags */
         {
            t11 = T3ln(r2, r1, r3, n, lag ,(lag-l));
            t12 = 90.0*sqrt(90.0)*(t11-1.0/216.0);   /* normalized */
            T123[sum3] = t12 ;
            p =  1.0-F(n,t12,3);
            PV123[sum3] = p;
            sum1 += t11;
            sum2 += -2.0*log(p);
            sum3++;

         }

         for(l=lag+1; l <= lag3;l++)   /*negative lags */
         {
            t11 = T3ln(r3, r1, r2, n, l,(l-lag));
            t12 = 90.0*sqrt(90.0)*(t11-1.0/216.0);   /* normalized */
            T123[sum3] = t12 ;
            p =  1.0-F(n,t12,3);
            PV123[sum3] = p;
            sum1 += t11;
            sum2 += -2.0*log(p);
            sum3++;

         }
      }


      W123[0]   = pi2*(sum1 -((double) sum3)*bias[0]);


      chi123[0] = sum2;

      W[0] = W12[0]+W13[0]+W23[0]+W123[0];

      chi[0] = chi12[0]+chi13[0]+chi23[0]+chi123[0];


      PVF12[0]    = 1.0-edgeworth_cdf(chi12,cum_F12);
      PVF13[0]    = 1.0-edgeworth_cdf(chi13,cum_F12);
      PVF23[0]    = 1.0-edgeworth_cdf(chi23,cum_F12);
      PVF123[0]   = 1.0-edgeworth_cdf(chi123,cum_F123);
      PVF[0]      = 1.0-edgeworth_cdf(chi,cum_F);

      PVW12[0]    = 1.0-edgeworth_cdf(W12,cum_W12);
      PVW13[0]    = 1.0-edgeworth_cdf(W13,cum_W12);
      PVW23[0]    = 1.0-edgeworth_cdf(W23,cum_W12);
      PVW123[0]   = 1.0-edgeworth_cdf(W123,cum_W123);
      PVW[0]      = 1.0-edgeworth_cdf(W,cum_W);



   }

   double crosscor1(double *r1, double *r2, int n, int lag)
   {

      int i;
      double  m1,m2, mat1, mat2;
      mat1=0.0; mat2=0.0;
      m1 = mean(r1,n);
      m2 = mean(r2,n);

      for(i=0;i<n;i++)
      {
         mat1=mat1+(r1[i]-m1)*(r1[i]-m1);
         mat2=mat2+(r2[i]-m2)*(r2[i]-m2);
      }
      mat1=pow(mat1*mat2,0.5);
      mat2=0.0;
      for(i=0;i<n ;i++)
      {
         mat2=mat2+(r1[i]-m1)*(r2[i+lag]-m2);
      }
      mat2=mat2/mat1;
      return mat2;
   }
   void crosscor2d(double *x, double *y, int *pn, int *pl_max, double *T,  double *W)
   {
      int i,lag, sum3,n,l_max;
      double  sum1,  t11;
      double *r1, *r2;

      l_max = pl_max[0];
      n=pn[0];
      r1=calloc(n+l_max,sizeof(double));
      r2=calloc(n+l_max,sizeof(double));
      for(i=0;i<n ;i++){r1[i]=x[i]; r2[i]=y[i];}
      for(i=0;i<l_max;i++){r1[n+i]=r1[i]; r2[n+i]=r2[i];}

      sum1 = 0.0; sum3 = 0;
      for(lag = 0; lag <= l_max;lag++)   /*positive lags */
      {
         t11 = crosscor1(r1, r2, n, lag);
         T[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }

      for(lag = 1;lag <= l_max;lag++)   /*negative lags */
      {
         t11 = crosscor1(r2, r1, n, lag);
         T[sum3] = t11;
         sum1=sum1+t11*t11;
         sum3++;
      }
      W[0]   = pn[0]*sum1;
   }

   double crosscor3(double *r1, double *r2, double *r3, int n, int lag2, int lag3)

   {

      int i;
      double  m1,m2,m3,mat1, mat2,mat3;
      mat1=0.0; mat2=0.0, mat3=0.0;
      m1=0.0; m2=0.0,m3=0.0;
      for(i=0;i<n;i++)
      {
         m1=m1+r1[i];
         m2=m2+r2[i];
         m3=m3+r3[i];
      }
      m1=m1/n;m2=m2/n;m3=m3/n;
      for(i=0;i<n;i++)
      {
         mat1=mat1+(r1[i]-m1)*(r1[i]-m1);
         mat2=mat2+(r2[i]-m2)*(r2[i]-m2);
         mat3=mat3+(r3[i]-m3)*(r3[i]-m3);
      }
      mat1=mat1/n; mat2=mat2/n; mat3=mat3/n;
      mat1=pow(mat1*mat2*mat3,0.5);
      mat2=0.0;
      for(i=0;i<n ;i++)
      {
         mat2=mat2+(r1[i]-m1)*(r2[i+lag2]-m2)*(r3[i+lag3]-m3);
      }
      mat2=mat2/n;
      mat2=mat2/mat1;
      return mat2;
   }

   void crosscor3d(double *x, double *y, double *z, int *pn, int *plag2, int *plag3, double *T12, double *T13,
                  double *T23, double *T123,  double *W12,
                  double *W13, double *W23, double *W123,double *W)
   {
      int i, l, lag, l_max, sum3,n,lag2,lag3;
      double sum1,  t11;
      double *r1, *r2, *r3;

      n = pn[0];
      lag2 = plag2[0];
      lag3 = plag3[0];
      if(lag2 > 2*lag3)
         l_max = lag2;
      else l_max = 2*lag3;

      r1=calloc(n+l_max,sizeof(double));
      r2=calloc(n+l_max,sizeof(double));
      r3=calloc(n+l_max,sizeof(double));
      for(i=0;i<n ;i++){r1[i]=x[i]; r2[i]=y[i]; r3[i]=z[i];}
      for(i=0;i<l_max;i++){r1[n+i]=r1[i]; r2[n+i]=r2[i]; r3[n+i]=r3[i];}






      /*  A = {1,2}  */
      sum1 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = crosscor1(r1, r2, n, lag);
         T12[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }

      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = crosscor1(r2, r1, n, lag);
         T12[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }

      W12[0]   = n*sum1 ;

    /*  A = {1,3}  */
      sum1 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = crosscor1(r1, r3, n, lag);
         T13[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }

      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = crosscor1(r3, r1, n, lag);
         T13[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }
      W13[0]   = n*sum1 ;
     /*  A = {2,3}  */
      sum1 = 0.0; sum3 = 0;
      for(lag=0;lag <=lag2;lag++)   /*positive lags */
      {
         t11 = crosscor1(r2, r3, n, lag);
         T23[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }

      for(lag=1;lag <=lag2;lag++)   /*negative lags */
      {
         t11 = crosscor1(r3, r2, n, lag);
         T23[sum3] = t11 ;
         sum1=sum1+t11*t11;
         sum3++;
      }
      W23[0]   = n*sum1 ;


   /*  A = {1,2,3}  */

      sum1 = 0.0; sum3 = 0;
      for(lag = 0;lag <= lag3; lag++)   /*positive lags */
      {
         for(l=0;l <=lag3;l++)           /*positive lags */
         {
            t11 = crosscor3(r1, r2, r3, n, lag,l);
            T123[sum3] = t11 ;
            sum1 += t11*t11;
            sum3++;

         }
      }

      for(lag = 0; lag <=lag3;lag++)   /*positive lags */
      {
         for(l=1;l <=lag3;l++)   /*negative lags */
         {
            t11 = crosscor3(r3, r1, r2, n, l, (l+lag));
            T123[sum3] = t11 ;
            sum1 += t11*t11;
            sum3++;

         }
      }

      for(lag=1;lag <=lag3;lag++)   /*negative lags */
      {
         for(l=0;l <=lag3;l++)         /*positive lags */
         {
            t11 = crosscor3(r2, r1, r3, n, lag,(l+lag));
            T123[sum3] = t11 ;
            sum1 += t11*t11;
            sum3++;

         }
      }

      for(lag=1;lag <=lag3;lag++)   /*negative lags */
      {
         for(l=1;l <=lag;l++)   /*negative lags */
         {
            t11 = crosscor3(r2, r1, r3, n, lag ,(lag-l));
            T123[sum3] = t11 ;
            sum1 += t11*t11;
            sum3++;

         }

         for(l=lag+1; l <= lag3;l++)   /*negative lags */
         {
            t11 = crosscor3(r3, r1, r2, n, l,(l-lag));
            T123[sum3] = t11 ;
            sum1 += t11*t11;
            sum3++;

         }
      }


      W123[0]   = n*sum1;

      W[0] = W12[0]+W13[0]+W23[0]+W123[0];
   }




