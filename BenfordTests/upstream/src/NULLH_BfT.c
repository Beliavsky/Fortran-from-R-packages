/*
	Do not delete!
	File name		NULLH_BfT.c
	Part of:		   BenfordTests (GNU R contributed package)
	Author:			Dieter William Joenssen
	Copyright:		Dieter William Joenssen
	Email:			Dieter.Joenssen@TU-Ilmenau.de
	Created:		   07 May 2013
	Last Update: 	16 July 2015
	Description:	C Translation of R code to calculate the distribution
					   of certain test statistics for Benford's law
					   under the null-hypothesis
*/
#include <R.h>
#include <Rmath.h>

// gives random frequencies of digits under Null
double rpbenf(double *r_pbenf, int *combfdigits, double *qbenfvals, int *n);
//one wraper for rpbenf to use under R
void c_wraper_rpbenf(double *r_pbenf, int *combfdigits, double *qbenfvals, int *n);

//Calculate test statistics
double compute_chi_square(double *f_obs, double *f_exp, int *n, int *length_f);
double compute_KSD(double *f_obs, double *q_exp, int *n, int *length_f);
double compute_mstar(double *f_obs, double *f_exp, int *n, int *length_f);
double compute_dstar(double *f_obs, double *f_exp, int *n, int *length_f);
double compute_U_square(double *f_obs, double *q_exp, int *n, int *length_f);
double compute_astar(double *f_obs, double mu_benf,int start_digit, int *n, int *length_f);
double compute_J_stat(double *f_obs, double var_benf, double *pbenf_cent, int start_digit, int *n, int *length_f);

//Calculate H0 for the teststatistics
void compute_H0_chi_square(double *H0_chi_square,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_KSD(double *H0_KSD,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_mstar(double *H0_mstar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_dstar(double *H0_dstar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_U_square(double *H0_U_square,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_astar(double *H0_astar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);
void compute_H0_J_stat(double *H0_J_stat,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim);

// Function that computes the relative frequencies of the first digits of a random number satisfying benfords law
double rpbenf(double *r_pbenf, int *combfdigits, double *qbenfvals, int *n)
{
   int i,j;
   double random_x;
   
//   GetRNGstate();
   //set r_pbenf to zeros
   for (j = 0; j< combfdigits[0]; j++)
   {
      r_pbenf[j] = 0;
   }
   for (i = 0; i < n[0]; i++)
   {
      random_x = runif(0,1);
      
      for (j = 0; j< combfdigits[0]; j++)
      {
         if(random_x<=qbenfvals[j])
         {
            r_pbenf[j] = r_pbenf[j]+1;
            break;
         }
      }
   }
   for (j = 0; j< combfdigits[0]; j++)
   {
      r_pbenf[j] = r_pbenf[j]/n[0];
   }
//   PutRNGstate();
   return(*r_pbenf);
}

// wraper for rpbenf so it may be directly accessed from R
void c_wraper_rpbenf(double *r_pbenf, int *combfdigits, double *qbenfvals, int *n)
{
   GetRNGstate();
   *r_pbenf=rpbenf(r_pbenf, combfdigits, qbenfvals, n);
   PutRNGstate();
}

// computes the chi_square test statistic given a certain frequency distribution
double compute_chi_square(double *f_obs, double *f_exp, int *n, int *length_f)
{
   int i;
   double f_dev;
   double chi_square_cum=0;
   
   for(i = 0; i< length_f[0]; i++)
   {
      f_dev=f_obs[i]-f_exp[i];
      chi_square_cum = chi_square_cum + ((f_dev*f_dev)/f_exp[i]);
   }
   chi_square_cum= chi_square_cum*n[0];
   return(chi_square_cum);
}


// computes the H0 distribution for the chi square statistic for the benford distribution
void compute_H0_chi_square(double *H0_chi_square,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {combfdigits= combfdigits*10;}
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_chi_square[i]= compute_chi_square(rpbenf_val,pbenf,n,&combfdigits);
   }
   PutRNGstate();
   
   free(rpbenf_val);
}


// computes the KSD test statistic given a certain frequency distribution
double compute_KSD(double *f_obs, double *q_exp, int *n, int *length_f)
{
   int i;
   //make f_obs the cum sums of f_obs
   // disregarding the last one b/C allways = 1
   for(i = 1;i< (length_f[0]-1);i++)
   {
      f_obs[i]=f_obs[i]+f_obs[i-1];
   }

   // make f_obs the deviations between f_obs and q_exp
   // disregarding the last because allways 1-1=0
   for(i = 0;i< (length_f[0]-1);i++)
   {
      f_obs[i]=f_obs[i]-q_exp[i];
   }
   // find both max and min deviation
   // first value initializes deviations; thus it is skipped in loop
   double max_dev=f_obs[0];
   double min_dev= f_obs[0];
   for(i = 1;i< (length_f[0]-1);i++ )
   {
      if(f_obs[i]>max_dev)
      {max_dev=f_obs[i];continue;}
      if(f_obs[i]<min_dev)
      {min_dev=f_obs[i];}
   }
   //looking for maximum absolute deviation, thus min deveation may be max

   min_dev=fabs(min_dev);

   if(min_dev>max_dev)
   {max_dev=min_dev;}

   //KSD is sqrt(n)* times maximum absolute deviation
   max_dev=max_dev*sqrt((double)n[0]);
   
   return(max_dev);
}

// computes the H0 distribution for the KSD statistic for the benford distribution
void compute_H0_KSD(double *H0_KSD,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {combfdigits= combfdigits*10;}
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_KSD[i]= compute_KSD(rpbenf_val,qbenf,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}


// computes the mstar test statistic given a certain frequency distribution
double compute_mstar(double *f_obs, double *f_exp, int *n, int *length_f)
{
   int i;
   double dev,max_dev,min_dev;
   max_dev=min_dev=f_obs[0]-f_exp[0];
   for(i =1; i<length_f[0];i++)
   {
      dev=f_obs[i]-f_exp[i];
      if(dev>max_dev)
      {max_dev=dev;continue;}
      if(dev<min_dev)
      {min_dev=dev;}
   }
   
   min_dev = fabs(min_dev);
   
   if(min_dev>max_dev)
   {max_dev=min_dev;}
   
   max_dev=max_dev*sqrt(n[0]);   
   
   return(max_dev);
}

// computes the H0 distribution for the mstar statistic for the benford distribution
void compute_H0_mstar(double *H0_mstar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {combfdigits= combfdigits*10;}
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_mstar[i]= compute_mstar(rpbenf_val,pbenf,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}


// computes the dstar test statistic given a certain frequency distribution
double compute_dstar(double *f_obs, double *f_exp, int *n, int *length_f)
{
   int i;
   double dev;
   double cum_dev=0.0;

   for(i =0; i<length_f[0];i++)
   {
      dev=f_obs[i]-f_exp[i];
      cum_dev= cum_dev+(dev*dev);
   }
   
   cum_dev = sqrt(cum_dev)*sqrt(n[0]);
   
   return(cum_dev);
}

// computes the H0 distribution for the mstar statistic for the benford distribution
void compute_H0_dstar(double *H0_dstar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {combfdigits= combfdigits*10;}
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_dstar[i]= compute_dstar(rpbenf_val,pbenf,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}


// computes the U^2 test statistic given a certain frequency distribution
double compute_U_square(double *f_obs, double *q_exp, int *n, int *length_f)
{
   int i;
   double nodigit;
   
   nodigit = length_f[0];
   
   //make f_obs the cum sums of f_obs
   // disregarding the last one b/C allways = 1
   for(i = 1;i< length_f[0];i++)
   {
      f_obs[i]=f_obs[i]+f_obs[i-1];
   }

   // make f_obs the deviations between f_obs and q_exp
   // disregarding the last because allways 1-1=0
   for(i = 0;i< length_f[0];i++)
   {
      f_obs[i]=f_obs[i]-q_exp[i];
   }
   //square the deviations then sum them
   double sum_of_squares = 0;
   double square_of_sums = 0;
   double U_square;
   for(i = 0;i< length_f[0];i++)
   {
      sum_of_squares = sum_of_squares + (f_obs[i]*f_obs[i]);
      //printf("sum_of_squares = %f\n",sum_of_squares);
      square_of_sums = square_of_sums + f_obs[i];
      //printf("square_of_sums = %f\n",square_of_sums);
   }
   sum_of_squares = sum_of_squares * nodigit;
   square_of_sums = square_of_sums*square_of_sums;
   
   //printf("sum_of_squares = %f\n",sum_of_squares);
   //printf("square_of_sums = %f\n",square_of_sums);
   //printf("n[0]/(length_f[0]) = %f\n",(double)n[0]/((double)length_f[0]));
   U_square = ((double)n[0]/(nodigit * nodigit)) * (sum_of_squares-square_of_sums);
   //printf("U_square = %f\n",U_square);
   return(U_square);
}

// computes the H0 distribution for the U^2 statistic for the benford distribution
void compute_H0_U_square(double *H0_U_square,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {combfdigits= combfdigits*10;}
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_U_square[i]= compute_U_square(rpbenf_val,qbenf,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}


// computes the a* test statistic given a certain frequency distribution
double compute_astar(double *f_obs, double mu_benf,int start_digit, int *n, int *length_f)
{
   int i;
   double mu_emp=0.0;
   double astar;
   for(i=0;i<length_f[0];i++)
   {
      mu_emp = mu_emp + (f_obs[i]*(start_digit+i));
   }
   
   astar= fabs(mu_emp-mu_benf)/(length_f[0]+start_digit-1-mu_benf);

   return(astar);
}

// computes the H0 distribution for the a* statistic for the benford distribution
void compute_H0_astar(double *H0_astar,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   int start_digit=1;
   double mu_benf=0.0;
   
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {
      combfdigits= combfdigits*10;
      start_digit = start_digit*10;
   }
   
   for(i=0;i<combfdigits;i++)
   {
      mu_benf = mu_benf + (pbenf[i]*(start_digit+i));
   }
   
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      
      H0_astar[i]= compute_astar(rpbenf_val,mu_benf,start_digit,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}


// computes the J_P^2 test statistic given a certain frequency distribution
double compute_J_stat(double *f_obs, double var_benf, double *pbenf_cent, int start_digit, int *n, int *length_f)
{
   int i;
   double mu_emp=0.0;
   double var_emp=0.0;
   double corr_val =0.0;
   double J_stat;
   double lenfless = (double)(length_f[0]-1);
   
   //compute empirical mean of frequencies
   for(i=0;i<length_f[0];i++)
   {
      mu_emp = mu_emp + f_obs[i];
   }
   mu_emp = mu_emp/length_f[0];
   
   //center empirical frequencies
   for(i=0;i<length_f[0];i++)
   {
      f_obs[i] = f_obs[i] - mu_emp;
   }
   //calculate empirical variance
   for(i=0;i<length_f[0];i++)
   {
      var_emp = var_emp + (f_obs[i]*f_obs[i]);
   }
   
   //calculate empirical covariance
   for(i=0;i<length_f[0];i++)
   {
      corr_val = corr_val + (f_obs[i]*pbenf_cent[i]);
   }
   //calculate correlation
   corr_val = corr_val/ sqrt(lenfless*var_benf*var_emp);

   if(corr_val<0)
   {J_stat=-1;}
   else
   {J_stat=1;}
   J_stat= J_stat*corr_val*corr_val;
   return(J_stat);
}

// computes the H0 distribution for the Jstat statistic for the benford distribution
void compute_H0_J_stat(double *H0_J_stat,int *digits ,double *pbenf,double *qbenf ,int *n, int *n_sim)
{
   int i;
   int combfdigits=9;
   int start_digit=1;
   double mu_benf=0.0;
   double var_benf = 0.0;
   GetRNGstate();
   
   for(i=1; i<digits[0];i++)
   {
      combfdigits= combfdigits*10;
      start_digit = start_digit*10;
   }
   
   //compute empirical mean of frequencies
   for(i=0;i<combfdigits;i++)
   {
      mu_benf = mu_benf + pbenf[i];
   }
   mu_benf = mu_benf/combfdigits;
   //printf("mu_benf = %f\n",mu_benf);
   //center empirical frequencies
   for(i=0;i<combfdigits;i++)
   {
      pbenf[i] = pbenf[i] - mu_benf;
   }
   //calculate empirical variance
   for(i=0;i<combfdigits;i++)
   {
      var_benf = var_benf + (pbenf[i]*pbenf[i]);
   }
   var_benf = var_benf/(double)(combfdigits-1);
   //printf("var_benf = %f\n",var_benf);
   double *rpbenf_val = malloc(combfdigits * sizeof(double));
   
   for(i=0;i<n_sim[0];i++)
   {
      *rpbenf_val=rpbenf(rpbenf_val,&combfdigits,qbenf, n);
      //(double *f_obs, double var_benf, double *pbenf_cent, int start_digit, int *n, int *length_f)
      H0_J_stat[i]= compute_J_stat(rpbenf_val,var_benf,pbenf,start_digit,n,&combfdigits);
   }
   PutRNGstate();
   free(rpbenf_val);
}
