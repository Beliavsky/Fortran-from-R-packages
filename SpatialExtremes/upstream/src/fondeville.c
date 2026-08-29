#include "header.h"

/* 
The two following functions are tiny modifications of a code fomr yi
Xuan
*/

double fastpnorm_pos(double x)
{
    #include "fastpnorm_data.h"

    if (x >= fastpnorm_max)
      return 1.0;

    const int i = (int)(x * fastpnorm_hinv);
    const double w = (x - fastpnorm_x[i]) * fastpnorm_hinv;
    return w * fastpnorm_y[i + 1] + (1.0 - w) * fastpnorm_y[i];
}

double fastpnorm(double x)
{
    if (x < 0)
        return 1.0 - fastpnorm_pos(-x);

    return fastpnorm_pos(x);
}


double pointEstimate(int j, int *dim, double *shift, double *upper, double *chol){

  #include "lattice_data.h"

  double value = 0, valuea = 0,
    *x = malloc(*dim * sizeof(double)),
    *xa = malloc(*dim * sizeof(double)),//a for antithetic
    *e = malloc(*dim * sizeof(double)),
    *ea = malloc(*dim * sizeof(double)),
    *y = malloc(*dim * sizeof(double)),
    *ya = malloc(*dim * sizeof(double));
    //x[*dim], xa[*dim], e[*dim], ea[*dim], y[*dim], ya[*dim];

  for (int k = 0; k < *dim; k++){
    double dummy = j * sqrtPrimeNumbers[k] + shift[k];
    x[k] = fabs(2 * (dummy - ftrunc(dummy)) - 1);
    xa[k] = 1 - x[k];
  }

  value = valuea = e[0] = ea[0] = fastpnorm(upper[0]);//pnorm(upper[0], 0, 1, 1, 0);
  for (int k=1; k<*dim; k++){
    y[k-1] = qnorm(e[k - 1] * x[k - 1], 0, 1, 1, 0);
    ya[k-1] = qnorm(ea[k - 1] * xa[k - 1], 0, 1, 1, 0);

    //DEAL WITH INFINTE BOUNDS
    /*
      if (!(R_FINITE(y[k-1]))) {
      value = y[k-1] > 0;
      break;
      }
    */

    double tmp = 0, tmpa = 0;
    for (int l=0; l<k; l++){
      tmp += chol[k + l * *dim] * y[l];
      tmpa += chol[k + l * *dim] * ya[l];
    }    

    e[k] = fastpnorm((upper[k] - tmp) / chol[k * (*dim + 1)]);//pnorm((upper[k] - tmp) / chol[k * (*dim + 1)], 0, 1, 1, 0);
    ea[k] = fastpnorm((upper[k] - tmpa) / chol[k * (*dim + 1)]);//pnorm((upper[k] - tmpa) / chol[k * (*dim + 1)], 0, 1, 1, 0);
    value *= e[k];
    valuea *= ea[k];
  }
  
  free(x); free(e); free(y); free(xa); free(ea); free(ya);
  return(0.5 * (value + valuea));
}

void pmvnorm2(int *nMC, int *dim, double *covmat, double *upper, double *est, double *err) {

  /*
    Watch out!!! covmat is assumed to be a correlation matrix
  */  

  //PARAMETERS FOR VARIABLE REORDERING
  double *y = malloc(*dim * sizeof(double)),
    *chol = malloc(*dim * *dim * sizeof(double));
  //double y[*dim], chol[*dim * *dim];
  int pos = 0;
  double min = LONG_MAX;

  for (int i=0;i<(*dim * *dim);i++)
    chol[i] = 0;

  //--------------------------------------------------- VARIABLE REORDERING --------------------------------------------
  //FIND THE FIRST POSITION
  for (int itPos=0; itPos<*dim; itPos++){
    double phi_upper = dnorm(upper[itPos], 0, 1, 0),
      Phi_upper = fastpnorm(upper[itPos]),//pnorm(upper[itPos], 0, 1, 1, 0),
      v = 1 - phi_upper / Phi_upper * (upper[itPos] + phi_upper / Phi_upper);

    if (v < min) {
      min = v;
      pos = itPos;
      y[0] = - phi_upper / Phi_upper;
    }
  }

  //SWITCH POSITIONS 0 and in upper bound
  if (pos != 0){
    {
      double tmp = upper[0];
      upper[0] = upper[pos];
      upper[pos] = tmp;
    }
    
    //SWITCH ROWS AND COLUMNS
    for (int i=0; i<*dim; i++){
      //row
      double tmp = covmat[i * *dim];
      covmat[i * *dim] = covmat[pos + i * *dim];
      covmat[pos + i * *dim] = tmp;
      
      //column
      //double
      tmp = covmat[i];
      covmat[i] = covmat[i + pos * *dim];
      covmat[i + pos * *dim] = tmp;
    }
  }

  //COMPUTE CHOLESKY (due to the normalization it is the 1st column of
  //covmat because it is a correlation matrix)
  for (int i=0 ; i<*dim ; i++)
    chol[i] = covmat[i];

  //ITERATE THE PROCESS
  double *upper_copy = malloc(*dim * sizeof(double));
  //double upper_copy[*dim];
  for (int rec=1; rec<*dim; rec++){
    min = LONG_MAX;
    pos = rec;
    
    //------------------COMPUTE CANDDIDATE FOR NEW UPPERBOUND
    
    for (int itPos=rec; itPos<*dim; itPos++){
      
      double conditionalSum = 0, squaredSum = 0;
      for (int i=0; i<rec; i++){
	conditionalSum += chol[itPos + i * *dim] * y[i];
	squaredSum += chol[itPos + i * *dim] * chol[itPos + i * *dim];
      }
      
      upper_copy[itPos] = (upper[itPos] - conditionalSum) / sqrt(1 - squaredSum);
    }
    
    //FIND THE MINUMUM
    for (int itPos=rec; itPos<*dim; itPos++){
      double phi_upper = dnorm(upper_copy[itPos], 0, 1, 0),
	Phi_upper = fastpnorm(upper_copy[itPos]),//pnorm(upper_copy[itPos], 0, 1, 1, 0),
	v = 1 - phi_upper / Phi_upper * (upper_copy[itPos] + phi_upper / Phi_upper);
      if (v < min) {
	min = v;
	pos = itPos;
	y[rec] = - phi_upper / Phi_upper;
      }
    }

    //------------------SWITCH POSITIONS
    //SWITCH POSITIONS 
    if (rec != pos){
      {
	double tmp = upper[rec];
	upper[rec] = upper[pos];
	upper[pos] = tmp;
      }
      
      //SWITCH ROWS AND COLUMNS
      for (int i=0; i<*dim; i++){
	//row
	double tmp = covmat[rec + i * *dim];
	covmat[rec + i * *dim] = covmat[pos + i * *dim];
	covmat[pos + i * *dim] = tmp;

	//column
	//double
	tmp = covmat[i + rec * *dim];
	covmat[i + rec * *dim] = covmat[i + pos * *dim];
	covmat[i + pos * *dim] = tmp;
	
	//  double
	tmp = chol[rec + i * *dim];
	chol[rec + i * *dim] = chol[pos + i * *dim];
	chol[pos + i * *dim] = tmp;
      }
    } 

    //----------------COMPUTE CHOLESKY TERMS
    double rowSum = 0;
    for (int j=0; j<rec; j++)
      rowSum += chol[rec + j * *dim] * chol[rec + j * *dim];

    chol[rec * (1 + *dim)] = sqrt(1 - rowSum);
    
    for (int i=(rec+1) ; i<*dim ; i++){
      double lProd = 0;
      for (int k=0; k<rec; k++)
	lProd += chol[rec + k * *dim] * chol[i + k * *dim];
      
      chol[i + rec * *dim] = (covmat[i + rec * *dim] - lProd) / chol[rec * (1 + *dim)];
    }
  }
  free(upper_copy);

  //--------------------------------------------------- INTEGRATION ROUTINE --------------------------------------------

  //PARAMETERS FOR INTEGRATION
  int nRep = 12;
  double diff = 0, p = 0, error = 0,
    *shift = malloc(*dim * sizeof(double));
  //shift[*dim];

  for (int i=0; i<nRep; i++){

    GetRNGstate();
    for(int i=0; i<*dim;i++)
      shift[i] = unif_rand();
    PutRNGstate();
  
    double est = 0;
#pragma omp parallel for reduction(+:est)
    for (int j=0; j<*nMC; j++)
      est += pointEstimate(j, dim, shift, upper, chol);

    diff = (est / ((double) *nMC) - p) / ((double) i + 1);
    p += diff;
    error = (i - 1) * error / ((double) i + 1) + diff * diff;
  }
  
  error = 3 * sqrt(error);
  *est = p;
  *err = error;

  free(y); free(chol); free(shift);
  return;
}
