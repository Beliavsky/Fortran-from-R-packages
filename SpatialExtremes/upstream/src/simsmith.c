#include "header.h"

void rsmith1d(double *coord, double *center, double *edge, int *nObs,
	      int *nSites, double *var, double *ans){
  /* This function generates random fields for the 1d smith model

     coord: the coordinates of the locations
    center: the center of the compact set - here I use an interval
      edge: the length of the interval
      nObs: the number of observations to be generated
    nSites: the number of locations
       var: the variance of the univariate normal density
       ans: the generated random field */

  const double uBound = M_1_SQRT_2PI / sqrt(*var);

  if (*var <= 0)
    error("The variance should be strictly positive!\n");

  /* We first center the coordinates to avoid repetition of
    unnecessary operations in the while loop */
  for (int i=0;i<*nSites;i++)
    coord[i] -= center[0];

  /* Simulation according to the Schlather methodology. The compact
     set need to be inflated first */
  *edge += 6.92 * sqrt(*var);
  const double lebesgue = *edge;

  GetRNGstate();
  
  for (int i=0;i<*nObs;i++){
    double poisson = 0;
    int nKO = *nSites;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      
      poisson += exp_rand();
      double ipoisson = 1 / poisson, thresh = uBound * ipoisson;
      
      //We simulate points uniformly in [-r/2, r/2]
      double u = *edge * runif(-0.5, 0.5);
            
      nKO = *nSites;
      for (int j=0;j<*nSites;j++){
	//This is the normal density with 0 mean and variance var
	double y = exp(-(coord[j] - u) * (coord[j] - u) / (2 * *var)) * thresh;	  
	ans[i + j * *nObs] = fmax2(y, ans[i + j * *nObs]);
	nKO -= (thresh <= ans[i + j * *nObs]);
      }
    }
  }
 
  PutRNGstate();

  /* Lastly, we multiply by the Lebesgue measure of the dilated
    compact set */
  for (int i=0;i<(*nSites * *nObs);i++)
    ans[i] *= lebesgue;

  return;
}

void rsmith2d(double *coord, double *center, double *edge, int *nObs,
	      int *nSites, int *grid, double *cov11, double *cov12,
	      double *cov22, double *ans){
  /* This function generates random fields for the 2d smith model

     coord: the coordinates of the locations
    center: the center of the compact set - here I use a square
      edge: the length of the edge of the square
      nObs: the number of observations to be generated
      grid: Does coord specifies a grid?
    nSites: the number of locations
     covXX: the parameters of the bivariate normal density
       ans: the generated random field */

  const double det = *cov11 * *cov22 - *cov12 * *cov12,
    uBound = 1 / (M_2PI * sqrt(det)), itwiceDet = 1 / (2 * det);

  if ((det <= 0) || (*cov11 <= 0))
    error("The covariance matrix isn't semi-definite positive!\n");

  /* We first center the coordinates to avoid repetition of
    unnecessary operations in the while loop */
  for (int i=0;i<*nSites;i++){
    coord[i] -= center[0];
    coord[*nSites + i] -= center[1];
  }

  /* Simulation according to the Schlather methodology. The compact
     set need to be inflated first */
  *edge += 6.92 * sqrt(fmax2(*cov11, *cov22));
  double lebesgue = *edge * *edge;

  GetRNGstate();
  
  if (*grid){
    //Simulation part if a grid is used
    for (int i=0;i<*nObs;i++){
      double poisson = 0;
      int nKO = *nSites * *nSites;
    
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */

	poisson += exp_rand();
	double ipoisson = 1 / poisson, thresh = uBound * ipoisson;

	//We simulate points uniformly in [-r/2, r/2]^2
	double u1 = *edge * runif(-0.5, 0.5),
	  u2 = *edge * runif(-0.5, 0.5);
      
	nKO = *nSites * *nSites;
	for (int j=0;j<*nSites;j++){
	  for (int k=0;k<*nSites;k++){
	    /* This is the bivariate normal density with 0 mean and
	       cov. matrix [cov11, cov12; cov12, cov22] */
	    double y = exp((-*cov22 * (coord[j] - u1) * (coord[j] - u1) + 2 * *cov12 *
			    (coord[j] - u1) * (coord[*nSites + k] - u2) - *cov11 *
			    (coord[*nSites + k] - u2) * (coord[*nSites + k] - u2)) *
			   itwiceDet) * thresh;
	      
	    ans[j + k * *nSites + i * *nSites * *nSites] = 
	      fmax2(y, ans[j + k * *nSites + i * *nSites * *nSites]);
	    
	    nKO -= (thresh <= ans[j + k * *nSites + i * *nSites * *nSites]);
	  }
	}
      }
    }
  }

  else{
    //Simulation part if a grid isn't used
    for (int i=0;i<*nObs;i++){
      double poisson = 0;
      int nKO = *nSites;
    
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */

	poisson += exp_rand();
	double ipoisson = 1 / poisson, thresh = uBound * ipoisson;

	//We simulate points uniformly in [-r/2, r/2]^2
	double u1 = *edge * runif(-0.5, 0.5),
	  u2 = *edge * runif(-0.5, 0.5);
      
	nKO = *nSites;
	for (int j=0;j<*nSites;j++){
	  /* This is the bivariate normal density with 0 mean and
	     cov. matrix [cov11, cov12; cov12, cov22] */
	  double y = exp((-*cov22 * (coord[j] - u1) * (coord[j] - u1) + 2 * *cov12 *
			  (coord[j] - u1) * (coord[*nSites + j] - u2) - *cov11 *
			  (coord[*nSites + j] - u2) * (coord[*nSites + j] - u2)) *
			 itwiceDet) * thresh;
	  
	  ans[i + j * *nObs] = fmax2(y, ans[i + j * *nObs]);
	  
	  nKO -= (thresh <= ans[i + j * *nObs]);
	}
      }
    }
  }
 
  PutRNGstate();

  /* Lastly, we multiply by the Lebesgue measure of the dilated
    compact set */
  if (*grid){
    for (int i=0;i<(*nSites * *nSites * *nObs);i++)
      ans[i] *= lebesgue;
  }
  
  else{
    for (int i=0;i<(*nSites * *nObs);i++)
      ans[i] *= lebesgue;
  }

  return;
}
