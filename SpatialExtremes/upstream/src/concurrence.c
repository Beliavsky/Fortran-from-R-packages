#include "header.h"

void empiricalConcProb(double *data, int *nSite, int *nObs, int *blockSize,
		       int *nBlock, double *concProb){
  /* This function compputes the pairwise concurence probability
     assuming the data are in the max domain of attraction */

  double *dataBlock = malloc(*blockSize * *nSite * sizeof(double));
  // will store the data corresponding to a given block

  for (int k=0;k<*nBlock;k++){//k is the block number
    double blockMax[2];//will store the pairwise block maxima
    int currentPair = 0,
      blockMaxIdx[2];//will store the indices for the pairwise block maxima

    for (int i=0;i<*blockSize;i++)
      for (int j=0;j<*nSite;j++)
	dataBlock[i + j * *blockSize] = data[k * *blockSize + i + j * *nObs];

    for (int i=0;i<(*nSite-1);i++){
      // Compute the block maxima for site i (and its index) and block number k
      blockMax[0] = dataBlock[i * *blockSize];
      blockMaxIdx[0] = 0;

      for (int l=1;l<*blockSize;l++){
	if (blockMax[0] < dataBlock[l + i * *blockSize]){
	  blockMax[0] = dataBlock[l + i * *blockSize];
	  blockMaxIdx[0] = l;
	}
      }

      for (int j=i+1;j<*nSite;j++){
	// Compute the block maxima for site j (and its index) and block number k
	blockMax[1] = dataBlock[j * *blockSize];
	blockMaxIdx[1] = 0;

	for (int l=1;l<*blockSize;l++){
	  if (blockMax[1] < dataBlock[l + j * *blockSize]){
	    blockMax[1] = dataBlock[l + j * *blockSize];
	    blockMaxIdx[1] = l;
	  }
	}

	if (blockMaxIdx[0] == blockMaxIdx[1])
	  // Concurent extremes
	  concProb[currentPair]++;

	currentPair++;
      }
    }
  }

  for (int i=0;i<(*nSite * (*nSite - 1) / 2);i++)
    concProb[i] /= (double) *nBlock;

  free(dataBlock);
  return;
}

void concProbKendall(double *data, int *nSite, int *nObs, double *concProb,
		     double *jackKnife, int *computeStdErr){

  const int nPair = *nSite * (*nSite - 1) / 2;

  if (!*computeStdErr){
#pragma omp parallel for          
    for (int currentPair=0;currentPair<nPair;currentPair++){
      int i,j;
      getSiteIndex(currentPair, *nSite, &i, &j);

      concProb[currentPair]=0;
      int nEffObs = 0;//the number of effective contribution (because of possible missing values)
      for (int k=0;k<(*nObs-1);k++){

	if (ISNA(data[k + i * *nObs]) || ISNA(data[k + j * *nObs]))
	  continue;

	for (int l=k+1;l<*nObs;l++){

	  if (ISNA(data[l + i * *nObs]) || ISNA(data[l + j * *nObs]))
	    continue;

	  nEffObs++;
	  concProb[currentPair] += sign(data[k + i * *nObs] - data[l + i * *nObs]) *
	    sign(data[k + j * *nObs] - data[l + j * *nObs]);
	}
      }

      concProb[currentPair] *= nEffObs == 0 ? NA_REAL : 1.0 / ((double) nEffObs);
    }
  }

  else {
    //the number of effective observation for the jackknife estimates
    int *nEffObsJack = malloc(*nObs * sizeof(int));

#pragma omp parallel for       
    for (int currentPair=0;currentPair<nPair;currentPair++){
      int i,j;
      getSiteIndex(currentPair, *nSite, &i, &j);

      // Reinitialization for the new pair

      //the number of effective contribution (because of possible missing values)
      int nEffObs = 0;

      for (int k=0;k<*nObs;k++){
	nEffObsJack[k] = 0;


	if (ISNA(data[k + i * *nObs]) || ISNA(data[k + j * *nObs])){
	  // If the discarded observation is missing, jackknife
	  // estimate is meaningless...
	  jackKnife[k + currentPair * *nObs] = NA_REAL;
	  continue;
	}

	for (int l=0;l<*nObs;l++){
	  if ((k == l) || ISNA(data[l + i * *nObs]) || ISNA(data[l + j * *nObs]))
	    continue;

	  nEffObs++; nEffObsJack[k]++;

	  double tmp = sign(data[k + i * *nObs] - data[l + i * *nObs]) *
	    sign(data[k + j * *nObs] - data[l + j * *nObs]);

	  concProb[currentPair] += tmp;
	  jackKnife[k + currentPair * *nObs] += tmp;
	}
      }

      /*
	Rmk: To compute the jackknife estimates (for a given pair of
	stations), we will use the following formula

	tau_(-k) = 2 / ((n - 1) (n-2)) * (n (n-1) / 2 * tau - sum_{l=1}^n
	sign(data_l(s_1) - data_k(s_1)) * sign(data_l(s_2) -
	data_k(s_2))
      */

      for (int k = 0; k<*nObs; k++)
	jackKnife[k + currentPair * *nObs] = nEffObsJack[k] == 0 ? NA_REAL :
	  (concProb[currentPair] - 2.0 * jackKnife[k + currentPair * *nObs]) /
	  ((double) (nEffObs - 2 * nEffObsJack[k]));

      concProb[currentPair] *= nEffObs == 0 ? NA_REAL : 1.0 / ((double) nEffObs);
    }
    free(nEffObsJack);
  }

  return;
}


void empiricalBootConcProb(double *data, int *nSite, int *nObs, int *blockSize,
			   double *concProb){

  const double normCst = lchoose(*nObs, *blockSize);
  const int nPair = *nSite * (*nSite - 1) / 2;

#pragma omp parallel for
  for (int currentPair=0;currentPair<nPair;currentPair++){
    int i, j;
    getSiteIndex(currentPair, *nSite, &i, &j);

    // For each pair compute the estimator
    concProb[currentPair] = 0;

    for (int k=0;k<*nObs;k++){
      int d = 0;

      for (int l=0;l<*nObs;l++)
	d += (data[l + i * *nObs] < data[k + i * *nObs]) && (data[l + j * *nObs] < data[k + j * *nObs]);

      concProb[currentPair] += exp(lchoose(d, *blockSize - 1) - normCst);
    }
  }

  return;
}
