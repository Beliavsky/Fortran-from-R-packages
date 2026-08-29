#include "header.h"

void rcondMaxLin(double *data, double *dsgnMat, int *p, int *nSite, int *nSim,
		 double *Z){
  /* This function performs conditional simulations from a max-linear
     unit Frechet process, well actually only the Zs */

  /* 1. Compute the hitting bounds, i.e. the upper bound for the unit
     Frechet r.v. */
  double *hitBounds = (double *) R_alloc(*p, sizeof(double));
  
  for (int j=*p;j--;){
    hitBounds[j] = R_PosInf;
   
    for (int i=*nSite;i--;)
      if (dsgnMat[i + j * *nSite] != 0)
	hitBounds[j] = fmin2(hitBounds[j], data[i] / dsgnMat[i + j * *nSite]);
  }

  /* 2. Compute the hitting matrix */
  int *hitMat = (int *) R_alloc(*p * *nSite, sizeof(int));
  for (int i=*nSite;i--;)
    for (int j=*p;j--;)
      hitMat[j * *nSite + i] = (fabs(dsgnMat[j * *nSite + i] * hitBounds[j] -
				     data[i]) <= 1e-14);

  /* 3. Partitioning {0, ...., nSite - 1} i.e. the sets I_s of Stoev
     and Wang */

  /* We first check if there is no connection between site --- this
     will arise in many situations*/
  int flag = 0;
  for (int j=*p;j--;){
    int dummy = 0;

    for (int i=*nSite;i--;)
      dummy += hitMat[i + j * *nSite];

    flag += (dummy == 1);
  }

  GetRNGstate();
  
  if (flag == *p){
    /* No connection i.e. I_s = {0}, {1}, ..., {nSite - 1}, r = nSite
       and J_s = Jbar_s */

    for (int i=*nSite;i--;){
      // We need to know the size of the set J_s for all s
      int sizeJs = 0;
    
      for (int j=*p;j--;)
	sizeJs += hitMat[i + j * *nSite];

      if (sizeJs){
	/* i.e. J_s \neq \empyset. Should not arise in theory but
	   numerically do! */
	int *Js = (int *) R_alloc(sizeJs, sizeof(int));
	double *weights = (double *) R_alloc(sizeJs, sizeof(double));

	int current = 0;
	double dummy = 1.0;
	for (int j=*p;j--;){
	  if (hitMat[i + j * *nSite]){
	    Js[current] = j;
	    dummy *= exp(-1 / hitBounds[j]);
	    current++;
	  }
	}
	
	for (int j=sizeJs;j--;)
	  weights[j] = dummy / hitBounds[Js[j]];
	
	dummy = 0.0;
	for (int j=sizeJs;j--;)
	  dummy += weights[j];
	
	weights[0] /= dummy;
	
	for (int j=1;j<sizeJs;j++)
	  weights[j] = weights[j-1] + weights[j] / dummy;
	
	for (int iter=*nSim;iter--;){
	  for (int j=sizeJs;j--;)
	    Z[iter * *p + Js[j]] = 1 / (1 / hitBounds[Js[j]] - log(unif_rand()));
	
	  double u = unif_rand();
	  int idxMixt = 0;
	  for (int j=0;j<sizeJs;j++)
	    idxMixt += u > weights[j];
	  
	  Z[iter * *p + Js[idxMixt]] = hitBounds[Js[idxMixt]];
	}
      }
    }
  }

  else {
    /* Some rows do cover, i.e., I_s \neq {0}, {1}, ..., {nSite -
       1} */
    
    // Initializing I_s
    int *Is = (int *) R_alloc(*nSite, sizeof(int));
    for (int i=*nSite;i--;)
      Is[i] = i;

    /* I_s will be stored in the following way. We reorder {1, ..., n}
       and stack the cardinal of each I_s into cardIs. cumCardIs is
       the sum of the cardIs */
    int remains = *nSite, cardIs = 0, cumCardIs = 0;
    
    while (remains){
      cumCardIs++;
      cardIs = 1;
      remains--;

      int current = Is[cumCardIs - 1];
      int *unionCol = (int *) R_alloc(*p, sizeof(int));
      /* unionCol corresponds to the sum of the rows that cover each
	 other i.e. which belongs to the same equivalence class I_s. It
	 allows us to "easily" know if a candidate row will belong or
	 not to the equivalence class */
      
      for (int j=*p;j--;)
	unionCol[j] = hitMat[current + j * *nSite];

      /* endClass is a logical indicating wether the end of an
	 equivalence class occurred */
      int endClass = remains;

      
      while (endClass){

	endClass = remains;

	for (int i=cumCardIs;i<*nSite;i++){
	  int doCover = 0;
	  
	  /*The i-th row will belong to the current equivalence class
	    iff <unionCol, hitMat[i,]> is positive i.e., doCover >
	    0 */

	  for (int j=*p;j--;)
	    doCover += unionCol[j] * hitMat[Is[i] + j * *nSite];
	  
	  if (doCover){
	    for (int j=*p;j--;)
	      unionCol[j] += hitMat[Is[i] + j * *nSite];
	    
	    Is[cumCardIs] = current = Is[i];
	    Is[i] = cumCardIs;
	    cumCardIs++;
	    cardIs++;
	    remains--;
	  }
	  
	  else
	    endClass--;
	}
      }

      // Now we need to compute the sets J_s and Jbar_s
      int *dummy1 = (int *) R_alloc(*p, sizeof(int)),
	*dummy2 = (int *) R_alloc(*p, sizeof(int));

      for (int i=0;i<*p;i++){
	dummy1[i] = 1;
	dummy2[i] = 0;
      }

      for (int i=cardIs;i--;){
	for (int j=*p;j--;){
	  dummy1[j] *= hitMat[Is[cumCardIs - i - 1] + j * *nSite];
	  dummy2[j] += hitMat[Is[cumCardIs - i - 1] + j * *nSite];
	}
      }

      int cardJs = 0, cardJbars = 0;
      for (int i=*p;i--;){
	cardJs += dummy1[i];
	cardJbars += (dummy2[i] > 0);
      }      

      if (cardJs){
	int *Js = (int *) R_alloc(cardJs, sizeof(int)),
	  *Jbars = (int *) R_alloc(cardJbars, sizeof(int));

	int current1 = 0, current2 = 0;
	double dummy = 1;
	for (int i=*p;i--;){
	  if (dummy1[i]){
	    Js[current1] = i;
	    current1++;
	  }
	  
	  if (dummy2[i]){
	    Jbars[current2] = i;
	    dummy *= exp(-1 / hitBounds[i]);
	    current2++;
	  }
	}

	double *weights = (double *) R_alloc(cardJs, sizeof(double));
	for (int j=cardJs;j--;)
	  weights[j] = dummy / hitBounds[Jbars[j]];
	
	dummy = 0.0;
	for (int j=1;j<cardJs;j++)
	  dummy += weights[j];
	
	weights[0] /= dummy;
	for (int j=1;j<cardJs;j++)
	  weights[j] = weights[j-1] + weights[j] / dummy;
	
	for (int iter=*nSim;iter--;){
	  for (int j=cardJbars;j--;)
	    Z[iter * *p + Jbars[j]] = 1 / (1 / hitBounds[Jbars[j]] - log(unif_rand()));
	
	  double u = unif_rand();
	  int idxMixt = 0;
	  for (int j=0;j<cardJs;j++)
	    idxMixt += u > weights[j];
	  
	  Z[iter * *p + Js[idxMixt]] = hitBounds[Js[idxMixt]];
	}
      }
    }
  }

  PutRNGstate();
  return;
}


void maxLinear(int *nSim, double *dsgnMat, double *Z, int *nSite, int *p,
	       int *grid, double *sim){
  /* This function generates realisation from a (unit Frechet)
     max-linear model */

  if (*grid){
    for (int i=*nSim;i--;){
      for (int j=*nSite;j--;){
	sim[j + i * *nSite] = R_NegInf;
	
	for (int k=*p;k--;)
	  if (dsgnMat[j + *nSite * k] != 0)
	    sim[j + i * *nSite] = fmax2(sim[j + i * *nSite],
					dsgnMat[j + *nSite * k] * Z[i * *p + k]);
      }
    }
  }

  else{
    for (int i=*nSim;i--;){
      for (int j=*nSite;j--;){
	sim[i + j * *nSim] = R_NegInf;
	
	for (int k=*p;k--;)
	  if (dsgnMat[j + *nSite * k] != 0)
	    sim[i + j * *nSim] = fmax2(sim[i + j * *nSim],
				       dsgnMat[j + *nSite * k] * Z[i * *p + k]);
      }
    }
  }

  return;
}
  
  
void maxLinDsgnMat(double *coord, double *grid, int *nSite, int *p,
		   double *areaPixel, int *dim, double *param,
		   double *dsgnMat){
  /* This function computes the design matrix for a max-linear
     model */

  if (*dim == 1){
    // Unidimensional discretized Smith process

    double iVar = 1 / param[0],
      cst = *areaPixel * M_1_SQRT_2PI * sqrt(iVar);

    for (int i=*nSite;i--;)
      for (int j=*p;j--;)
	dsgnMat[i + j * *nSite] = exp(- 0.5 * (coord[i] - grid[j]) *
				      (coord[i] - grid[j]) * iVar) * cst;
  }

  else if (*dim == 2){
    // Two dimensional discretized Smith process
    double idet = 1 / (param[0] * param[2] - param[1] * param[1]),
      cst = *areaPixel / M_2PI * sqrt(idet);

    for (int i=*nSite;i--;)
      for (int j=*p;j--;){
	double dummy1 = coord[i] - grid[j],
	  dummy2 = coord[i + *nSite] - grid[j + *p];

	dsgnMat[i + j * *nSite] = exp(-0.5 * (param[2] * dummy1 * dummy1 -
					      2 * param[1] * dummy1 * dummy2 +
					      param[0] * dummy2 * dummy2) *
				      idet) * cst;
      }

  }

  else
    error("not implemented yet!");

  for (int i=(*p * *nSite);i--;){
    //Set to zero all entries <= 1e-8
    if (dsgnMat[i] <= 1e-8)
      dsgnMat[i] = 0;
  }

  return;  
}
