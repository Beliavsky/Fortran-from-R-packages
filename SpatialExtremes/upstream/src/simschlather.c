#include "header.h"

void rschlathertbm(double *coord, int *nObs, int *nSite, int *dim,
		   int *covmod, int *grid, double *nugget, double *range,
		   double *smooth, double *uBound, int *nlines,
		   double *ans){
  /* This function generates random fields from the Schlather model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      nugget: the nugget parameter
     range: the range parameter
    smooth: the smooth parameter
    uBound: the uniform upper bound for the stoch. proc.
    nlines: the number of lines used for the TBM algo
       ans: the generated random field */

  /* lagi, lagj are integers useful to fill in the output depending on
     if locations are on a grid or not */
  int neffSite, lagi = 1, lagj = 1;
  double sill = 1 - *nugget;
  const double irange = 1 / *range;

  //rescale the coordinates
  for (int i= 0; i<(*nSite * *dim);i++)    
    coord[i] = coord[i] * irange;

  double *lines = malloc(3 * *nlines * sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  GetRNGstate();

  double *gp = malloc(neffSite * sizeof(double));

  for (int i = 0; i < *nObs;i++){
    int nKO = neffSite;
    double poisson = 0;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
            
      /* ------- Random rotation of the lines ----------*/
      double u = unif_rand() - 0.5,
	v = unif_rand() - 0.5,
	w = unif_rand() - 0.5,
	angle = runif(0, M_2PI),	
	inorm = 1 / sqrt(u * u + v * v + w * w);
      
      u *= inorm;
      v *= inorm;
      w *= inorm;
      
      rotation(lines, nlines, &u, &v, &w, &angle);
      /* -------------- end of rotation ---------------*/
      
      poisson += exp_rand();
      double ipoisson = 1 / poisson,
	thresh = *uBound * ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      for (int j = 0; j < neffSite; j++)
	gp[j] = 0;
      tbmcore(nSite, &neffSite, dim, covmod, grid, coord, nugget,
	      &sill, range, smooth, nlines, lines, gp);
      
      nKO = neffSite;
      for (int j = 0; j < neffSite; j++){
	ans[j * lagj + i * lagi] = fmax2(gp[j] * ipoisson, ans[j * lagj + i * lagi]);
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
	
      }
    }
  }
  
  PutRNGstate();

  //Lastly we multiply by 1 / E[max(0, Y)]
  const double imean = M_SQRT2 * M_SQRT_PI;
  for (int i=0; i < (neffSite * *nObs); i++)    
    ans[i] *= imean;

  free(lines); free(gp);
  return;
}

void rschlatherdirect(double *coord, int *nObs, int *nSite, int *dim,
		      int *covmod, int *grid, double *nugget, double *range,
		      double *smooth, double *uBound, double *ans){
  /* This function generates random fields for the Schlather model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      nugget: the nugget parameter
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);
  
  /* Compute the Cholesky decomposition of the covariance matrix */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info FCONE);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
 
  for (int i = 0; i < *nObs; i++){
    double poisson = 0;
    int nKO = neffSite;
        
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      poisson += exp_rand();
      double ipoisson = 1 / poisson, thresh = *uBound * ipoisson;
            
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      for (int j=0;j<neffSite;j++)
	gp[j] = norm_rand();
      
      F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt
		      FCONE FCONE FCONE);

      nKO = neffSite;
      for (int j = 0; j < neffSite; j++){
	ans[j * lagj + i * lagi] = fmax2(gp[j] * ipoisson, ans[j * lagj + i * lagi]);
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
      }
    }
  }

  PutRNGstate();
  //Lastly we multiply by 1 / E[max(0, Y)]
  const double imean = M_SQRT2 * M_SQRT_PI;
  for (int i = 0; i < (neffSite * *nObs); i++)    
    ans[i] *= imean;

  free(covmat); free(gp);
  return;
}

void rschlathercirc(int *nObs, int *ngrid, double *steps, int *dim,
		    int *covmod, double *nugget, double *range,
		    double *smooth, double *uBound, double *ans){
  /* This function generates random fields from the Schlather model

     nObs: the number of observations to be generated
    ngrid: the number of locations along one axis
      dim: the random field is generated in R^dim
   covmod: the covariance model
     nugget: the nugget parameter
    range: the range parameter
   smooth: the smooth parameter
   uBound: the uniform upper bound for the stoch. proc.
      ans: the generated random field */

  int i, j, k = -1, nbar = R_pow_di(*ngrid, *dim), r, m;
  const double sill = 1 - *nugget, zero = 0;
  double *rho, *irho;
    //Below is a table of highly composite numbers
  int HCN[39] = {1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240,
		 360, 720, 840, 1260, 1680, 2520, 5040, 7560,
		 10080, 15120, 20160, 25200, 27720, 45360, 50400,
		 55440, 83160, 110880, 166320, 221760, 277200,
		 332640, 498960, 554400, 665280, 720720, 1081080};

    
  /* Find the smallest size m for the circulant embedding matrix */
  {
    int dummy = 2 * (*ngrid - 1);
    do {
      k++;
      m = HCN[k];
    } while (m < dummy);
  }
  
  /* ---------- beginning of the embedding stage ---------- */
  int mbar = m * m, halfM = m / 2, notPosDef = 0;
  do {
    double *dist;
    dist = (double *)R_alloc(mbar, sizeof(double));

    notPosDef = 0;
    //Computation of the distance
    for (r=mbar;r--;){
      i = r % m;
      j = r / m;
      
      if (i > halfM)
	i -= m;
      
      if (j > halfM)
	j -= m;
      
      dist[r] = hypot(steps[0] * i, steps[1] * j);
    }

    //Computations of the covariances
    rho = (double *)R_alloc(mbar, sizeof(double));
    irho = (double *)R_alloc(mbar, sizeof(double));
    for (i=mbar;i--;)
      irho[i] = 0;

    switch (*covmod){
    case 1:
      whittleMatern(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 2:
      cauchy(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 3:
      powerExp(dist, mbar, zero, sill, *range, *smooth, rho);
      break;
    case 4:
      bessel(dist, mbar, *dim, zero, sill, *range, *smooth, rho);
      break;
    }

    /* Compute the eigen values to check if the circulant embbeding
       matrix is positive definite */

    /* Note : The next lines is only valid for 2d random fields. I
       need to change if there are m_1 \neq m_2 as I suppose that m_1
       = m_2 = m */
    int maxf, maxp, *iwork;
    double *work;

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, m, m, 1, -1, work, iwork);

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, 1, m, m, -1, work, iwork);

    //Check if the eigenvalues are all positive
    for (i=mbar;i--;){
      notPosDef |= (rho[i] <= 0) || (fabs(irho[i]) > 0.001);
    }

    if (notPosDef){
      k++;
      m = HCN[k];
      halfM = m / 2;
      mbar = m * m;
    }

    if (k > 30)
      error("Impossible to embbed the covariance matrix");
    
  } while (notPosDef);
  /* --------- end of the embedding stage --------- */

  /* Computation of the square root of the eigenvalues */
  for (i=mbar;i--;){
    rho[i] = sqrt(rho[i]);
    irho[i] = 0;//No imaginary part
  }

  int mdag = m / 2 + 1, mdagbar = mdag * mdag;
  double *a, *ia, isqrtMbar = 1 / sqrt(mbar);

  a = (double *)R_alloc(mbar, sizeof(double));
  ia = (double *)R_alloc(mbar, sizeof(double));
  
  GetRNGstate();
  for (i=*nObs;i--;){
    int nKO = nbar;
    double poisson = 0;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      int j;
      double ipoisson, thresh, *gp;
      
      gp = (double *)R_alloc(nbar, sizeof(double));
      
      poisson += exp_rand();
      ipoisson = 1 / poisson;
      thresh = *uBound * ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      circcore(rho, a, ia, m, halfM, mdag, mdagbar, *ngrid, nbar, isqrtMbar, *nugget, gp);
      
      nKO = nbar;
      for (j=nbar;j--;){
	ans[j + i * nbar] = fmax2(gp[j] * ipoisson, ans[j + i * nbar]);
	nKO -= (thresh <= ans[j + i * nbar]);
	
      }
    }
  }
  
  PutRNGstate();
  
  //Lastly we multiply by 1 / E[max(0, Y)]
  const double imean = M_SQRT2 * M_SQRT_PI;
  for (i=(nbar * *nObs);i--;)    
    ans[i] *= imean;
  
  return;
}

void tbmcore(int *nsite, int *neffSite, int *dim, int *covmod,
	     int *grid, double *coord, double *nugget, double *sill,
	     double *range, double *smooth, int *nlines, double *lines,
	     double *ans){
  /* This function is the same as the tbm function except that it
     generates only one realisation of the random field. This is only
     useful for the call by rschlathertbm - for CPU reasons. */

  int j;
  
  if (*grid){
    //coord defines a grid
    int k, l, m;
    double freq, eucProd, u1, u2, G, phase;
    double cl, lj, cl1, cl2;

    switch (*dim){
    case 2:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd =  cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -2 / *smooth) *
		   R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			 (2 - *smooth) / *smooth));
	  
	  freq = sqrt(2 * rchisq(3) * G);
	  phase = M_2PI * unif_rand();
	
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      }
      break;
      
    case 3:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth ,1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -2 / *smooth) *
		   R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			 (2 - *smooth) / *smooth));
	
	  freq = sqrt(2 * rchisq(3) * G);
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5));
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}    
      }
    }
  }  

  else{
    //coord doesn't define a grid
    int k;
    double freq, eucProd, u1, u2, G, phase;
    switch (*dim){
    case 2:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 2:
	 //Cauchy
	 for (j=*nlines;j--;){
	   freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 3:
	 //Powered exponential
	 for (j=*nlines;j--;){
	   u1 = rexp(1);
	   u2 = runif(-M_PI_2, M_PI_2);
	   G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -2 / *smooth) *
		    R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			  (2 - *smooth) / *smooth));
	
	   freq = sqrt(2 * rchisq(3) * G);
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 4:
	 //Bessel
	 for (j=*nlines;j--;){
	   freq = sqrt(beta(1.5, *smooth - 0.5));
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 5:
	 //Gaussian
	 for (j=*nlines;j--;){
	   freq = sqrt(2 * rchisq(3));
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
      }
      break;    
    case 3:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth ,1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -2 / *smooth) *
		   R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			 (2 - *smooth) / *smooth));
	  
	  freq = sqrt(2 * rchisq(3) * G);
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}    
      }
    }
  }  

  const double normConst = sqrt(2 * *sill / *nlines);
  for (j=*neffSite;j--;)    
    ans[j] *= normConst;

  if (*nugget != 0){
    const double sdnugget = sqrt(*nugget);
    for (j=*neffSite;j--;)
      ans[j] += sdnugget * norm_rand();
  }
  
  return;
}

void circcore(double *rho, double *a, double *ia, int m, int halfM, int mdag,
	      int mdagbar, int ngrid, int nbar, double isqrtMbar, double nugget,
	      double *ans){
  /* This function is the same as the circemb function except that it
     generates only one realisation of the random field. This is only
     useful for the call by rschlathercirc - for CPU reasons. */
  for (int r=0;r<mdagbar;r++){
    /* Below is the procedure 5.2.4 in Wood and Chan */

    //Computation of the cardinality of A(j)
    int j1, j2, i = r % mdag, j = r / mdag;
    double u, v;

    int card = (i != 0) * (i != halfM) + 2 * (j != 0) * (j != halfM);
      
    switch (card){
    case 3:
      //B(1) = {1}, B^c(1) = {2}
      j1 = (m - i) + m * j;
      j2 = i + m * (m - j);
      u = norm_rand();
      v = norm_rand();
      a[j1] = ia[j1] = M_SQRT1_2 * rho[j1];
      a[j1] *= u;
      ia[j1] *= v;
      a[j2] = ia[j2] = M_SQRT1_2 * rho[j2];
      a[j2] *= u;
      ia[j2] *= -v;
	
      //B(2) = {1,2}, B^c(2) = {0}
      j1 = (m - i) + m * (m - j);
      j2 = i + m * j;
      u = norm_rand();
      v = norm_rand();
      a[j1] = ia[j1] = M_SQRT1_2 * rho[j1];
      a[j1]*= u;
      ia[j1] *= v;
      a[j2] = ia[j2] = M_SQRT1_2 * rho[j2];
      a[j2]*= u;
      ia[j2] *= -v;      
      break;
    case 1:
      //B(1) = 0, B^c(1) = {1}
      j1 = i + m * j;
      j2 = m - i + m * j;
      u = norm_rand();
      v = norm_rand();
      a[j1] = ia[j1] = M_SQRT1_2 * rho[j1];
      a[j1] *= u;
      ia[j1] *= v;
      a[j2] = ia[j2] = M_SQRT1_2 * rho[j2];
      a[j2] *= u;
      ia[j2] *= -v;
      break;
    case 2:
      //B(1) = 0, B^c(1) = {2}
      j1 = i + m * j;
      j2 = i + m * (m - j);
      u = norm_rand();
      v = norm_rand();
      a[j1] = ia[j1] = M_SQRT1_2 * rho[j1];
      a[j1] *= u;
      ia[j1] *= v;
      a[j2] = ia[j2] = M_SQRT1_2 * rho[j2];
      a[j2] *= u;
      ia[j2] *= -v;
      break;
    case 0:
      j1 = i + m * j;
      a[j1] = rho[j1] * norm_rand();
      ia[j1] = 0;
      break;      
    }
  }

  /* ---------- Computation of Q \Lambda^1/2 Q* Z ------------ */
  int maxf, maxp, *iwork;
  double *work;
    
  /* The next lines is only valid for 2d random fields. I need to
     change if m_1 \neq m_2 as here I suppose that m_1 = m_2 = m */
  fft_factor(m, &maxf, &maxp);
  work = (double *)R_alloc(4 * maxf, sizeof(double));
  iwork = (int *)R_alloc(maxp, sizeof(int));
  fft_work(a, ia, m, m, 1, -1, work, iwork);
    
  fft_factor(m, &maxf, &maxp);
  work = (double *)R_alloc(4 * maxf, sizeof(double));
  iwork = (int *)R_alloc(maxp, sizeof(int));
  fft_work(a, ia, 1, m, m, -1, work, iwork);
      
  int i;
  for (i=nbar;i--;)
    ans[i] = isqrtMbar * a[i % ngrid + m * (i / ngrid)];
    
  if (nugget > 0){
    double sqrtNugget = sqrt(nugget);
    
    for (i=nbar;i--;)
      ans[i] += sqrtNugget * norm_rand();
  
  }

  return;
}
