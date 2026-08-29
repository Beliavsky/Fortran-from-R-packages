#include "header.h"

void circemb(int *nsim, int *ngrid, double *steps, int *dim, int *covmod,
	     double *nugget, double *sill, double *range, double *smooth,
	     double *ans){

  int i, j, k = -1, r, nbar = *ngrid * *ngrid, m;
  //irho is the imaginary part of the covariance -> 0
  double *rho, *irho;
  const double zero = 0;
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
    double *dist = malloc(mbar * sizeof(double));

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
      whittleMatern(dist, mbar, zero, *sill, *range, *smooth, rho);
      break;
    case 2:
      cauchy(dist, mbar, zero, *sill, *range, *smooth, rho);
      break;
    case 3:
      powerExp(dist, mbar, zero, *sill, *range, *smooth, rho);
      break;
    case 4:
      bessel(dist, mbar, *dim, zero, *sill, *range, *smooth, rho);
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
      notPosDef |= (rho[i] < 0) || (fabs(irho[i]) > 0.001);
    }

    if (notPosDef){
      k++;
      m = HCN[k];
      halfM = m / 2;
      mbar = m * m;
    }

    if (k > 30)
      error("Impossible to embbed the covariance matrix");

    free(dist);
    
  } while (notPosDef);
  /* --------- end of the embedding stage --------- */

  /* Computation of the square root of the eigenvalues */
  for (i=mbar;i--;){
    rho[i] = sqrt(rho[i]);
    irho[i] = 0;//No imaginary part
  }

  int mdag = m / 2 + 1, mdagbar = mdag * mdag;
  double isqrtMbar = 1 / sqrt(mbar);

  double *a = malloc(mbar * sizeof(double)),
    *ia = malloc(mbar * sizeof(double));
    
  GetRNGstate();
  for (k=*nsim;k--;){
    
    /* ---------- Simulation from \Lambda^1/2 Q* Z ------------ */
    for (r=mdagbar;r--;){
      /* Below is the procedure 5.2.4 in Wood and Chan */

      //Computation of the cardinality of A(j)
      int j1, j2,i = r % mdag, j = r / mdag;
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
        
    for (i=nbar;i--;)
      ans[i + k * nbar] = isqrtMbar * a[i % *ngrid + m * (i / *ngrid)];
  }
  PutRNGstate();  
    
  if (*nugget > 0){
    int dummy = *nsim * nbar;
    double sqrtNugget = sqrt(*nugget);
    
    GetRNGstate();
    for (i=dummy;i--;)
      ans[i] += sqrtNugget * norm_rand();

    PutRNGstate();
  }

  free(a); free(ia);

  return;
}
