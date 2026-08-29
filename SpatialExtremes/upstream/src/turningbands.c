#include "header.h"
/* WARNING: There's a bias in the simulation of fbm!!! */

void tbm(int *nobs, int *nsite, int *dim, int *covmod, int *grid, 
	 double *coord, double *nugget, double *sill, double *range,
	 double *smooth, int *nlines, double *ans){
  
  double normConst = sqrt(*sill / *nlines), sdnugget = sqrt(*nugget);
  int i, j, k, l, m, neffSite = *nsite;
  double freq, eucProd, u1, u2, G, phase, irange = 1 / *range, u, v, w,
    inorm, angle, r, theta;
 
  //rescale the coordinates
  for (i=(*nsite * *dim);i--;)
    coord[i] = coord[i] * irange;

  if (*grid)
    neffSite = R_pow_di(neffSite, *dim);

  for (i=(*nobs * neffSite);i--;)
    ans[i] = 0;
    
  double *lines = (double *)R_alloc(3 * *nlines, sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);

  GetRNGstate();
  if (*grid){
    //coord defines a grid
    
    for (i=*nobs;i--;){
      //Random rotation of the lines
      u = unif_rand() - 0.5;
      v = unif_rand() - 0.5;
      w = unif_rand() - 0.5;
      angle = runif(0, M_2PI);

      inorm = 1 / sqrt(u * u + v * v + w * w);
    
      u *= inorm;
      v *= inorm;
      w *= inorm;
    
      rotation(lines, nlines, &u, &v, &w, &angle);

      //Turning bands part
      if (*dim == 2){
	double cl;
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd =  cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * neffSite + k  + l * *nsite] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	  break;
	case 2:
	  //Cauchy
	  for (j=*nlines;j--;){
	    freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
	      ans[i * neffSite + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * neffSite + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * neffSite + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * neffSite + k  + l * *nsite] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	  break;
	case 6:
	  //Fractional brownian motion
	  for (j=*nlines;j--;){
	    r = rgamma(1 - 0.5 * *smooth, 1) / rgamma(0.5 * *smooth, 1);
	    theta = sqrt((1 + r) / R_pow(r, 0.5 * *smooth + 1));
	    freq = M_2PI * r;
	    phase = M_2PI * unif_rand();

	    for (k=*nsite;k--;){
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * neffSite + k + l * *nsite] += theta * cos(freq * eucProd + phase);
	      }
	    }

	    if (r < 1e-3){
	      double cosPhase = cos(phase);
	      for (k=*nsite;k--;){
		cl = coord[k] * lines[*nlines + j];
		for (l=*nsite;l--;){
		  eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + l * *nsite] -= theta * cosPhase;
		}
	      }
	      }
	  }
	  break;
	}
      }

      else if (*dim == 3){
	double cl1, cl2;
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
		}
	      }
	    }
	  }
	  break;
	  case 6:
	  //Fractional brownian motion
	  for (j=*nlines;j--;){
	    r = rgamma(1 - 0.5 * *smooth, 1) / rgamma(0.5 * *smooth, 1);
	    theta = sqrt((1 + r) / R_pow(r, 0.5 * *smooth + 1));
	    freq = M_2PI * r;
	    phase = M_2PI * unif_rand();

	    for (k=*nsite;k--;){
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * neffSite + k + *nsite * (l + m * *nsite)] += theta * cos(freq * eucProd + phase);
		}
	      }
	    }

	    if (r < 1e-3){
	      double cosPhase = cos(phase);
	      for (k=*nsite;k--;){
		cl1 = coord[k] * lines[*nlines + j];
		for (l=*nsite;l--;){
		  cl2 = coord[*nsite + l] * lines[*nlines + j];
		  for (m=*nsite;m--;){
		    eucProd = cl1 + cl2 + coord[2 * *nsite + l] * lines[2 * *nlines + j];
		    ans[i * neffSite + k + *nsite * (l + m * *nsite)] -= theta * cosPhase;
		  }
		}
	      }
	    }
	  }
	  break;
	}
      }
    }  
  }

  else{
    //coord doesn't define a grid
    for (i=*nobs;i--;){
      //Random rotation of the lines
      u = unif_rand() - 0.5;
      v = unif_rand() - 0.5;
      w = unif_rand() - 0.5;
      angle = runif(0, M_2PI);

      inorm = 1 / sqrt(u * u + v * v + w * w);
    
      u *= inorm;
      v *= inorm;
      w *= inorm;
    
      rotation(lines, nlines, &u, &v, &w, &angle);

      //Turning bands part
      if (*dim == 2){
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
	    }
	  }
	  break;
	  case 6:
	  //Fractional brownian motion
	  for (j=*nlines;j--;){
	    r = rgamma(1 - 0.5 * *smooth, 1) / rgamma(0.5 * *smooth, 1);
	    theta = sqrt((1 + r) / R_pow(r, 0.5 * *smooth + 1));
	    freq = M_2PI * r;
	    phase = M_2PI * unif_rand();

	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += theta * cos(freq * eucProd + phase);
	      }

	    if (r < 1e-3){
	      double cosPhase = cos(phase);
	      for (k=*nsite;k--;){
		eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
		ans[i + k * *nobs] -= theta * cosPhase;
	      }
	    }
	  }
	  break;
	}
      }

      else if (*dim == 3){
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
		coord[2 * *nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	
	    freq = sqrt(2 * rchisq(3)* G);
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
		coord[2 * *nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
	    }
	  }
	  break;
	case 6:
	  //Fractional brownian motion
	  for (j=*nlines;j--;){
	    r = rgamma(1 - 0.5 * *smooth, 1) / rgamma(0.5 * *smooth, 1);
	    theta = sqrt((1 + r) / R_pow(r, 0.5 * *smooth + 1));
	    freq = M_2PI * r;
	    phase = M_2PI * unif_rand();

	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
		coord[2 * *nsite + k] * lines[2 * *nlines + j];
	      ans[i * k * *nobs] += theta * cos(freq * eucProd + phase);
	      }

	    if (r < 1e-3){
	      double cosPhase = cos(phase);
	      for (k=*nsite;k--;){
		eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
		  coord[2 * *nsite + k] * lines[2 * *nlines + j];
		ans[i + k * *nobs] -= theta * cosPhase;
	      }
	    }
	  }
	  break;
	}
      }
    }  
  }
  
  if (*covmod != 6)
    normConst *= M_SQRT2;

  else
    normConst *= sqrt(4 * gammafn(0.5 * *smooth + 1) * 
		      gammafn(0.5 * (*dim + *smooth)) /
		      (R_pow(M_PI, *smooth) * gammafn(0.5 * *dim)));

  for (i=(neffSite * *nobs);i--;)
    ans[i] *= normConst;

  if (*nugget != 0)
    for (i=(neffSite * *nobs);i--;)
      ans[i] += sdnugget * norm_rand();
  
  PutRNGstate();
  return;
}
