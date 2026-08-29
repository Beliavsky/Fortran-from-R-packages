#include "header.h"

void vandercorput(int *n, double *coord){
  //This function computes the Van der Corput sequence in R^3
  int i, k, r;
  double base, u, v;
  
  for (i=*n;i--;){
    //Binary decomposition
    k = i;
    u = 0;
    base = 2.;

    while (k){
      r = k % 2;
      u += r / base;
      base *= 2;
      k = k / 2;
    }

    //Ternary decomposition
    k = i;
    v = 0;
    base = 3;

    while (k){
      r = k % 3;
      v += r / base;
      base *= 3;
      k = k / 3;
    }

    coord[i] = cos(M_2PI * u) * sqrt(1 - v * v);
    coord[*n + i] = sin(M_2PI * u) * sqrt(1 - v * v);
    coord[2 * *n + i] = v;
  }

  return;
}

void rotation(double *coord, int *n, double *u, double *v,
	      double *w, double *angle){
  /* This function performs a rotation of coord w.r.t. to a rotation
     axis (u, v, w) and a rotation angle */
  
  int i;
  double r, px, py, pz, bx, by, bz, cx, cy, cz,
    cosAngle = cos(*angle), sinAngle = sin(*angle);


  for (i=*n;i--;){
    //Projection of the points onto the rotation axis
    r = coord[i] * *u + coord[*n + i] * *v + coord[2 * *n + i] * *w;
    
    px = r * *u;
    py = r * *v;
    pz = r * *w;

    r = sqrt((coord[i] - px) * (coord[i] - px) + 
	     (coord[*n + i] - py) * (coord[*n + i] - py) +
	     (coord[2 * *n + i] - pz) * (coord[2 * *n + i] - pz));

    bx = (coord[i] - px) / r;
    by = (coord[*n + i] - py) / r;
    bz = (coord[2 * *n + i] - pz) / r;
    
    cx = *v * bz - *w * by;
    cy = *w * bx - *u * bz;
    cz = *u * by - *v * bx;

    coord[i] = px + r * cosAngle * bx + r * sinAngle * cx;
    coord[*n + i] = py + r * cosAngle * by + r * sinAngle * cy;
    coord[2 * *n + i] = pz + r * cosAngle * bz + r * sinAngle * cz;
  }

  return;
}
