/*
    kza.c
    Copyright (C) 2005  Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <errno.h>
#include <R.h>
#include <Rdefines.h>
#include "kz.h"

typedef struct {
    int x,y;
} point;

/* upper left is 0,0 */
typedef struct{
    point upperLeft, lowerRight;
} boundingBox ;

/* adaptive curvature function */
static double adaptive(double d, double m)
{
	/* m = max of the difference metric; a constant input gives m == 0 and
	   d/m = 0/0 = NaN, and (int)floor(q*NaN) downstream is undefined
	   behavior (it only survived via the min-window clamp by accident).
	   m <= 0 means no structure anywhere: keep the full window. */
	if (m <= 0.0) return 1.0;
	return( 1 - (d/m) );
}

static int cmp_dbl(const void *a, const void *b)
{
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/* Normalizer for the adaptive shrink factor: the literal maximum of the
   difference metric when prob >= 1 (the published behavior), otherwise a
   type-7 quantile of its finite values -- so a handful of extreme pixels
   don't single-handedly set the scale for the whole image -- with a
   fallback to the maximum when the quantile collapses to zero: structure
   occupying less than (1 - prob) of an otherwise flat field would
   otherwise read as "no structure anywhere" and be smoothed away. */
static double normalizer(SEXP d, double prob)
{
    long i, nf, n = LENGTH(d);
    double h, q, *buf;
    long lo;

    if (!(prob > 0.0 && prob < 1.0)) return R_maximum(d);
    buf = (double *) R_alloc(n, sizeof(double));
    for (i = 0, nf = 0; i < n; i++)
        if (R_FINITE(REAL(d)[i])) buf[nf++] = REAL(d)[i];
    if (nf == 0) return R_maximum(d);
    qsort(buf, nf, sizeof(double), cmp_dbl);
    h = (nf - 1) * prob;
    lo = (long) floor(h);
    q = buf[lo] + (h - lo) * (buf[lo + 1 < nf ? lo + 1 : lo] - buf[lo]);
    if (q <= 0.0) return R_maximum(d);
    return q;
}

static double mavg1d(SEXP x, int a, int b)
{
	double s=0.00;
	long i, z;
  
	for(i=a, z=0; i<b; i++)
	{
		if (R_FINITE(REAL(x)[i])) {
			z++;
			s += REAL(x)[i];
		}
	}
	if (z==0) return R_NaN;
	return s/z;
}

SEXP kza1d(SEXP v, SEXP y, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob)
{
	int i;
	int n,q;
	SEXP d, dprime;
	double m;
	long t, qh, qt;
	int min_window_length;
	SEXP tmp, ans;
	double eps;

    eps = REAL(tolerance)[0];
	q = INTEGER_VALUE(window);
	min_window_length = INTEGER_VALUE(minimum_window_length);
    n = LENGTH(y);
    
	PROTECT(d = allocVector(REALSXP, n));
	PROTECT(dprime = allocVector(REALSXP, n));
    R_differenced(y, d, dprime, q);
 
    m = normalizer(d, REAL(normalize_prob)[0]);

    PROTECT(tmp = allocVector(REALSXP, n));
    copyVector(tmp, v);
    
    PROTECT(ans = allocVector(REALSXP, n));
    for(i=0; i<INTEGER_VALUE(iterations); i++) {
    	for (t=0; t<n; t++) {
		    if (fabs(REAL(dprime)[t]) < eps) { /* dprime[t] = 0 */
			    qh = (int) floor(q*adaptive(REAL(d)[t], m));
			    qt = (int) floor(q*adaptive(REAL(d)[t], m));
		    } else if (REAL(dprime)[t] < 0) {
		    	qh = q;
    		    qt = (int) floor(q*adaptive(REAL(d)[t], m));
	    	} else {
		    	qh = (int) floor(q*adaptive(REAL(d)[t], m));
			    qt = q;
			}
			qt = ((qt) < min_window_length) ? min_window_length : qt;
			qh = ((qh) < min_window_length) ? min_window_length : qh;
	
	        /* check bounds */
        	qh = (qh > n-t-1) ? n-t-1 : qh; /* head past end of series */
            qt = (qt > t) ? t : qt;  	        		
   		    REAL(ans)[t] = mavg1d(tmp, t-qt, t+qh+1); 
	    }
	    copyVector(tmp, ans);
	}
	UNPROTECT(4);
	return ans;
}

static double mavg2d(SEXP v, boundingBox b)
{
	double s;
	int i, j, z;
	int nr;

    nr = nrows(v);
    s=0.0;
    /* i spans the box's x (column) range and j its y (row) range, so the
       column-major element is [row + nr*col] = [j + nr*i]. The shipped
       code indexed [i + nr*j] -- transposed: self-consistently wrong on
       square matrices and out-of-bounds reads on non-square ones. */
	for(i=b.upperLeft.x, z=0; i<=b.lowerRight.x; i++) {
	    for(j=b.upperLeft.y; j<=b.lowerRight.y; j++) {
		    if (R_FINITE(REAL(v)[j+nr*i])) {
			    z++;
			    s += REAL(v)[j+nr*i];
		    }
		}
	}
	if (z==0) return R_NaN;
	return s/z;
}

SEXP kza2d(SEXP v, SEXP kz, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob)
{
    int min_window_length;
    int q1, q2;
	int i, j, k;
	SEXP dx, dy, dprimex, dprimey;
	double mx, my;
	int qx_head, qx_tail, qy_head, qy_tail;
	int col_head, col_tail, row_head, row_tail;
	SEXP dim, tmp, ans;
	int nr, nc;
	boundingBox b;
	double epsilon;

    epsilon = REAL(tolerance)[0];
	if (length(window)<2) {q1 = q2 = INTEGER_VALUE(window);}
	else {q1 = INTEGER(window)[0]; q2 = INTEGER(window)[1];}
	
	min_window_length = INTEGER_VALUE(minimum_window_length);
	
    dim = GET_DIM(v);
	nr = INTEGER(dim)[0]; nc = INTEGER(dim)[1];
	
	PROTECT(tmp = allocMatrix(REALSXP, nr, nc));
    copyMatrix(tmp, v, 0);

	/* calculate d1 = |Z(x+q,y) - Z(x-q,y)| , d2 = |Z(x,y+q) - Z(x,y-q)| */
	PROTECT(dx = allocMatrix(REALSXP, nr, nc));
	PROTECT(dy = allocMatrix(REALSXP, nr, nc));

    /*  origin is upper left (0,0)
        positive directions are down (y) and right (x) 
     */
	for(i=0; i<nr; i++) { /* ith row (y) */
	    for(j=0;j<nc;j++) { /* jth column (x) */
	        row_tail = (i-q1<0 ? 0 : i-q1);
	        col_tail = (j-q2<0 ? 0 : j-q2);
	        row_head = (i+q1>=nr ? nr-1 : i+q1);
	        col_head = (j+q2>=nc ? nc-1 : j+q2);
	        REAL(dx)[i+nr*j] = fabs(REAL(kz)[i+nr*col_head] - REAL(kz)[i+nr*col_tail]);
	        REAL(dy)[i+nr*j] = fabs(REAL(kz)[row_head+nr*j] - REAL(kz)[row_tail+nr*j]);
	    }
	}

	/* d'(t) = d(t+1)-d(t) */
	PROTECT(dprimex = allocMatrix(REALSXP, nr, nc));
	PROTECT(dprimey = allocMatrix(REALSXP, nr, nc));

	for(i=0; i<nr; i++) { /* ith row  (y)*/
	    for(j=0;j<nc-1;j++) { /* jth column (x) */
	        REAL(dprimex)[i+j*nr] = REAL(dx)[i+nr*(j+1)]-REAL(dx)[i+nr*j];
	    }
	    REAL(dprimex)[i+nr*j] = REAL(dprimex)[i+nr*(j-1)];
	}
    for(j=0;j<nc;j++) {
    	for(i=0; i<nr-1; i++) {
	        REAL(dprimey)[i+j*nr] = REAL(dy)[i+1+nr*j]-REAL(dy)[i+nr*j];
	    }
	    REAL(dprimey)[nr-1+nr*j] = REAL(dprimey)[nr-2+nr*j];
    }
    
    mx = normalizer(dx, REAL(normalize_prob)[0]);
    my = normalizer(dy, REAL(normalize_prob)[0]);

	PROTECT(ans = allocMatrix(REALSXP, nr, nc));

    for(k=0; k<INTEGER_VALUE(iterations); k++) {

        /* determine bounding box for i,j coordinate */
    	for (i=0; i<nr; i++) {
    	    for(j=0; j<nc; j++) { /* coordinate (x (col),y (row) is [i+nr*j] */
		        if (fabs(REAL(dprimex)[i+nr*j]) < epsilon) { /* dprimex[i,j] = 0 */
			        qx_head = (int) floor(q2*adaptive(REAL(dx)[i+j*nr], mx));
			        qx_tail = (int) floor(q2*adaptive(REAL(dx)[i+j*nr], mx));
		        } else if (REAL(dprimex)[i+nr*j] < 0) {
		        	qx_head = q2;
    		        qx_tail = (int) floor(q2*adaptive(REAL(dx)[i+j*nr], mx));
	    	    } else {
		        	qx_head = (int) floor(q2*adaptive(REAL(dx)[i+j*nr], mx));
			        qx_tail = q2;
			    }

		        if (fabs(REAL(dprimey)[i+nr*j]) < epsilon) { /* dprimey[i,j] = 0 */
			        qy_head = (int) floor(q1*adaptive(REAL(dy)[i+j*nr], my));
			        qy_tail = (int) floor(q1*adaptive(REAL(dy)[i+j*nr], my));
		        } else if (REAL(dprimey)[i+nr*j] < 0) {
		        	qy_head = q1;
    		        qy_tail = (int) floor(q1*adaptive(REAL(dy)[i+j*nr], my));
	    	    } else {
		        	qy_head = (int) floor(q1*adaptive(REAL(dy)[i+j*nr], my));
			        qy_tail = q1;
			    }

			    qx_tail = ((qx_tail) < min_window_length) ? min_window_length : qx_tail;
			    qx_head = ((qx_head) < min_window_length) ? min_window_length : qx_head;
			    qy_tail = ((qy_tail) < min_window_length) ? min_window_length : qy_tail;
			    qy_head = ((qy_head) < min_window_length) ? min_window_length : qy_head;
	            
	            /* bounds are enforced by the box clamps below; the old
	               "check bounds" pair here clamped qy_head against nc (the
	               column count) and qy_tail against the linear index --
	               both wrong, and the first shrank the vertical window on
	               tall narrow matrices whenever q > nc */
                b.upperLeft.x = (j-qx_tail<0 ? 0 : j-qx_tail);
                b.upperLeft.y = (i-qy_tail<0 ? 0 : i-qy_tail);
                b.lowerRight.x = (j+qx_head>=nc ? nc-1 : j+qx_head);
                b.lowerRight.y = (i+qy_head>=nr ? nr-1 : i+qy_head);
                if ((b.lowerRight.y - b.upperLeft.y) + (b.lowerRight.x-b.upperLeft.x))
                    REAL(ans)[i+j*nr]=mavg2d(tmp, b);
                else /* zero sized box - just use avg of one item */
   		            REAL(ans)[i+j*nr]=REAL(tmp)[i+j*nr];
   		    }
	    }
	    copyMatrix(tmp, ans, 0); /* feed this iteration's result into the next, as kza1d/kza3d do */
	}

    UNPROTECT(6);
	return ans;
}

static double averaged(SEXP x, SEXP box_center, SEXP width)
{
	SEXP dim;
	int ndim;
	int count=0;
	double s=0.00;
    SEXP index, box;
    int offset, offset_a;
    int i,j;
    
    dim = GET_DIM(x);
    ndim = LENGTH(dim);
    PROTECT(index = allocVector(INTSXP, LENGTH(dim)));
    PROTECT(box = allocMatrix(INTSXP, LENGTH(dim), 2));
    
    if (isArray(x)) {
        /* bounding box */
        for(i=0; i<ndim; i++) {
            /* box edges */
            INTEGER(box)[i] = (INTEGER(box_center)[i] - INTEGER(width)[i]>0 ? INTEGER(box_center)[i] - INTEGER(width)[i] : 0);
            INTEGER(box)[i+ndim] = (INTEGER(box_center)[i] + INTEGER(width)[i+ndim] <INTEGER(dim)[i] ? INTEGER(box_center)[i] + INTEGER(width)[i+ndim] : INTEGER(dim)[i]-1);
            /* A width-1 extent in one dimension is a legitimate box; the
               loops below handle it (start == end iterates once). The old
               early return here handed back the raw center value whenever
               ANY dimension collapsed -- e.g. kza3d was the identity on
               single-slice (n x m x 1) arrays. */
        }
        for(INTEGER(index)[0]=INTEGER(box)[0]; INTEGER(index)[0]<=INTEGER(box)[ndim]; INTEGER(index)[0]++) {
            for(INTEGER(index)[1]=INTEGER(box)[1]; INTEGER(index)[1]<=INTEGER(box)[1+ndim]; INTEGER(index)[1]++) {
              for(INTEGER(index)[2]=INTEGER(box)[2]; INTEGER(index)[2]<=INTEGER(box)[2+ndim]; INTEGER(index)[2]++) {
/*
        offset += INTEGER(index)[0];
        offset += INTEGER(dim)[0]*INTEGER(index)[1];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(index)[2];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(dim)[2]*INTEGER(index)[3]; 
*/
                    /* find offset into array */
                    for(i=0, offset=0; i<ndim; i++) {
                        for(j=0, offset_a=1; j<i; j++) {
                            offset_a *= INTEGER(dim)[j];
                        }
                        offset += offset_a * INTEGER(index)[i];
                    }
    	            if (R_FINITE(REAL(x)[offset])) {
    	                count++;
    		            s += REAL(x)[offset];
	                }
	            }
            }
        }
	} else error("Input is not a vector or Matrix.");
	
	UNPROTECT(2);
	if (count == 0) return R_NaN;
	return s/count;
}

SEXP kza3d(SEXP v, SEXP kz, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob)
{
    int min_window_length;
	int i, j, k, l;
	SEXP dx, dy, dprimex, dprimey;
	double mx, my, mz;
	SEXP dim, tmp, ans;
	double epsilon;
	SEXP dz, dprimez;
    SEXP box;
    int ndim;
    int offset;
    SEXP index;
    SEXP filter;
    int head, tail;
    int q1, q2, q3;
    
    epsilon = REAL(tolerance)[0];
	if (length(window)<3) {q1 = q2 = q3 = INTEGER_VALUE(window);}
	else {q1 = INTEGER(window)[0]; q2 = INTEGER(window)[1]; q3 = INTEGER(window)[2];}
	
	min_window_length = INTEGER_VALUE(minimum_window_length);
	
    dim = GET_DIM(v);
    ndim = LENGTH(dim);
    if (ndim > 3) error("Too many dimensions, this has not been implemented yet.");
    
	PROTECT(tmp = allocArray(REALSXP, dim));
    copyArray(tmp, v);

	PROTECT(dx = allocArray(REALSXP, dim));
	PROTECT(dy = allocArray(REALSXP, dim));
	PROTECT(dz = allocArray(REALSXP, dim));

	PROTECT(index = allocVector(INTSXP, LENGTH(dim)));
    PROTECT(box = allocMatrix(INTSXP, LENGTH(dim), 2));
	
    /* calculate 
	        dx = |Z(x+q,y,z) - Z(x-q,y,z)| , 
	        dy = |Z(x,y+q,z) - Z(x,y-q),z| , 
	        dz = |Z(x,y,z+q) - Z(x,y,z-q)|
	*/
	for(i=0; i<INTEGER(dim)[0]; i++) { 
	    for(j=0;j<INTEGER(dim)[1];j++) {
	        for(k=0;k<INTEGER(dim)[2];k++) {
	            INTEGER(box)[0] = (i-q1<0 ? 0 : i-q1); /* row tail */
	            INTEGER(box)[1] = (j-q2<0 ? 0 : j-q2); /* col tail */
	            INTEGER(box)[2] = (k-q3<0 ? 0 : k-q3); /* z tail */
	            INTEGER(box)[ndim] = (i+q1>=INTEGER(dim)[0] ? INTEGER(dim)[0]-1 : i+q1); /* row head */
	            INTEGER(box)[ndim+1] = (j+q2>=INTEGER(dim)[1] ? INTEGER(dim)[1]-1 : j+q2); /* col head */
	            INTEGER(box)[ndim+2] = (k+q3>=INTEGER(dim)[2] ? INTEGER(dim)[2]-1 : k+q3); /* z head */
	            offset = i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k;
	            /* head = offset+INTEGER(box)[ndim] 
	            (i+INTEGER(box)[ndim])+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k; */
	            REAL(dx)[offset] = fabs(REAL(kz)[INTEGER(box)[ndim]+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k] - REAL(kz)[INTEGER(box)[0]+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k]);
	            REAL(dy)[offset] = fabs(REAL(kz)[i+INTEGER(dim)[0]*INTEGER(box)[ndim+1]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] - REAL(kz)[i+INTEGER(dim)[0]*INTEGER(box)[1]+INTEGER(dim)[0]*INTEGER(dim)[1]*k]);
	            REAL(dz)[offset] = fabs(REAL(kz)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(box)[ndim+2]] - REAL(kz)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(box)[2]]);
	        }
	    }
	}

	/* d'(t) = d(t+1)-d(t) */
	PROTECT(dprimex = allocArray(REALSXP, dim));
	PROTECT(dprimey = allocArray(REALSXP, dim));
	PROTECT(dprimez = allocArray(REALSXP, dim));

	/* d'(x) = d(x+1,y,z)-d(x,y,z) */
	for(i=0; i<INTEGER(dim)[0]-1; i++) { /* ith row  (y)*/
	    for(j=0;j<INTEGER(dim)[1];j++) { /* jth column (x) */
    	    for(k=0;k<INTEGER(dim)[2];k++) {
	            REAL(dprimex)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = REAL(dx)[(i+1)+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k]-REAL(dx)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k];
	        }
        }
	}
    /* edge -- when dim[0]==1 there is no second-to-last row to copy
       from (the index below underflows to -1); a single-row volume has
       no x-direction structure to measure, so d'(x) is 0 there. */
    if (INTEGER(dim)[0] > 1) {
        for(j=0;j<INTEGER(dim)[1];j++) {
     	    for(k=0;k<INTEGER(dim)[2];k++) {
                REAL(dprimex)[(INTEGER(dim)[0]-1)+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = REAL(dprimex)[(INTEGER(dim)[0]-2)+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k];
            }
        }
    } else {
        for(j=0;j<INTEGER(dim)[1];j++) {
     	    for(k=0;k<INTEGER(dim)[2];k++) {
                REAL(dprimex)[j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = 0.0;
            }
        }
    }
    
	/* d'(y) = d(x,y+1,z)-d(x,y,z) */
   	for(j=0; j<INTEGER(dim)[1]-1; j++) {
        for(i=0;i<INTEGER(dim)[0];i++) {
    	    for(k=0;k<INTEGER(dim)[2];k++) {
    	        REAL(dprimey)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = REAL(dy)[i+INTEGER(dim)[0]*(j+1)+INTEGER(dim)[0]*INTEGER(dim)[1]*k]-REAL(dy)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k];
    	    }
	    }
    }
    /* edge -- source is one COLUMN back (dim[0]*(dim[1]-2)), not the
       literal element two back (dim[0]*dim[1]-2); those only coincide
       when dim[0]==1. For dim[0]>1 the old formula read from the wrong
       row/column (and, on the last z-slice, past the end of the array).
       When dim[1]==1 there is no second-to-last column to copy from
       (the index underflows); a single-column volume has no y-direction
       structure to measure, so d'(y) is 0 there. */
    if (INTEGER(dim)[1] > 1) {
        for(i=0;i<INTEGER(dim)[0];i++) {
            for(k=0;k<INTEGER(dim)[2];k++) {
                REAL(dprimey)[i+(INTEGER(dim)[0]*(INTEGER(dim)[1]-1))+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = REAL(dprimey)[i+(INTEGER(dim)[0]*(INTEGER(dim)[1]-2))+INTEGER(dim)[0]*INTEGER(dim)[1]*k];
            }
        }
    } else {
        for(i=0;i<INTEGER(dim)[0];i++) {
            for(k=0;k<INTEGER(dim)[2];k++) {
                REAL(dprimey)[i+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = 0.0;
            }
        }
    }

	/* d'(z) = d(x,y,z+1)-d(x,y,z) */
    for(k=0;k<INTEGER(dim)[2]-1;k++) {
        for(i=0;i<INTEGER(dim)[0];i++) {
        	for(j=0; j<INTEGER(dim)[1]; j++) {
    	        REAL(dprimez)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k] = REAL(dz)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*(k+1)]-REAL(dz)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k];
    	    }
	    }
    }
    /* edge -- when dim[2]==1 there is no second-to-last slice to copy
       from (the index below underflows to negative offsets, reading
       before the start of the array); a single-slice volume has no
       z-direction structure to measure, so d'(z) is 0 there. */
    if (INTEGER(dim)[2] > 1) {
        for(i=0;i<INTEGER(dim)[0];i++) {
          	for(j=0; j<INTEGER(dim)[1]; j++) {
                REAL(dprimez)[i+INTEGER(dim)[0]*j+(INTEGER(dim)[0]*INTEGER(dim)[1]*(INTEGER(dim)[2]-1))] = REAL(dprimez)[i+INTEGER(dim)[0]*j+(INTEGER(dim)[0]*INTEGER(dim)[1]*(INTEGER(dim)[2]-2))];
            }
        }
    } else {
        for(i=0;i<INTEGER(dim)[0];i++) {
          	for(j=0; j<INTEGER(dim)[1]; j++) {
                REAL(dprimez)[i+INTEGER(dim)[0]*j] = 0.0;
            }
        }
    }
    mx = normalizer(dx, REAL(normalize_prob)[0]);
    my = normalizer(dy, REAL(normalize_prob)[0]);
    mz = normalizer(dz, REAL(normalize_prob)[0]);

	PROTECT(ans = allocArray(REALSXP, dim));
	PROTECT(filter = allocMatrix(INTSXP, ndim, 2));

    for(l=0; l<INTEGER_VALUE(iterations); l++) {

        /* determine bounding box for i,j,k coordinate */
        /*  tail is INTEGER(filter)[0]
            head is INTEGER(filter)[ndim] */
    	for (i=0; i<INTEGER(dim)[0]; i++) {
    	    for(j=0; j<INTEGER(dim)[1]; j++) {
        	    for(k=0; k<INTEGER(dim)[2]; k++) { 
        	        /* x */
                    tail = 0;
                    head = ndim;
		            if (fabs(REAL(dprimex)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k]) < epsilon) { /* dprimex[i,j] = 0 */
			            INTEGER(filter)[head++] = (int) floor(q1*adaptive(REAL(dx)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mx)); /* filter x head */
			            INTEGER(filter)[tail++] = (int) floor(q1*adaptive(REAL(dx)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mx)); /* filter x tail */
		            } else if (REAL(dprimex)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k] < 0) {
		            	INTEGER(filter)[head++] = q1;
    		            INTEGER(filter)[tail++] = (int) floor(q1*adaptive(REAL(dx)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mx));
	    	        } else {
		            	INTEGER(filter)[head++] = (int) floor(q1*adaptive(REAL(dx)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mx));
			            INTEGER(filter)[tail++] = q1;
			        }
			        /* y */
		            if (fabs(REAL(dprimey)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k]) < epsilon) { /* dprimey[i,j,k] = 0 */
		                INTEGER(filter)[head++] = (int) floor(q2*adaptive(REAL(dy)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], my));
			            INTEGER(filter)[tail++] = (int) floor(q2*adaptive(REAL(dy)[i+j*INTEGER(dim)[0]+INTEGER(dim)[1]*INTEGER(dim)[0]*k], my));
		            } else if (REAL(dprimey)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k] < 0) {
		            	INTEGER(filter)[head++] = q2;
    		            INTEGER(filter)[tail++] = (int) floor(q2*adaptive(REAL(dy)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], my));
	    	        } else {
		            	INTEGER(filter)[head++] = (int) floor(q2*adaptive(REAL(dy)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], my));
			            INTEGER(filter)[tail++] = q2;
			        }
			        /* z */
		            if (fabs(REAL(dprimez)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k]) < epsilon) { /* dprimey[i,j,k] = 0 */
		                INTEGER(filter)[head++] = (int) floor(q3*adaptive(REAL(dz)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mz));
			            INTEGER(filter)[tail++] = (int) floor(q3*adaptive(REAL(dz)[i+j*INTEGER(dim)[0]+INTEGER(dim)[1]*INTEGER(dim)[0]*k], mz));
		            } else if (REAL(dprimez)[i+INTEGER(dim)[0]*j+INTEGER(dim)[0]*INTEGER(dim)[1]*k] < 0) {
		            	INTEGER(filter)[head++] = q3;
    		            INTEGER(filter)[tail++] = (int) floor(q3*adaptive(REAL(dz)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mz));
	    	        } else {
		            	INTEGER(filter)[head++] = (int) floor(q3*adaptive(REAL(dz)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k], mz));
			            INTEGER(filter)[tail++] = q3;
			        }
                    
                    for( head=ndim, tail=0; tail<ndim; head++, tail++) {
                        INTEGER(filter)[tail] = ((INTEGER(filter)[tail]) < min_window_length) ? min_window_length : INTEGER(filter)[tail];
                        INTEGER(filter)[head] = ((INTEGER(filter)[head]) < min_window_length) ? min_window_length : INTEGER(filter)[head];
                    }
	                
                    INTEGER(index)[0]=i; INTEGER(index)[1]=j; INTEGER(index)[2]=k;
                    /* EXP x, SEXP box_center, int width */
                    REAL(ans)[i+j*INTEGER(dim)[0]+INTEGER(dim)[0]*INTEGER(dim)[1]*k]=averaged(tmp, index, filter);
                }
   		    }
	    }
        copyArray(tmp, ans);
	}

    UNPROTECT(11);
	return ans;
}

SEXP kza(SEXP x, SEXP kz, SEXP window, SEXP iterations, SEXP min_window, SEXP tol, SEXP normalize_prob)
{
    SEXP ans=R_NilValue, dim;
       
    dim = getAttrib(x, R_DimSymbol);
    //class = getAttrib(x, R_ClassSymbol);
    
    if (LENGTH(x) != LENGTH(kz)) {
        error("The size of the first two arguments do not match.");
        return (ans);
    }    
    if (isArray(x) && LENGTH(dim) >= 3) {
        if (LENGTH(dim) > 3) {
            error("Too many dimensions -- not yet implemented, please contact the author for more info.");
            return(ans);
        }
        PROTECT(ans = kza3d(x, kz, window, iterations, min_window, tol, normalize_prob)); 
        UNPROTECT(1); 
    } else if (isMatrix(x)) {
        PROTECT(ans = kza2d(x, kz, window, iterations, min_window, tol, normalize_prob)); 
        UNPROTECT(1); 
    } else if (isVector(x)) {
        PROTECT(ans = kza1d(x, kz, window, iterations, min_window, tol, normalize_prob)); 
        UNPROTECT(1); 
    }
    else error("Input is not a matrix or vector.");
    
    return(ans);
}

