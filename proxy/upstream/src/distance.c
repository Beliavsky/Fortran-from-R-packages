#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

// extends the code from stats/src/distance.c
// to handle auto- and cross-distances.
//
// note: the runtime in the symmetric case is
//       not always optimal.
//
// ceeboo 2007, 2008, 2009, 2014, 2016

#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))

typedef double (* DFUN)(double *, double *, int, int, int);

static double dfp = 1;

static double minkowski(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            dev = (*x - *y);
            if (!ISNAN(dev)) {
                dist += R_pow(fabs(dev), dfp);
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return R_pow(dist, 1.0/dfp);
}

static double euclidean(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    dev = (*x - *y);
            if (!ISNAN(dev)) {
		dist += dev * dev;
		count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}

static double maximum(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist = -DBL_MAX;
    for (j = 0 ; j < nc ; j++) {
        if (both_non_NA(*x, *y)) {
            dev = fabs(*x - *y);
            if (!ISNAN(dev)) {
                if (dev > dist)
                    dist = dev;
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    //if (count != nc) dist /= ((double)count/nc);
    return dist;
}

static double manhattan(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            dev = fabs(*x - *y);
            if (!ISNAN(dev)) {
                dist += dev;
                count++;
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return dist;
}

// fixme: NA for two all-zero vectors

static double canberra(double *x, double *y, int nx, int ny, int nc)
{
    double dev, dist, sum, diff;
    int count, j;

    count = 0;
    dist  = 0;
    for (j = 0 ;j < nc; j++) {
        if (both_non_NA(*x, *y)) {
            sum = fabs(*x + *y);
            diff = fabs(*x - *y);
            if (sum > DBL_MIN || diff > DBL_MIN) {
                dev = diff/sum;
                if (!ISNAN(dev) ||
                   (!R_FINITE(diff) && diff == sum &&
                    /* use Inf = lim x -> oo */ (dev = 1.))) {
                    dist += dev;
                    count++;
                }
            }
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return dist;
}

// FIXME why treat not both finite as NA?
static double binary(double *x, double *y, int nx, int ny, int nc)
{
    int total, count, dist;
    int j;

    total = count = dist = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x || *y) {
		count++;
		if (!(*x && *y))
		    dist++;
            }
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (count == 0) return 0;
    return (double) dist / count;
}

static double matching(double *x, double *y, int nx, int ny, int nc)
{
    int total, count;
    int j;

    total = count = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x != *y) 
		count++;
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (count == 0) return 0;
    return (double) count / total;
}

static double fuzzy(double *x, double *y, int nx, int ny, int nc)
{
    double dist, smax, smin;
    int count, j;

    count = 0;
    smax = smin = 0;
    for (j = 0; j < nc; j++) {
	if (both_FINITE(*x, *y)) {
	    if (*x > *y) {
		smax += *x;
		smin += *y;
	    }
	    else {
		smax += *y;
		smin += *x;
	    }
	    count++;
        }
        x += nx;
        y += ny;
    }
    if (count == 0) return NA_REAL;
    if (!R_FINITE(smin)) return NA_REAL;
    dist = smin / smax;
    if (ISNAN(dist)) return 0;
    return 1-dist;
}

static double mutual(double *x, double *y, int nx, int ny, int nc)
{
    double dist;
    int total, count, cx, cy, j;

    total = count = cx = cy = 0;
    for (j = 0; j < nc; j++) {
	if (both_non_NA(*x, *y)) {
	    if (*x && *y) 
		count++;
	    cx += (*x && 1);
	    cy += (*y && 1);
            total++;
        }
        x += nx;
        y += ny;
    }
    if (total == 0) return NA_REAL;
    if (cx == 0 || cy == 0 || cx == total || cy == total) return 0;
    dist = 0;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total);
    cy = total - cy;
    count = cx - count;
    if (count > 0) 
	dist += (double) count / total * log((double) count / cx / cy * total);
    cx = total - cx;
    count = cy - count;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total);
    cy = total - cy;
    count = cx - count;
    if (count > 0)
	dist += (double) count / total * log((double) count / cx / cy * total); 
    if (total != nc) dist /= ((double)total/nc);
    return dist ;
}


// wrapper

static SEXP dists(SEXP R_x, SEXP R_y, SEXP R_d, DFUN f, SEXP R_p) {
    if (!isMatrix(R_x))
	error("'x' not of class matrix");
    if (!isNull(R_y) && !isMatrix(R_y))
	error("'y' not of class matrix");
    if (TYPEOF(R_d) != LGLSXP)
	error("'d' not of type logical");
    int i, j, nx, ny, nc, m = 0;
    SEXP x = R_x, y = R_y, r;
    R_xlen_t l;

    if (!isNull(R_p))			    // fixme: check?
	dfp = *REAL(R_p);
    
    if (isNull(y)) 
	y = x;
    else
    if (LOGICAL(R_d)[0] == TRUE)
	m = 2;
    else
	m = 1;				    // return matrix
    
    nc = INTEGER(GET_DIM(x))[1];
    
    if (INTEGER(GET_DIM(y))[1] != nc)
	error("invalid number of columns");

    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];
   
    if (m == 2 && nx != ny)
	error("invalid number of rows for pairwise mode");

    if (TYPEOF(x) != REALSXP) {
	PROTECT(x = coerceVector(R_x, REALSXP));

	if (isNull(R_y) || R_x == R_y)
	    y = x;
    }

    if (TYPEOF(y) != REALSXP) 
	PROTECT(y = coerceVector(R_y, REALSXP));

    if (m == 0) {
	SEXP d;

	l = nx;	
	PROTECT(r = allocVector(REALSXP, l*(l-1)/2));
	setAttrib(r, install("Size"), PROTECT(ScalarInteger(nx)));
	UNPROTECT(1);
	
	if (!isNull(d = getAttrib(x, R_DimNamesSymbol)))
	    setAttrib(r, install("Labels"), VECTOR_ELT(d, 0));
	// fixme: package?
	setAttrib(r, R_ClassSymbol, PROTECT(mkString("dist")));
	UNPROTECT(1);
    } else
    if (m == 1) {
	SEXP d1, d2;

	PROTECT(r = allocMatrix(REALSXP, nx, ny));

	d1 = getAttrib(x, R_DimNamesSymbol);
	d2 = getAttrib(y, R_DimNamesSymbol);
	if (!isNull(d1) || !isNull(d2)) {
	    SEXP d;
	    
	    setAttrib(r, R_DimNamesSymbol, PROTECT(d = allocVector(VECSXP, 2)));
	    UNPROTECT(1);
	    SET_VECTOR_ELT(d, 0, isNull(d1) ? d1 : VECTOR_ELT(d1, 0));
	    SET_VECTOR_ELT(d, 1, isNull(d2) ? d2 : VECTOR_ELT(d2, 0));
	}
    } else
	PROTECT(r = allocVector(REALSXP, nx));

    l = 0;
    for (j = 0; j < ny; j++) {
	if (m == 2)
	    REAL(r)[l++] = f(REAL(x)+j, REAL(y)+j, nx, ny, nc);
	else
	    for (i = (m == 0) ? j+1 : 0; i < nx; i++)
		REAL(r)[l++] = f(REAL(x)+i, REAL(y)+j, nx, ny, nc);
	R_CheckUserInterrupt();
    }

    UNPROTECT(1);
    if (x != R_x)
	UNPROTECT(1);
    if (!isNull(R_y) && R_y != R_x && y != R_y)
	UNPROTECT(1);

    return r;
}

// R wrappers

SEXP R_minkowski_dist(SEXP x, SEXP y, SEXP d, SEXP p) {
    if (isNull(p))
	error("'p' invalid");
    return dists(x, y, d, minkowski, p);
}

SEXP R_euclidean_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, euclidean, R_NilValue);
}

SEXP R_maximum_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, maximum, R_NilValue);
}

SEXP R_manhattan_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, manhattan, R_NilValue);
}

SEXP R_canberra_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, canberra, R_NilValue);
}

SEXP R_binary_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, binary, R_NilValue);
}

SEXP R_matching_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, matching, R_NilValue);
}

SEXP R_fuzzy_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, fuzzy, R_NilValue);
}

SEXP R_mutual_dist(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, mutual, R_NilValue);
}

// 2017/10
SEXP R_bjaccard(SEXP x, SEXP y, SEXP d) {
    SEXP r = dists(x, y, d, binary, R_NilValue);
    for (R_xlen_t k = 0; k < XLENGTH(r); k++) {
	double z = REAL(r)[k];
	if (!ISNAN(z))
	    REAL(r)[k] = 1 - z;
    }
    return r;
}

static double ebinary(double *x, double *y, int nx, int ny, int nc)
{
    double dev, prod, dist, xy;
    int count;
    int j;

    dist = xy = 0;
    count = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
	    dev = (*x - *y);
	    prod = *x * *y;
	    if (!ISNAN(dev) && !ISNAN(prod)) {
		dist += dev * dev;
		xy += prod;
		count++;
	    }
        }
	x += nx;
	y += ny;
    }
    if (count == 0) return NA_REAL;
    if (!R_FINITE(xy)) return NA_REAL;
    dist = dist / dfp + xy;
    xy /= dist;
    if (ISNAN(xy))
	if (dist < DBL_MIN)
	    return 1;
	else
	    return NA_REAL;
    return xy;
}

SEXP R_ess2(SEXP x, SEXP y, SEXP d) {
    dfp = .5;
    return dists(x, y, d, ebinary, R_NilValue);
}

SEXP R_ejaccard(SEXP x, SEXP y, SEXP d) {
    dfp = 1;
    return dists(x, y, d, ebinary, R_NilValue);
}

SEXP R_edice(SEXP x, SEXP y, SEXP d) {
    dfp = 2;
    return dists(x, y, d, ebinary, R_NilValue);
}

static double cosine(double *x, double *y, int nx, int ny, int nc)
{
    double prod, xy, xx, yy;
    int count;
    int j;

    xy = xx = yy = 0;
    count = 0;
    for (j = 0; j < nc; j++) {
        if (both_non_NA(*x, *y)) {
	    prod = *x * *y;
	    if (!ISNAN(prod)) {
		xy += prod;
		xx += *x * *x;
		yy += *y * *y;
		count++;
	    }
        }
	x += nx;
	y += ny;
    }
    if (count == 0) return NA_REAL;
    if (!R_FINITE(xy)) return NA_REAL;
    xy /= sqrt(xx) * sqrt(yy);
    if (ISNAN(xy))
	if (xx < DBL_MIN && yy < DBL_MIN)
	    return 1;
	else
	    if (xx < DBL_MIN || yy < DBL_MIN)
		return 0;
	    else
		return NA_REAL;
    return xy;
}

SEXP R_cosine(SEXP x, SEXP y, SEXP d) {
    return dists(x, y, d, cosine, R_NilValue);
}


//
