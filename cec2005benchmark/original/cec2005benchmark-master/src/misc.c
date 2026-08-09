/* Some auxiliary functions (not part of any algorithm or procedure) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "sub.h"
#include "rand.h"

/* Function to return the maximum of two variables */
double maximum(double a, double b)
{
    if (a > b) {
        return (a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double minimum(double a, double b)
{
    if (a < b) {
        return (a);
    }
    return (b);
}

/* Function to return the modulus of a vector */
double modulus(double *x, int n)
{
    int i;
    double res;
    res = 0.0;
    for (i = 0; i < n; i++) {
        res += x[i] * x[i];
    }
    return (sqrt(res));
}

/* Function to return the dot product of two vecors */
double dot(double *a, double *b, int n)
{
    int i;
    double res;
    res = 0.0;
    for (i = 0; i < n; i++) {
        res += a[i] * b[i];
    }
    return (res);
}

/* Function to return the mean of n variables */
double mean(double *x, int n)
{
    int i;
    double res;
    res = 0.0;
    for (i = 0; i < n; i++) {
        res += x[i];
    }
    return (res / (double) n);
}
