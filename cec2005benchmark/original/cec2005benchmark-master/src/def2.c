/* Function definitions of utility functions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "sub.h"
#include "rand.h"

/* Code to allocate memory to global variables being used in evaluation of functions */
void allocate_memory()
{
    int i, j, k;
    norm_x = (double *) malloc(nreal * sizeof(double));
    norm_f = (double *) malloc(nfunc * sizeof(double));
    trans_x = (double *) malloc(nreal * sizeof(double));
    basic_f = (double *) malloc(nfunc * sizeof(double));
    temp_x1 = (double *) malloc(nreal * sizeof(double));
    temp_x2 = (double *) malloc(nreal * sizeof(double));
    temp_x3 = (double *) malloc(nreal * sizeof(double));
    temp_x4 = (double *) malloc(nreal * sizeof(double));
    weight = (double *) malloc(nfunc * sizeof(double));
    sigma = (double *) malloc(nfunc * sizeof(double));
    lambda = (double *) malloc(nfunc * sizeof(double));
    bias = (double *) malloc(nfunc * sizeof(double));
    o = (double **) malloc(nfunc * sizeof(double));
    l = (double ***) malloc(nfunc * sizeof(double));
    g = (double **) malloc(nreal * sizeof(double));
    for (i = 0; i < nfunc; i++) {
        o[i] = (double *) malloc(nreal * sizeof(double));
        l[i] = (double **) malloc(nreal * sizeof(double));
    }
    for (i = 0; i < nreal; i++) {
        g[i] = (double *) malloc(nreal * sizeof(double));
    }
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            l[i][j] = (double *) malloc(nreal * sizeof(double));
        }
    }
    /* Do some trivial (common) initialization here itself */
    C = 2000.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 5.0;
        trans_x[i] = 0.0;
        temp_x1[i] = 0.0;
        temp_x2[i] = 0.0;
        temp_x3[i] = 0.0;
        temp_x4[i] = 0.0;
        for (j = 0; j < nreal; j++) {
            if (i == j) {
                g[i][j] = 1.0;
            } else {
                g[i][j] = 0.0;
            }
        }
    }
    for (i = 0; i < nfunc; i++) {
        basic_f[i] = 0.0;
        norm_f[i] = 0.0;
        weight[i] = 1.0 / (double) nfunc;
        sigma[i] = 1.0;
        lambda[i] = 1.0;
        bias[i] = 100.0 * (double) i;
        for (j = 0; j < nreal; j++) {
            o[i][j] = 0.0;
            for (k = 0; k < nreal; k++) {
                if (j == k) {
                    l[i][j][k] = 1.0;
                } else {
                    l[i][j][k] = 0.0;
                }
            }
        }
    }
    A_f5 = NULL;
    B_f5 = NULL;
    A_f12 = NULL;
    B_f12 = NULL;
    alpha_f12 = NULL;

    return;
}

/* Code to transform a variable vector based on function index 'count' */
void transform(double *x, int count)
{
    int i, j;
    for (i = 0; i < nreal; i++) {
        temp_x1[i] = x[i] - o[count][i];
    }
    for (i = 0; i < nreal; i++) {
        temp_x2[i] = temp_x1[i] / lambda[count];
    }
    for (j = 0; j < nreal; j++) {
        temp_x3[j] = 0.0;
        for (i = 0; i < nreal; i++) {
            temp_x3[j] += g[i][j] * temp_x2[i];
        }
    }
    for (j = 0; j < nreal; j++) {
        trans_x[j] = 0.0;
        for (i = 0; i < nreal; i++) {
            trans_x[j] += l[count][i][j] * temp_x3[i];
        }
    }
    return;
}

/* Code to transform a vector (with elements 5.0) based on function index 'count' */
void transform_norm(int count)
{
    int i, j;
    for (i = 0; i < nreal; i++) {
        temp_x2[i] = 5.0 / lambda[count];
    }
    for (j = 0; j < nreal; j++) {
        temp_x3[j] = 0.0;
        for (i = 0; i < nreal; i++) {
            temp_x3[j] += g[i][j] * temp_x2[i];
        }
    }
    for (j = 0; j < nreal; j++) {
        trans_x[j] = 0.0;
        for (i = 0; i < nreal; i++) {
            trans_x[j] += l[count][i][j] * temp_x3[i];
        }
    }
    return;
}

/* Code to compute the weights for a variable vector */
void calc_weight(double *x)
{
    int i, j;
    double sum;
    double max;
    max = -INF;
    for (i = 0; i < nfunc; i++) {
        sum = 0.0;
        for (j = 0; j < nreal; j++) {
            sum += (x[j] - o[i][j]) * (x[j] - o[i][j]);
        }
        weight[i] = exp(-(sum) / (2.0 * nreal * sigma[i] * sigma[i]));
        max = maximum(max, weight[i]);
    }
    sum = 0.0;
    for (i = 0; i < nfunc; i++) {
        if (weight[i] != max) {
            weight[i] *= (1.0 - pow(max, 10.0));
        }
        sum += weight[i];
    }
    if (sum == 0.0) {
        for (i = 0; i < nfunc; i++) {
            weight[i] = 1.0 / (double) nfunc;
        }
    } else {
        for (i = 0; i < nfunc; i++) {
            weight[i] /= sum;
        }
    }
    return;
}

/* Code to free the allocated memory */
void free_memory()
{
    int i, j;
    free(norm_x);
    free(norm_f);
    free(trans_x);
    free(basic_f);
    free(temp_x1);
    free(temp_x2);
    free(temp_x3);
    free(temp_x4);
    free(weight);
    free(sigma);
    free(lambda);
    free(bias);
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            free(l[i][j]);
        }
    }
    for (i = 0; i < nfunc; i++) {
        free(o[i]);
        free(l[i]);
    }
    for (i = 0; i < nreal; i++) {
        free(g[i]);
    }
    free(o);
    free(l);
    free(g);
    if (A_f5 != NULL) {
        for (i = 0; i < nreal; i++) {
            free(A_f5[i]);
        }
        free(A_f5);
    }
    if (B_f5 != NULL)free(B_f5);
    if (A_f12 != NULL) {
        for (i = 0; i < nreal; i++) {
            free(A_f12[i]);
        }
        free(A_f12);
    }
    if (B_f12 != NULL) {
        for (i = 0; i < nreal; i++) {
            free(B_f12[i]);
        }
        free(B_f12);
    }
    free(alpha_f12);
    return;
}
