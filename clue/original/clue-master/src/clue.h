#ifndef _CLUE_H
#define _CLUE_H

#include <R.h>

typedef int Sint;

void solve_LSAP(double *c, Sint *n, Sint *p);

double **clue_vector_to_square_matrix(double *x, Sint n);

void clue_dissimilarity_count_inversions(double *x, double *y, Sint *n, double *count);

void deviation_from_ultrametricity(double *x, int *n, double *v, int *max); 
void deviation_from_ultrametricity_gradient(double *x, int *n, double *out);
void deviation_from_additivity(double *x, int *n, double *v, int *max);
void deviation_from_additivity_gradient(double *x, int *n, double *out);

void ls_fit_ultrametric_by_iterative_reduction(double *d, int *n, int *order, int *maxiter, int *iter, double *tol, int *verbose);
void ls_fit_ultrametric_by_iterative_projection(double *d, int *n, int *order, int *maxiter, int *iter, double *tol, int *verbose);
void ls_fit_addtree_by_iterative_reduction(double *d, int *n, int *order, int *maxiter, int *iter, double *tol, int *verbose);
void ls_fit_addtree_by_iterative_projection(double *d, int *n, int *order, int *maxiter, int *iter, double *tol, int *verbose);

#endif
