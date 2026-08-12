#include <R.h>
#include <Rmath.h>

#include "clue.h"

static int iwork3[3];
static int iwork4[4];

static void
isort3(int *i, int *j, int *k)
{
    iwork3[0] = *i;
    iwork3[1] = *j;
    iwork3[2] = *k;
    R_isort(iwork3, 3);
    *i = iwork3[0];
    *j = iwork3[1];
    *k = iwork3[2];
}

static void
isort4(int *i, int *j, int *k, int *l)
{
    iwork4[0] = *i;
    iwork4[1] = *j;
    iwork4[2] = *k;
    iwork4[3] = *l;
    R_isort(iwork4, 4);
    *i = iwork4[0];
    *j = iwork4[1];
    *k = iwork4[2];
    *l = iwork4[3];
}

void
deviation_from_ultrametricity(double *x, int *n, double *v, int *max)
{
    double **D, p, delta, A, B, C;
    int i, j, k;

    D = clue_vector_to_square_matrix(x, *n);

    p = 0;
    for(i = 0; i < *n - 2; i++)
        for(j = i + 1; j < *n - 1; j++) {
	    A = D[i][j];
            for(k = j + 1; k < *n; k++) {
		B = D[i][k];
		C = D[j][k];
		if((A <= B) && (A <= C))
		    delta = C - B;
		else if(B <= C)
		    delta = A - C;
		else
		    delta = B - A;
		if(*max)
		    p = fmax2(p, fabs(delta));
		else
		    p += delta * delta;
	    }
	}

    *v = p;
}

void
deviation_from_ultrametricity_gradient(double *x, int *n, double *out)
{
    double **D, **G, A, B, C, delta;
    int i, j, k;

    D = clue_vector_to_square_matrix(x, *n);
    G = clue_vector_to_square_matrix(out, *n);

    for(i = 0; i < *n - 2; i++)
        for(j = i + 1; j < *n - 1; j++) {
	    A = D[i][j];
            for(k = j + 1; k < *n; k++) {
		B = D[i][k];
		C = D[j][k];
		if((A <= B) && (A <= C)) {
		    delta = 2 * (B - C);
		    G[i][k] += delta;
		    G[j][k] -= delta;
		}
		else if(B <= C) {
		    delta = 2 * (C - A);
		    G[j][k] += delta;
		    G[i][j] -= delta;
		}
		else {
		    delta = 2 * (A - B);
		    G[i][j] += delta;
		    G[i][k] -= delta;
		}
	    }
	}

    for(i = 0; i < *n; i++)
	for(j = 0; j < *n; j++)
	    *out++ = G[i][j];

}
    
void
deviation_from_additivity(double *x, int *n, double *v, int *max)
{
    double **D, p, delta, A, B, C;
    int i, j, k, l;

    D = clue_vector_to_square_matrix(x, *n);

    p = 0;
    for(i = 0; i < *n - 3; i++)
	for(j = i + 1; j < *n - 2; j++)
	    for(k = j + 1; k < *n - 1; k++)
		for(l = k + 1; l < *n; l++) {
		    A = D[i][j] + D[k][l];
		    B = D[i][k] + D[j][l];
		    C = D[i][l] + D[j][k];
		    if((A <= B) && (A <= C))
			delta = (C - B);
		    else if(B <= C)
			delta = (A - C);
		    else
			delta = (B - A);
		    if(*max)
			p = fmax2(p, fabs(delta));
		    else
			p += delta * delta;
		}

    *v = p;
}
		    
void
deviation_from_additivity_gradient(double *x, int *n, double *out)
{
    double **D, **G, A, B, C, delta;
    int i, j, k, l;

    D = clue_vector_to_square_matrix(x, *n);
    G = clue_vector_to_square_matrix(out, *n);

    for(i = 0; i < *n - 3; i++)
	for(j = i + 1; j < *n - 2; j++)
	    for(k = j + 1; k < *n - 1; k++)
		for(l = k + 1; l < *n; l++) {
		    A = D[i][j] + D[k][l];
		    B = D[i][k] + D[j][l];
		    C = D[i][l] + D[j][k];
		    if((A <= B) && (A <= C)) {
			delta = 2 * (B - C);
			G[i][l] -= delta;
			G[j][k] -= delta;
			G[i][k] += delta;
			G[j][l] += delta;
		    }
		    else if(B <= C) {
			delta = 2 * (C - A);
			G[i][l] += delta;
			G[j][k] += delta;
			G[i][j] -= delta;
			G[k][l] -= delta;
		    }
		    else {
			delta = 2 * (A - B);
			G[i][k] -= delta;
			G[j][l] -= delta;
			G[i][j] += delta;
			G[k][l] += delta;
		    }
		}

    for(i = 0; i < *n; i++)
	for(j = 0; j < *n; j++)
	    *out++ = G[i][j];

}

void
ls_fit_ultrametric_by_iterative_reduction(double *d, int *n, int *order,
					  int *maxiter, int *iter,
					  double *tol, int *verbose)
{
    double A, B, C, **D, DQ, delta, tmp;
    int i, i1, j, j1, k, k1, N3;

    D = clue_vector_to_square_matrix(d, *n);
    /* And initialize the upper half of D ("work array") to 0.
       (Yes, this could be done more efficiently by just propagating the
       veclh dist representation.)
    */
    for(i = 0; i < *n - 1; i++)
	for(j = i + 1; j < *n; j++)
	    D[i][j] = 0;

    N3 = (*n - 2);

    for(*iter = 0; *iter < *maxiter; (*iter)++) {
	if(*verbose) Rprintf("Iteration: %d, ", *iter);
	for(i1 = 0; i1 < *n - 2; i1++)
	    for(j1 = i1 + 1; j1 < *n - 1; j1++)
		for(k1 = j1 + 1; k1 < *n; k1++) {
		    i = order[i1];
		    j = order[j1];
		    k = order[k1];
		    isort3(&i, &j, &k);
		    A = D[j][i];
		    B = D[k][i];
		    C = D[k][j];
		    /*
		      <FIXME>
		      B & G have a divisor of 2 for case 1 and 4 for
		      cases 2 and 3 ... clearly, we should use the same
		      in all cases, but should it be 2 or 4?
		      </FIXME>
		    */
		    if((A <= B) && (A <= C)) {
			/* Case 1: 5080 */
			DQ = (C - B) / 2;
			D[i][k] += DQ;
			D[j][k] -= DQ;
		    }
		    else if(B <= C) {
			/* Case 2: 5100 */
			DQ = (C - A) / 2;
			D[i][j] += DQ;
			D[j][k] -= DQ;
		    }
		    else {
			/* Case 3: 5120 */
			DQ = (B - A) / 2;
			D[i][j] += DQ;
			D[i][k] -= DQ;
		    }
		}

	delta = 0;
	for(i = 0; i < *n - 1; i++)
	    for(j = i + 1; j < *n; j++) {
		tmp = D[i][j] / N3;
		D[j][i] += tmp;
		D[i][j] = 0;
		delta += fabs(tmp);
	    }

	if(*verbose) Rprintf("change: %f\n", delta);
	if(delta < *tol)
	    break;
	   
    }

    /* And now write results back.
       Could make this more efficient, of course ...
    */
    for(j = 0; j < *n; j++)
	for(i = 0; i < *n; i++)
	    *d++ = D[i][j];
}

void
ls_fit_ultrametric_by_iterative_projection(double *d, int *n, int *order,
					   int *maxiter, int *iter,
					   double *tol, int *verbose)
{
    double A, B, C, **D, delta;
    int i, i1, j, j1, k, k1;    

    D = clue_vector_to_square_matrix(d, *n);

    for(*iter = 0; *iter < *maxiter; (*iter)++) {
	if(*verbose) Rprintf("Iteration: %d, ", *iter);
	delta = 0;
	for(i1 = 0; i1 < *n - 2; i1++)
	    for(j1 = i1 + 1; j1 < *n - 1; j1++)
		for(k1 = j1 + 1; k1 < *n; k1++) {
		    i = order[i1];
		    j = order[j1];
		    k = order[k1];
		    isort3(&i, &j, &k);
		    A = D[i][j];		    
		    B = D[i][k];
		    C = D[j][k];
		    if((A <= B) && (A <= C)) {
			D[i][k] = D[j][k] = (B + C) / 2;
			delta += fabs(B - C);
		    }
		    else if(B <= C) {
			D[i][j] = D[j][k] = (C + A) / 2;
			delta += fabs(C - A);
		    }
		    else {
			D[i][j] = D[i][k] = (A + B) / 2;
			delta += fabs(A - B);
		    }

		}
	
	if(*verbose) Rprintf("change: %f\n", delta);
	if(delta < *tol)
	    break;
    }
	
    for(i = 0; i < *n - 1; i++)
	for(j = i + 1; j < *n; j++)
	    D[j][i] = D[i][j];

    /* And now write results back.
       Could make this more efficient, of course ...
    */
    for(j = 0; j < *n; j++)
	for(i = 0; i < *n; i++)
	    *d++ = D[i][j];
}

void
ls_fit_addtree_by_iterative_reduction(double *d, int *n, int *order,
				      int *maxiter, int *iter,
				      double *tol, int *verbose)
{
    /* Once we have ls_fit_ultrametric_by_iterative_reduction() we can
       always do this as well ...

       See page 67f in Barthelemy and Guenoche.
     
     */
    
    double A, B, C, **D, DQ, delta, tmp, N3;
    int i, i1, j, j1, k, k1, l, l1;

    D = clue_vector_to_square_matrix(d, *n);
    /* And initialize the upper half of D ("work array") to 0.
       (Yes, this could be done more efficiently by just propagating the
       veclh dist representation.)
    */
    for(i = 0; i < *n - 1; i++)
	for(j = i + 1; j < *n; j++)
	    D[i][j] = 0;

    N3 = (*n - 2) * (*n - 3) / 2;

    for(*iter = 0; *iter < *maxiter; (*iter)++) {
	if(*verbose) Rprintf("Iteration: %d, ", *iter);
	for(i1 = 0; i1 < *n - 3; i1++)
	    for(j1 = i1 + 1; j1 < *n - 2; j1++)
		for(k1 = j1 + 1; k1 < *n - 1; k1++)
		    for(l1 = k1 + 1; l1 < *n; l1++) {
			i = order[i1];
			j = order[j1];
			k = order[k1];
			l = order[l1];
			isort4(&i, &j, &k, &l);
			A = D[j][i] + D[l][k];
			B = D[k][i] + D[l][j];
			C = D[l][i] + D[k][j];
			if((A <= B) && (A <= C)) {
			    /* Case 1: 5090 */
			    DQ = (C - B) / 4;
			    D[i][l] -= DQ;
			    D[j][k] -= DQ;
			    D[i][k] += DQ;
			    D[j][l] += DQ;
			}
			else if(B <= C) {
			    /* Case 2: 5120 */
			    DQ = (A - C) / 4;
			    D[i][l] += DQ;
			    D[j][k] += DQ;
			    D[i][j] -= DQ;
			    D[k][l] -= DQ;
			}
			else {
			    /* Case 3: 5150 */
			    DQ = (B - A) / 4;
			    D[i][k] -= DQ;
			    D[j][l] -= DQ;
			    D[i][j] += DQ;
			    D[k][l] += DQ;
			}
		    }
	
	delta = 0;
	for(i = 0; i < *n - 1; i++)
	    for(j = i + 1; j < *n; j++) {
		tmp = D[i][j] / N3;
		D[j][i] += tmp;
		D[i][j] = 0;
		delta += fabs(tmp);
	    }

	if(*verbose) Rprintf("change: %f\n", delta);	
	if(delta < *tol)
	    break;
	   
    }

    /* And now write results back.
       Could make this more efficient, of course ...
    */
    for(j = 0; j < *n; j++)
	for(i = 0; i < *n; i++)
	    *d++ = D[i][j];
}

void
ls_fit_addtree_by_iterative_projection(double *d, int *n, int *order,
				       int *maxiter, int *iter,
				       double *tol, int *verbose)
{
    double A, B, C, **D, DQ, delta;
    int i, i1, j, j1, k, k1, l, l1;

    D = clue_vector_to_square_matrix(d, *n);

    for(*iter = 0; *iter < *maxiter; (*iter)++) {
	delta = 0;
	if(*verbose) Rprintf("Iteration: %d, ", *iter);
	for(i1 = 0; i1 < *n - 3; i1++)
	    for(j1 = i1 + 1; j1 < *n - 2; j1++)
		for(k1 = j1 + 1; k1 < *n - 1; k1++)
		    for(l1 = k1 + 1; l1 < *n; l1++) {
			i = order[i1];
			j = order[j1];
			k = order[k1];
			l = order[l1];
			isort4(&i, &j, &k, &l);
			A = D[i][j] + D[k][l];
			B = D[i][k] + D[j][l];
			C = D[i][l] + D[j][k];
			if((A <= B) && (A <= C)) {
			    DQ = (C - B) / 4;
			    D[i][l] -= DQ;
			    D[j][k] -= DQ;
			    D[i][k] += DQ;
			    D[j][l] += DQ;
			    delta += fabs(C - B);
			}
			else if(B <= C) {
			    DQ = (A - C) / 4;
			    D[i][l] += DQ;
			    D[j][k] += DQ;
			    D[i][j] -= DQ;
			    D[k][l] -= DQ;
			    delta += fabs(A - C);
			}
			else {
			    DQ = (B - A) / 4;
			    D[i][k] -= DQ;
			    D[j][l] -= DQ;
			    D[i][j] += DQ;
			    D[k][l] += DQ;
			    delta += fabs(B - A);
			}
		    }
	
	if(*verbose) Rprintf("change: %f\n", delta);	
	if(delta < *tol)
	    break;
	   
    }

    for(i = 0; i < *n - 1; i++)
	for(j = i + 1; j < *n; j++)
	    D[j][i] = D[i][j];

    /* And now write results back.
       Could make this more efficient, of course ...
    */
    for(j = 0; j < *n; j++)
	for(i = 0; i < *n; i++)
	    *d++ = D[i][j];
}
