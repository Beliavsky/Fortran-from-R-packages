#include <R.h>
#include <Rmath.h>

#include "clue.h"

double **clue_vector_to_square_matrix(double *x, Sint n)
{
    double **data, *val;
    Sint i, j;
    data = (double **) R_alloc(n, sizeof(double));
    for(i = 0; i < n; i++) {
        data[i] = (double *) R_alloc(n, sizeof(double));
        val = x + i;
        for(j = 0; j < n; j++, val += n)
            data[i][j] = *val;
    }
    return(data);
}

static int
clue_sign(double x)
{
    if(x == 0) return(0);
    return((x > 0) ? 1 : -1);
}

void
clue_dissimilarity_count_inversions(double *x, double *y, Sint *n,
				    double *count)
{
    Sint i, j;
    for(i = 0; i < *n; i++)
	for(j = 0; j < *n; j++)
	    if((clue_sign(x[i] - x[j]) * clue_sign(y[i] - y[j])) < 0)
		(*count)++;
}
