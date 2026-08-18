#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

SEXP gcrp_c(SEXP n_, SEXP alpha_, SEXP theta_) {
    int n = INTEGER(n_)[0];
    double alpha = REAL(alpha_)[0];
    double theta = REAL(theta_)[0];
    
    if (alpha < 0.0 || alpha >= 1.0)
        error("alpha must be in [0,1)");
    if (theta < -alpha)
        error("theta must be >= -alpha");
    
    SEXP assignment_ = PROTECT(allocVector(INTSXP, n));
    int *assignment = INTEGER(assignment_);

    SEXP table_sizes_ = PROTECT(allocVector(INTSXP, n));
    int *table_sizes = INTEGER(table_sizes_);

    int n_tables = 1;
    assignment[0] = 1;
    table_sizes[0] = 1;

    GetRNGstate();

    for (int i = 1; i < n; i++) {
        double total = theta + i;
        double u = unif_rand() * total;

        double cum = 0.0;
        int chosen = n_tables;

        for (int k = 0; k < n_tables; k++) {
	  cum += (table_sizes[k] - alpha);
	  if (u <= cum) {
	    chosen = k;
	    break;
	  }
        }

        if (chosen == n_tables) {
            table_sizes[n_tables] = 1;
            n_tables++;
        } else {
            table_sizes[chosen]++;
        }

        assignment[i] = chosen + 1;
    }

    PutRNGstate();

    UNPROTECT(2);
    return assignment_;
}
