
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP R_minkowski_dist(SEXP x, SEXP y, SEXP d, SEXP p);
extern SEXP R_euclidean_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_maximum_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_manhattan_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_canberra_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_binary_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_matching_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_fuzzy_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_mutual_dist(SEXP x, SEXP y, SEXP d);
extern SEXP R_bjaccard(SEXP R_x, SEXP R_y, SEXP R_d);
extern SEXP R_ejaccard(SEXP R_x, SEXP R_y, SEXP R_d);
extern SEXP R_edice(SEXP R_x, SEXP R_y, SEXP R_d);
extern SEXP R_cosine(SEXP R_x, SEXP R_y, SEXP R_d);
extern SEXP R_subset_dist(SEXP R_x, SEXP s);
extern SEXP R_rowSums_dist(SEXP R_x, SEXP na_rm);
extern SEXP R_row_dist(SEXP x, SEXP col);

extern SEXP R_apply_dist_matrix(SEXP p);
extern SEXP R_apply_dist_list(SEXP p);
extern SEXP R_apply_dist_binary_matrix(SEXP p);
extern SEXP R_apply_dist_data_frame(SEXP p);

static const R_CallMethodDef CallEntries[] = {
    {"R_minkowski_dist", (DL_FUNC) R_minkowski_dist, 4},
    {"R_euclidean_dist", (DL_FUNC) R_euclidean_dist, 3},
    {"R_maximum_dist",	 (DL_FUNC) R_maximum_dist,   3},
    {"R_manhattan_dist", (DL_FUNC) R_manhattan_dist, 3},
    {"R_canberra_dist",  (DL_FUNC) R_canberra_dist,  3},
    {"R_binary_dist",	 (DL_FUNC) R_binary_dist,    3},
//  {"R_matching_dist",	 (DL_FUNC) R_matching_dist,  3},
    {"R_fuzzy_dist",	 (DL_FUNC) R_fuzzy_dist,     3},
//  {"R_mutual_dist",	 (DL_FUNC) R_mutual_dist,    3},
    {"R_bjaccard",	 (DL_FUNC) R_bjaccard,       3},
    {"R_ejaccard",	 (DL_FUNC) R_ejaccard,       3},
    {"R_edice",		 (DL_FUNC) R_edice,          3},
    {"R_cosine",	 (DL_FUNC) R_cosine,         3},
    {"R_subset_dist",	 (DL_FUNC) R_subset_dist,    2},
    {"R_rowSums_dist",   (DL_FUNC) R_rowSums_dist,   2},
    {"R_row_dist",	 (DL_FUNC) R_row_dist,	     2},
    {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalEntries[] = {
    {"R_apply_dist_matrix",	   (DL_FUNC) R_apply_dist_matrix,	 -1},
    {"R_apply_dist_list",	   (DL_FUNC) R_apply_dist_list,		 -1},
    {"R_apply_dist_binary_matrix", (DL_FUNC) R_apply_dist_binary_matrix, -1},
    {"R_apply_dist_data_frame",	   (DL_FUNC) R_apply_dist_data_frame,	 -1},
    {NULL, NULL, 0}
};

void R_init_proxy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
}

