#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <Rinternals.h>

void tvGarch1ab(int *nsample, double *a, double *b, double *C, double *P, double *X, double *F, double *K, double *y) {
    //    double a=*A;
    //    double b=*B;
        double c = *C;
    int n = *nsample;
    for (int i = 0; i < n; i++) {
        F[i] = P[i] + 1;//nuevo 
        K[i] = ((a[i] + b[i]) * P[i] + a[i]) / F[i];
        if (i < (n - 1)) {
            if (y[i] == 100000000) {
                P[i + 1] = (a[i] + b[i]) * P[i] * (a[i] + b[i]) + a[i] * (a[i]);
                X[i + 1] = (a[i] + b[i]) * X[i];
            }
            else {
                P[i + 1] = (a[i] + b[i]) * P[i] * (a[i] + b[i] - K[i]) + a[i] * (a[i] - K[i]);
                X[i + 1] = (a[i] + b[i] - K[i]) * X[i] + K[i] * y[i];
            }
            if (X[i + 1] + (c / (1 - a[i] - b[i])) < 0) {
                X[i + 1] = 0;
                //warning("Estimada varianza negativa")                               
            }
        }
    }
}
void tvGarch1ab2(int *nsample, double *a, double *b, double *c, double *P, double *X, double *F, double *K, double *y) {
    //    double a=*A;
    //    double b=*B;
    //    double c=*C;
    int n = *nsample;
    for (int i = 0; i < n; i++) {
        F[i] = P[i] + 1;//nuevo 
        K[i] = ((a[i] + b[i]) * P[i] + a[i]) / F[i];
        if (i < (n - 1)) {
            if (y[i] == 100000000) {
                P[i + 1] = (a[i] + b[i]) * P[i] * (a[i] + b[i]) + a[i] * (a[i]);
                X[i + 1] = (a[i] + b[i]) * X[i];
            } 
            else {
            P[i + 1] = (a[i] + b[i]) * P[i] * (a[i] + b[i] - K[i]) + a[i] * (a[i] - K[i]);
            X[i + 1] = (a[i] + b[i] - K[i]) * X[i] + K[i] * y[i];
            }
            if (X[i + 1] + (c[i] / (1 - a[i] - b[i])) < 0) {
                X[i + 1] = 0;
            //warning("Estimada varianza negativa")                               
            }
        }
    }
}
static const R_CMethodDef CEntries[] = {
    {"tvGarch1ab", (DL_FUNC)&tvGarch1ab, 9},
    {"tvGarch1ab2", (DL_FUNC)&tvGarch1ab2, 9},
    {NULL, NULL, 0}
};
void R_init_tvGarchKF(DllInfo* dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
