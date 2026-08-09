/* Based on the main.c example file */

#include <stdio.h>
#include <stdlib.h>

#include <R.h>

#include "global.h"
#include "sub.h"
#include "rand.h"

void
cec2005benchmark(char **extdatadir, int *i, double *X, int *row, int *col, double *f)
{
    int r, c;
    double *x;

    /* Assign nfunc and nreal in the begining */
    nfunc = 10;
    nreal = *col;

    /* nreal and nfunc need to be initialized before calling these routines */
    /* Routine to allocate memory to global variables */
    allocate_memory();

    /* Routine the initalize global variables */
    /* For test problems 15 to 25, we need to calculate a normalizing quantity */
    switch (*i) {
        case 1: initialize_f1(*extdatadir); break;
        case 2: initialize_f2(*extdatadir); break;
        case 3: initialize_f3(*extdatadir); break;
        case 4: initialize_f4(*extdatadir); break;
        case 5: initialize_f5(*extdatadir); break;
        case 6: initialize_f6(*extdatadir); break;
        case 7: initialize_f7(*extdatadir); break;
        case 8: initialize_f8(*extdatadir); break;
        case 9: initialize_f9(*extdatadir); break;
        case 10: initialize_f10(*extdatadir); break;
        case 11: initialize_f11(*extdatadir); break;
        case 12: initialize_f12(*extdatadir); break;
        case 13: initialize_f13(*extdatadir); break;
        case 14: initialize_f14(*extdatadir); break;
        case 15: initialize_f15(*extdatadir); calc_benchmark_norm_f15(); break;
        case 16: initialize_f16(*extdatadir); calc_benchmark_norm_f16(); break;
        case 17: initialize_f17(*extdatadir); calc_benchmark_norm_f17(); break;
        case 18: initialize_f18(*extdatadir); calc_benchmark_norm_f18(); break;
        case 19: initialize_f19(*extdatadir); calc_benchmark_norm_f19(); break;
        case 20: initialize_f20(*extdatadir); calc_benchmark_norm_f20(); break;
        case 21: initialize_f21(*extdatadir); calc_benchmark_norm_f21(); break;
        case 22: initialize_f22(*extdatadir); calc_benchmark_norm_f22(); break;
        case 23: initialize_f23(*extdatadir); calc_benchmark_norm_f23(); break;
        case 24: initialize_f24_f25(*extdatadir); calc_benchmark_norm_f24_f25(); break;
        case 25: initialize_f24_f25(*extdatadir); calc_benchmark_norm_f24_f25(); break;
        default: break;
    }

    /* Variable vector */
    x = (double *) malloc(nreal * sizeof(double));

    for (r = 0; r < *row; r++) {
        R_CheckUserInterrupt();

        for (c = 0; c < *col; c++) {
            x[c] = X[r + *row * c];
        }

        switch (*i) {
            case 1: f[r] = calc_benchmark_func_f1(x); break;
            case 2: f[r] = calc_benchmark_func_f2(x); break;
            case 3: f[r] = calc_benchmark_func_f3(x); break;
            case 4: f[r] = calc_benchmark_func_f4(x); break;
            case 5: f[r] = calc_benchmark_func_f5(x); break;
            case 6: f[r] = calc_benchmark_func_f6(x); break;
            case 7: f[r] = calc_benchmark_func_f7(x); break;
            case 8: f[r] = calc_benchmark_func_f8(x); break;
            case 9: f[r] = calc_benchmark_func_f9(x); break;
            case 10: f[r] = calc_benchmark_func_f10(x); break;
            case 11: f[r] = calc_benchmark_func_f11(x); break;
            case 12: f[r] = calc_benchmark_func_f12(x); break;
            case 13: f[r] = calc_benchmark_func_f13(x); break;
            case 14: f[r] = calc_benchmark_func_f14(x); break;
            case 15: f[r] = calc_benchmark_func_f15(x); break;
            case 16: f[r] = calc_benchmark_func_f16(x); break;
            case 17: f[r] = calc_benchmark_func_f17(x); break;
            case 18: f[r] = calc_benchmark_func_f18(x); break;
            case 19: f[r] = calc_benchmark_func_f19(x); break;
            case 20: f[r] = calc_benchmark_func_f20(x); break;
            case 21: f[r] = calc_benchmark_func_f21(x); break;
            case 22: f[r] = calc_benchmark_func_f22(x); break;
            case 23: f[r] = calc_benchmark_func_f23(x); break;
            case 24: f[r] = calc_benchmark_func_f24_f25(x); break;
            case 25: f[r] = calc_benchmark_func_f24_f25(x); break;
            default: break;
        }
    }

    /* Routine to free the memory allocated at run time */
    free_memory();
    free(x);
}
