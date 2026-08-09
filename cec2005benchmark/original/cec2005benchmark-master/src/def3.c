/* Source file for custom initialization */
/* Hard-coded for every function based on the type and nature of input files */
/* At present hard-coded for D=2, 10, 30 and 50 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "global.h"
#include "sub.h"
#include "rand.h"

FILE *open_input_data(char *extdata_dir, char *name)
{
    FILE *fpt = NULL;
    char buf[PATH_MAX];
    snprintf(buf, PATH_MAX, "%s/%s", extdata_dir, name);
    fpt = fopen(buf, "r");
    if (fpt == NULL) {
        error("cannot open input file for reading\n");
    }
    return fpt;
}

void initialize_f1(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "sphere_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -450.0;
    return;
}

void initialize_f2(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "schwefel_102_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -450.0;
    return;
}

void initialize_f3(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "elliptic_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "elliptic_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "elliptic_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "elliptic_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "high_cond_elliptic_rot_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -450.0;
    return;
}

void initialize_f4(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "schwefel_102_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -450.0;
    return;
}

void initialize_f5(char *extdata_dir)
{
    int i, j;
    int index;
    FILE *fpt = NULL;
    char c;
    A_f5 = (double **) malloc(nreal * sizeof(double));
    for (i = 0; i < nreal; i++) {
        A_f5[i] = (double *) malloc(nreal * sizeof(double));
    }
    B_f5 = (double *) malloc(nreal * sizeof(double));
    fpt = open_input_data(extdata_dir, "schwefel_206_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &A_f5[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal % 4 == 0) {
        index = nreal / 4;
    } else {
        index = nreal / 4 + 1;
    }
    for (i = 0; i < index; i++) {
        o[0][i] = -100.0;
    }
    index = (3 * nreal) / 4 - 1;
    for (i = index; i < nreal; i++) {
        o[0][i] = 100.0;
    }
    for (i = 0; i < nreal; i++) {
        B_f5[i] = 0.0;
        for (j = 0; j < nreal; j++) {
            B_f5[i] += A_f5[i][j] * o[0][j];
        }
    }
    bias[0] = -310.0;
    return;
}

void initialize_f6(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "rosenbrock_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
            o[i][j] -= 1.0;
        }
    }
    fclose(fpt);
    bias[0] = 390.0;
    return;
}

void initialize_f7(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "griewank_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "griewank_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "griewank_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "griewank_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "griewank_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -180.0;
    return;
}

void initialize_f8(char *extdata_dir)
{
    int i, j;
    int index;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "ackley_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "ackley_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "ackley_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "ackley_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "ackley_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    index = nreal / 2;
    for (i = 1; i <= index; i++) {
        o[0][2 * i - 2] = -32.0;
    }
    bias[0] = -140.0;
    return;
}

void initialize_f9(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "rastrigin_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -330.0;
    return;
}

void initialize_f10(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "rastrigin_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "rastrigin_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "rastrigin_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "rastrigin_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "rastrigin_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -330.0;
    return;
}

void initialize_f11(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "weierstrass_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "weierstrass_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "weierstrass_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "weierstrass_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "weierstrass_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = 90.0;
    return;
}

void initialize_f12(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    char c;
    A_f12 = (double **) malloc(nreal * sizeof(double));
    B_f12 = (double **) malloc(nreal * sizeof(double));
    alpha_f12 = (double *) malloc(nreal * sizeof(double));
    for (i = 0; i < nreal; i++) {
        A_f12[i] = (double *) malloc(nreal * sizeof(double));
        B_f12[i] = (double *) malloc(nreal * sizeof(double));
    }
    fpt = open_input_data(extdata_dir, "schwefel_213_data.txt");
    /* Reading A */
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &A_f12[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    if (i != 100) {
        for (i = nreal; i < 100; i++) {
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    /* Reading B */
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &B_f12[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    if (i != 100) {
        for (i = nreal; i < 100; i++) {
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    /* Reading alpha */
    for (i = 0; i < nreal; i++) {
        fscanf(fpt, "%lf", &alpha_f12[i]);
    }
    fclose(fpt);
    bias[0] = -460.0;
    return;
}

void initialize_f13(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    fpt = open_input_data(extdata_dir, "EF8F2_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
            o[i][j] -= 1.0;
        }
    }
    fclose(fpt);
    bias[0] = -130.0;
    return;
}

void initialize_f14(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    if (nreal == 2) fpt = open_input_data(extdata_dir, "E_ScafferF6_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "E_ScafferF6_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "E_ScafferF6_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "E_ScafferF6_M_D50.txt");
    for (i = 0; i < nreal; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &g[i][j]);
        }
    }
    fclose(fpt);
    fpt = open_input_data(extdata_dir, "E_ScafferF6_func_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
    }
    fclose(fpt);
    bias[0] = -300.0;
    return;
}

void initialize_f15(char *extdata_dir)
{
    int i, j;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func1_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0 / 12.0;
    lambda[5] = 1.0 / 12.0;
    lambda[6] = 5.0 / 32.0;
    lambda[7] = 5.0 / 32.0;
    lambda[8] = 1.0 / 20.0;
    lambda[9] = 1.0 / 20.0;
    global_bias = 120.0;
    return;
}

void initialize_f16(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func1_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    fclose(fpt);
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0 / 12.0;
    lambda[5] = 1.0 / 12.0;
    lambda[6] = 5.0 / 32.0;
    lambda[7] = 5.0 / 32.0;
    lambda[8] = 1.0 / 20.0;
    lambda[9] = 1.0 / 20.0;
    global_bias = 120.0;
    return;
}

void initialize_f17(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func1_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func1_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    fclose(fpt);
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0 / 12.0;
    lambda[5] = 1.0 / 12.0;
    lambda[6] = 5.0 / 32.0;
    lambda[7] = 5.0 / 32.0;
    lambda[8] = 1.0 / 20.0;
    lambda[9] = 1.0 / 20.0;
    global_bias = 120.0;
    return;
}

void initialize_f18(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func2_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    for (i = 0; i < nreal; i++) {
        o[nfunc - 1][i] = 0.0;
    }
    fclose(fpt);
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0 / 16.0;
    lambda[1] = 5.0 / 32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0 / 10.0;
    lambda[5] = 1.0 / 20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 6.0;
    lambda[9] = 1.0 / 12.0;
    global_bias = 10.0;
    return;
}

void initialize_f19(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func2_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    for (i = 0; i < nreal; i++) {
        o[nfunc - 1][i] = 0.0;
    }
    fclose(fpt);
    sigma[0] = 0.1;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 0.5 / 32.0;
    lambda[1] = 5.0 / 32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0 / 10.0;
    lambda[5] = 1.0 / 20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 6.0;
    lambda[9] = 1.0 / 12.0;
    global_bias = 10.0;
    return;
}

void initialize_f20(char *extdata_dir)
{
    int i, j, k;
    int index;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func2_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    index = nreal / 2;
    for (i = 1; i <= index; i++) {
        o[0][2 * i - 1] = 5.0;
    }
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func2_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    for (i = 0; i < nreal; i++) {
        o[nfunc - 1][i] = 0.0;
    }
    fclose(fpt);
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0 / 16.0;
    lambda[1] = 5.0 / 32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0 / 10.0;
    lambda[5] = 1.0 / 20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 6.0;
    lambda[9] = 1.0 / 12.0;
    global_bias = 10.0;
    return;
}

void initialize_f21(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func3_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    fclose(fpt);
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0 / 4.0;
    lambda[1] = 1.0 / 20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 8.0;
    lambda[9] = 1.0 / 40.0;
    global_bias = 360.0;
    return;
}

void initialize_f22(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func3_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func3_HM_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func3_HM_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func3_HM_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func3_HM_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    fclose(fpt);
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0 / 4.0;
    lambda[1] = 1.0 / 20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 8.0;
    lambda[9] = 1.0 / 40.0;
    global_bias = 360.0;
    return;
}

void initialize_f23(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func3_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func3_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    fclose(fpt);
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0 / 4.0;
    lambda[1] = 1.0 / 20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0 / 8.0;
    lambda[9] = 1.0 / 40.0;
    global_bias = 360.0;
    return;
}

void initialize_f24_f25(char *extdata_dir)
{
    int i, j, k;
    FILE *fpt = NULL;
    char c;
    fpt = open_input_data(extdata_dir, "hybrid_func4_data.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            fscanf(fpt, "%lf", &o[i][j]);
        }
        do {
            fscanf(fpt, "%c", &c);
        } while (c != '\n');
    }
    fclose(fpt);
    if (nreal == 2) fpt = open_input_data(extdata_dir, "hybrid_func4_M_D2.txt");
    if (nreal == 10) fpt = open_input_data(extdata_dir, "hybrid_func4_M_D10.txt");
    if (nreal == 30) fpt = open_input_data(extdata_dir, "hybrid_func4_M_D30.txt");
    if (nreal == 50) fpt = open_input_data(extdata_dir, "hybrid_func4_M_D50.txt");
    for (i = 0; i < nfunc; i++) {
        for (j = 0; j < nreal; j++) {
            for (k = 0; k < nreal; k++) {
                fscanf(fpt, "%lf", &l[i][j][k]);
            }
            do {
                fscanf(fpt, "%c", &c);
            } while (c != '\n');
        }
    }
    for (i = 0; i < nfunc; i++) {
        sigma[i] = 2.0;
    }
    fclose(fpt);
    lambda[0] = 10.0;
    lambda[1] = 1.0 / 4.0;
    lambda[2] = 1.0;
    lambda[3] = 5.0 / 32.0;
    lambda[4] = 1.0;
    lambda[5] = 1.0 / 20.0;
    lambda[6] = 1.0 / 10.0;
    lambda[7] = 1.0;
    lambda[8] = 1.0 / 20.0;
    lambda[9] = 1.0 / 20.0;
    global_bias = 260.0;
    return;
}
