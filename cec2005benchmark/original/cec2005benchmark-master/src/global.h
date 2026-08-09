/* Global variable and function definitions */

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <float.h>

/* Global Constants */
#define INF DBL_MAX
#define EPS 1.0e-10
#define E 2.7182818284590452353602874713526625
#ifdef PI
#undef PI
#endif
#define PI 3.1415926535897932384626433832795029

/* Global variables that you are required to initialize */
int nreal; /* number of real variables */
int nfunc; /* number of basic functions */
double bound; /* required for plotting the function profiles for nreal=2 */
int density; /* density of grid points for plotting for nreal=2 */

/* Global variables being used in evaluation of various functions */
/* These are initalized in file def2.c */
double C;
double global_bias;
double *trans_x;
double *basic_f;
double *temp_x1;
double *temp_x2;
double *temp_x3;
double *temp_x4;
double *weight;
double *sigma;
double *lambda;
double *bias;
double *norm_x;
double *norm_f;
double **o;
double **g;
double ***l;

/* Auxilary function declarations */
double maximum(double, double);
double minimum(double, double);
double modulus(double*, int);
double dot(double*, double*, int);
double mean(double*, int);

/* Basic funcion declarations */
double calc_ackley(double*);
double calc_rastrigin(double*);
double calc_weierstrass(double*);
double calc_griewank(double*);
double calc_sphere(double*);
double calc_schwefel(double*);
double calc_rosenbrock(double *x);
double nc_schaffer(double, double);
double nc_rastrigin(double*);

/* Utility function declarations */
void allocate_memory();
void initialize_f1(char *extdatadir);
void initialize_f2(char *extdatadir);
void initialize_f3(char *extdatadir);
void initialize_f4(char *extdatadir);
void initialize_f5(char *extdatadir);
void initialize_f6(char *extdatadir);
void initialize_f7(char *extdatadir);
void initialize_f8(char *extdatadir);
void initialize_f9(char *extdatadir);
void initialize_f10(char *extdatadir);
void initialize_f11(char *extdatadir);
void initialize_f12(char *extdatadir);
void initialize_f13(char *extdatadir);
void initialize_f14(char *extdatadir);
void initialize_f15(char *extdatadir);
void initialize_f16(char *extdatadir);
void initialize_f17(char *extdatadir);
void initialize_f18(char *extdatadir);
void initialize_f19(char *extdatadir);
void initialize_f20(char *extdatadir);
void initialize_f21(char *extdatadir);
void initialize_f22(char *extdatadir);
void initialize_f23(char *extdatadir);
void initialize_f24_f25(char *extdatadir);
void transform(double*, int);
void transform_norm(int);
void calc_weight(double*);
void free_memory();

/* Benchmark function declaration */
double calc_benchmark_func_f1(double *);
double calc_benchmark_func_f2(double *);
double calc_benchmark_func_f3(double *);
double calc_benchmark_func_f4(double *);
double calc_benchmark_func_f5(double *);
double calc_benchmark_func_f6(double *);
double calc_benchmark_func_f7(double *);
double calc_benchmark_func_f8(double *);
double calc_benchmark_func_f9(double *);
double calc_benchmark_func_f10(double *);
double calc_benchmark_func_f11(double *);
double calc_benchmark_func_f12(double *);
double calc_benchmark_func_f13(double *);
double calc_benchmark_func_f14(double *);
double calc_benchmark_func_f15(double *);
double calc_benchmark_func_f16(double *);
double calc_benchmark_func_f17(double *);
double calc_benchmark_func_f18(double *);
double calc_benchmark_func_f19(double *);
double calc_benchmark_func_f20(double *);
double calc_benchmark_func_f21(double *);
double calc_benchmark_func_f22(double *);
double calc_benchmark_func_f23(double *);
double calc_benchmark_func_f24_f25(double *);
void calc_benchmark_norm_f15();
void calc_benchmark_norm_f16();
void calc_benchmark_norm_f17();
void calc_benchmark_norm_f18();
void calc_benchmark_norm_f19();
void calc_benchmark_norm_f20();
void calc_benchmark_norm_f21();
void calc_benchmark_norm_f22();
void calc_benchmark_norm_f23();
void calc_benchmark_norm_f24_f25();

#endif
