/* Source file for various benchmark functions */
/* Hard-coded for every function. */
/* Some redundancy is present here and there */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "sub.h"
#include "rand.h"

double calc_benchmark_func_f1(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_sphere(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f2(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_schwefel(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f3(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = 0.0;
    for (i = 0; i < nreal; i++) {
        basic_f[0] += trans_x[i] * trans_x[i] * pow(1.0e6, i / (nreal - 1.0));
    }
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f4(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_schwefel(trans_x)
            * (1.0 + 0.4 * fabs(randomnormaldeviate()));
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f5(double *x)
{
    int i, j;
    double res;
    basic_f[0] = -INF;
    for (i = 0; i < nreal; i++) {
        res = 0.0;
        for (j = 0; j < nreal; j++) {
            res += A_f5[i][j] * x[j];
        }
        res = fabs(res - B_f5[i]);
        if (basic_f[0] < res) {
            basic_f[0] = res;
        }
    }
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f6(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_rosenbrock(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f7(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_griewank(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f8(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_ackley(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f9(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f10(double *x)
{
    double res;
    transform(x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f11(double *x)
{
    int i;
    double res;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 0);
    basic_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f12(double *x)
{
    double res;
    double sum1, sum2;
    int i, j;
    basic_f[0] = 0.0;
    for (i = 0; i < nreal; i++) {
        sum1 = 0.0;
        sum2 = 0.0;
        for (j = 0; j < nreal; j++) {
            sum1 += A_f12[i][j] * sin(alpha_f12[j])
                    + B_f12[i][j] * cos(alpha_f12[j]);
            sum2 += A_f12[i][j] * sin(x[j]) + B_f12[i][j] * cos(x[j]);
        }
        basic_f[0] += pow((sum1 - sum2), 2.0);
    }
    res = basic_f[0] + bias[0];
    return (res);
}

double calc_benchmark_func_f13(double *x)
{
    int i;
    double temp;
    double res;
    transform(x, 0);
    res = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0 * pow((trans_x[i] * trans_x[i] - trans_x[i + 1]), 2.0)+1.0
                * pow((trans_x[i] - 1.0), 2.0);
        res += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0 * pow((trans_x[nreal-1] * trans_x[nreal-1] - trans_x[0]), 2.0) + 1.0 * pow((trans_x[nreal-1] - 1.0), 2.0);
    res += (temp*temp)/4000.0 - cos(temp) + 1.0 + bias[0];
    return (res);
}

double calc_benchmark_func_f14(double *x)
{
    int i;
    double temp1, temp2;
    double res;
    transform(x,0);
    res = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i], 2.0) + pow(trans_x[i+1], 2.0)))), 2.0);
        temp2 = 1.0 + 0.001 * (pow(trans_x[i], 2.0) + pow(trans_x[i + 1], 2.0));
        res += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    res += 0.5 + (temp1-0.5)/(pow(temp2,2.0)) + bias[0];
    return (res);
}

void calc_benchmark_norm_f15()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm(1);
    norm_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(2);
    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(3);
    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(4);
    norm_f[4] = calc_griewank(trans_x);
    transform_norm(5);
    norm_f[5] = calc_griewank(trans_x);
    transform_norm(6);
    norm_f[6] = calc_ackley(trans_x);
    transform_norm(7);
    norm_f[7] = calc_ackley(trans_x);
    transform_norm(8);
    norm_f[8] = calc_sphere(trans_x);
    transform_norm(9);
    norm_f[9] = calc_sphere(trans_x);
    return;
}

double calc_benchmark_func_f15(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    transform(x, 1);
    basic_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 2);
    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 3);
    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 4);
    basic_f[4] = calc_griewank(trans_x);
    transform(x, 5);
    basic_f[5] = calc_griewank(trans_x);
    transform(x, 6);
    basic_f[6] = calc_ackley(trans_x);
    transform(x, 7);
    basic_f[7] = calc_ackley(trans_x);
    transform(x, 8);
    basic_f[8] = calc_sphere(trans_x);
    transform(x, 9);
    basic_f[9] = calc_sphere(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f16()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm(1);
    norm_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(2);
    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(3);
    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(4);
    norm_f[4] = calc_griewank(trans_x);
    transform_norm(5);
    norm_f[5] = calc_griewank(trans_x);
    transform_norm(6);
    norm_f[6] = calc_ackley(trans_x);
    transform_norm(7);
    norm_f[7] = calc_ackley(trans_x);
    transform_norm(8);
    norm_f[8] = calc_sphere(trans_x);
    transform_norm(9);
    norm_f[9] = calc_sphere(trans_x);
    return;
}

double calc_benchmark_func_f16(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    transform(x, 1);
    basic_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 2);
    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 3);
    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 4);
    basic_f[4] = calc_griewank(trans_x);
    transform(x, 5);
    basic_f[5] = calc_griewank(trans_x);
    transform(x, 6);
    basic_f[6] = calc_ackley(trans_x);
    transform(x, 7);
    basic_f[7] = calc_ackley(trans_x);
    transform(x, 8);
    basic_f[8] = calc_sphere(trans_x);
    transform(x, 9);
    basic_f[9] = calc_sphere(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f17()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_rastrigin(trans_x);
    transform_norm(1);
    norm_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(2);
    norm_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(3);
    norm_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(4);
    norm_f[4] = calc_griewank(trans_x);
    transform_norm(5);
    norm_f[5] = calc_griewank(trans_x);
    transform_norm(6);
    norm_f[6] = calc_ackley(trans_x);
    transform_norm(7);
    norm_f[7] = calc_ackley(trans_x);
    transform_norm(8);
    norm_f[8] = calc_sphere(trans_x);
    transform_norm(9);
    norm_f[9] = calc_sphere(trans_x);
    return;
}

double calc_benchmark_func_f17(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_rastrigin(trans_x);
    transform(x, 1);
    basic_f[1] = calc_rastrigin(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 2);
    basic_f[2] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 3);
    basic_f[3] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 4);
    basic_f[4] = calc_griewank(trans_x);
    transform(x, 5);
    basic_f[5] = calc_griewank(trans_x);
    transform(x, 6);
    basic_f[6] = calc_ackley(trans_x);
    transform(x, 7);
    basic_f[7] = calc_ackley(trans_x);
    transform(x, 8);
    basic_f[8] = calc_sphere(trans_x);
    transform(x, 9);
    basic_f[9] = calc_sphere(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = 0.0;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    res = res * (1.0 + 0.2 * fabs(randomnormaldeviate())) + global_bias;
    return (res);
}

void calc_benchmark_norm_f18()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_ackley(trans_x);
    transform_norm(1);
    norm_f[1] = calc_ackley(trans_x);
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = calc_sphere(trans_x);
    transform_norm(5);
    norm_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f18(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_ackley(trans_x);
    transform(x, 1);
    basic_f[1] = calc_ackley(trans_x);
    transform(x, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(x, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(x, 4);
    basic_f[4] = calc_sphere(trans_x);
    transform(x, 5);
    basic_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(x, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f19()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_ackley(trans_x);
    transform_norm(1);
    norm_f[1] = calc_ackley(trans_x);
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = calc_sphere(trans_x);
    transform_norm(5);
    norm_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f19(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_ackley(trans_x);
    transform(x, 1);
    basic_f[1] = calc_ackley(trans_x);
    transform(x, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(x, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(x, 4);
    basic_f[4] = calc_sphere(trans_x);
    transform(x, 5);
    basic_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(x, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f20()
{
    int i;
    transform_norm(0);
    norm_f[0] = calc_ackley(trans_x);
    transform_norm(1);
    norm_f[1] = calc_ackley(trans_x);
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = calc_sphere(trans_x);
    transform_norm(5);
    norm_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f20(double *x)
{
    int i;
    double res;
    transform(x, 0);
    basic_f[0] = calc_ackley(trans_x);
    transform(x, 1);
    basic_f[1] = calc_ackley(trans_x);
    transform(x, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(x, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(x, 4);
    basic_f[4] = calc_sphere(trans_x);
    transform(x, 5);
    basic_f[5] = calc_sphere(trans_x);
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(x, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f21()
{
    int i;
    double temp1, temp2, temp;
    transform_norm(0);
    norm_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i], 2.0) + pow(trans_x[i+1], 2.0)))), 2.0);
        temp2 = 1.0 + 0.001 * (pow(trans_x[i], 2.0) + pow(trans_x[i+1], 2.0));
        norm_f[0] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1], 2.0) + pow(trans_x[0], 2.0)))), 2.0);
    temp2 = 1.0 + 0.001 * (pow(trans_x[nreal - 1], 2.0) + pow(trans_x[0], 2.0));
    norm_f[0] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i], 2.0) + pow(trans_x[i+1], 2.0)))), 2.0);
        temp2 = 1.0 + 0.001 * (pow(trans_x[i], 2.0) + pow(trans_x[i + 1], 2.0));
        norm_f[1] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1], 2.0) + pow(trans_x[0], 2.0)))), 2.0);
    temp2 = 1.0 + 0.001 * (pow(trans_x[nreal - 1], 2.0) + pow(trans_x[0], 2.0));
    norm_f[1] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0 * pow((trans_x[i] * trans_x[i] - trans_x[i+1]), 2.0) + 1.0 * pow((trans_x[i] - 1.0), 2.0);
        norm_f[4] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0 * pow((trans_x[nreal-1] * trans_x[nreal-1]-trans_x[0]), 2.0) + 1.0 * pow((trans_x[nreal-1] - 1.0), 2.0);
    norm_f[4] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0 * pow((trans_x[i] * trans_x[i] - trans_x[i+1]), 2.0) + 1.0 * pow((trans_x[i] - 1.0), 2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f21(double *x)
{
    int i;
    double temp1, temp2, temp;
    double res;
    transform(x, 0);
    basic_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 1);
    basic_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform(x, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(x, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(x, 4);
    basic_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(x, 5);
    basic_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(x, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f22()
{
    int i;
    double temp1, temp2, temp;
    transform_norm(0);
    norm_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm(1);
    norm_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp =
            100.0
                    * pow((trans_x[nreal - 1] * trans_x[nreal - 1]
                                  - trans_x[0]), 2.0)
                          +1.0 * pow((trans_x[nreal - 1] - 1.0), 2.0);
    norm_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f22(double *x)
{
    int i;
    double temp1, temp2, temp;
    double res;
    transform(x, 0);
    basic_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (x, 1);
    basic_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform(x, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(x, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(x, 4);
    basic_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(x, 5);
    basic_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(x, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(x, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(x, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f23()
{
    int i;
    double temp1, temp2, temp;
    transform_norm (0);
    norm_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm(2);
    norm_f[2] = calc_rastrigin(trans_x);
    transform_norm(3);
    norm_f[3] = calc_rastrigin(trans_x);
    transform_norm(4);
    norm_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform_norm(5);
    norm_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(6);
    norm_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(7);
    norm_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm(8);
    norm_f[8] = calc_griewank(trans_x);
    transform_norm(9);
    norm_f[9] = calc_griewank(trans_x);
    return;
}

double calc_benchmark_func_f23(double *x)
{
    int i;
    double temp1, temp2, temp;
    double res;
    int a;
    double b;
    for (i = 0; i < nreal; i++) {
        if (fabs(x[i] - o[0][i]) >= 0.5) {
            res = 2.0 * x[i];
            a = res;
            b = fabs(res - a);
            if (b < 0.5) {
                temp_x4[i] = a / 2.0;
            } else {
                if (res <= 0.0) {
                    temp_x4[i] = (a - 1.0) / 2.0;
                } else {
                    temp_x4[i] = (a + 1.0) / 2.0;
                }
            }
        } else {
            temp_x4[i] = x[i];
        }
    }
    transform(temp_x4, 0);
    basic_f[0] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform (temp_x4, 1);
    basic_f[1] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform(temp_x4, 2);
    basic_f[2] = calc_rastrigin(trans_x);
    transform(temp_x4, 3);
    basic_f[3] = calc_rastrigin(trans_x);
    transform(temp_x4, 4);
    basic_f[4] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    transform(temp_x4, 5);
    basic_f[5] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform(temp_x4, 6);
    basic_f[6] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(temp_x4, 7);
    basic_f[7] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform(temp_x4, 8);
    basic_f[8] = calc_griewank(trans_x);
    transform(temp_x4, 9);
    basic_f[9] = calc_griewank(trans_x);
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(temp_x4);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}

void calc_benchmark_norm_f24_f25()
{
    int i;
    double temp1, temp2, temp;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    transform_norm(0);
    norm_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);
    transform_norm (1);
    norm_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    transform_norm(2);
    norm_f[2] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        norm_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    norm_f[2] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    transform_norm(3);
    norm_f[3] = calc_ackley(trans_x);
    transform_norm(4);
    norm_f[4] = calc_rastrigin(trans_x);
    transform_norm(5);
    norm_f[5] = calc_griewank(trans_x);
    transform_norm(6);
    norm_f[6] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        norm_f[6] += nc_schaffer(trans_x[i], trans_x[i + 1]);
    }
    norm_f[6] += nc_schaffer(trans_x[nreal - 1], trans_x[0]);
    transform_norm(7);
    norm_f[7] = nc_rastrigin(trans_x);
    transform_norm(8);
    norm_f[8] = 0.0;
    for (i = 0; i < nreal; i++) {
        norm_f[8] += trans_x[i] * trans_x[i] * pow(1.0e6, i / (nreal - 1.0));
    }
    transform_norm(9);
    norm_f[9] = calc_sphere(trans_x)
            * (1.0 + 0.1 * fabs(randomnormaldeviate()));
    return;
}

double calc_benchmark_func_f24_f25(double *x)
{
    int i;
    double temp1, temp2, temp;
    double res;
    for (i = 0; i < nreal; i++) {
        norm_x[i] = 0.0;
    }
    /* First function */
    transform(x, 0);
    basic_f[0] = calc_weierstrass(trans_x) - calc_weierstrass(norm_x);

    /* Second function */
    transform(x, 1);
    basic_f[1] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp1 = pow((sin(sqrt(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0)))),2.0);
        temp2 = 1.0 + 0.001*(pow(trans_x[i],2.0)+pow(trans_x[i+1],2.0));
        basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    }
    temp1 = pow((sin(sqrt(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0)))),2.0);
    temp2 = 1.0 + 0.001*(pow(trans_x[nreal-1],2.0)+pow(trans_x[0],2.0));
    basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));

    /* Third Function */
    transform(x, 2);
    basic_f[2] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        temp = 100.0*pow((trans_x[i]*trans_x[i]-trans_x[i+1]),2.0) + 1.0*pow((trans_x[i]-1.0),2.0);
        basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    temp = 100.0*pow((trans_x[nreal-1]*trans_x[nreal-1]-trans_x[0]),2.0) + 1.0*pow((trans_x[nreal-1]-1.0),2.0);
    basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;

    transform(x, 3);
    basic_f[3] = calc_ackley(trans_x);
    transform(x, 4);
    basic_f[4] = calc_rastrigin(trans_x);
    transform(x, 5);
    basic_f[5] = calc_griewank(trans_x);

    /* Seventh Function */
    transform (x, 6);
    basic_f[6] = 0.0;
    for (i = 0; i < nreal - 1; i++) {
        basic_f[6] += nc_schaffer(trans_x[i], trans_x[i + 1]);
    }
    basic_f[6] += nc_schaffer(trans_x[nreal - 1], trans_x[0]);

    transform(x, 7);
    basic_f[7] = nc_rastrigin(trans_x);

    transform(x, 8);
    basic_f[8] = 0.0;
    for (i = 0; i < nreal; i++) {
        basic_f[8] += trans_x[i]*trans_x[i]*pow(1.0e6,i/(nreal-1.0));
    }
    transform(x, 9);
    basic_f[9] = (calc_sphere(trans_x))
            * (1.0 + 0.1 * fabs(randomnormaldeviate()));
    for (i = 0; i < nfunc; i++) {
        basic_f[i] *= C / norm_f[i];
    }
    calc_weight(x);
    res = global_bias;
    for (i = 0; i < nfunc; i++) {
        res += weight[i] * (basic_f[i] + bias[i]);
    }
    return (res);
}
