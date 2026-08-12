#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>		/* INT_MAX */
#include <float.h>		/* DBL_MAX */
#include <assert.h>
#include <time.h>

/* constants used for improving readability of code */

#define COVERED       1
#define UNCOVERED     0
#define ASSIGNED      1
#define UNASSIGNED    0
#define TRUE          1
#define FALSE         0

#define MARKED        1
#define UNMARKED      0

#define REDUCE        1
#define NOREDUCE      0

typedef struct{
  int        n;            /* order of problem             */
  double   **C;            /* cost matrix		   */
  double   **c;            /* reduced cost matrix	   */
  int       *s;            /* assignment                   */
  int       *f;            /* column i is assigned to f[i] */
  int       na;            /* number of assigned items;	   */
  int     runs;            /* number of iterations	   */
  double  cost;            /* minimum cost		   */
  time_t rtime;            /* time                         */
} AP;

/* public interface */

/* constructors and destructor */
AP     *ap_create_problem(double *t, int n);
AP     *ap_create_problem_from_matrix(double **t, int n);
AP     *ap_read_problem(char *file);
void    ap_free(AP *p);

int     ap_assignment(AP *p, int *res);
int     ap_costmatrix(AP *p, double **m);
int     ap_datamatrix(AP *p, double **m);
int     ap_iterations(AP *p);
void    ap_hungarian(AP *p);
double  ap_mincost(AP *p);
void    ap_print_solution(AP *p);
void    ap_show_data(AP *p);
int     ap_size(AP *p);
int     ap_time(AP *p);

/* error reporting */
void ap_error(char *message);

/* private functions */
void    preprocess(AP *p);
void    preassign(AP *p);
int     cover(AP *p, int *ri, int *ci);
void    reduce(AP *p, int *ri, int *ci);
