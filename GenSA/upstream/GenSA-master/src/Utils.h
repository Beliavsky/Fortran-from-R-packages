// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vector>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
//#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <map>


using namespace std;

//#define GENSA_DBG

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define FALSE_ 0
#define TRUE_ 1

typedef std::vector<double> dVec;
typedef std::vector<double>::iterator dVecIt;
typedef std::vector<std::string> strVec;
typedef std::vector<std::string>::const_iterator strVecIt;
typedef std::map<std::string, dVec> Map;
typedef std::map<std::string, dVec>::iterator MapIt;

typedef struct opt_struct
{
    SEXP R_fn; /* function */
    SEXP R_jc; /* judge constraint function */
    SEXP R_env; /* where to evaluate the calls */
    SEXP xNames; /* names for x */
    int verbose;
} opt_struct, *OptStruct;


class Utils
{
    public:
        static double dSign(double *a, double *b);
        static double dMod(double *x, double *y);
        static double ran2(long int *idum);
        static double yyGas(long int *idum);
        static double yyGaml(double *xx);
        static void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
                double *f, double *g, double factr, double *pgtol, double *wa,
                int * iwa, char *task, int iprint, int *lsave, int *isave,
                double *dsave);
};

#endif // UTILS_H_
