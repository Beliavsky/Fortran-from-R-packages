/* Definition of random number generation routines */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "global.h"
#include "rand.h"

int randenabled = 1;

void disablerand()
{
    randenabled = 0;
}

/* Fetch a single random number between 0.0 and 1.0 */
double randomperc()
{
    double res = 0;

    if (randenabled) {
        GetRNGstate();
        res = unif_rand();
        PutRNGstate();
    }

    return res;
}

/* Compute the noise */
double randomnormaldeviate()
{
    double res = 0;

    if (randenabled) {
        GetRNGstate();
        res = norm_rand();
        PutRNGstate();
    }

    return res;
}
