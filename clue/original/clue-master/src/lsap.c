#include <R.h>
#include "assignment.h"
#include "clue.h"

void
solve_LSAP(double *c, Sint *n, Sint *p)
{
    AP *ap;
    ap = ap_create_problem(c, *n);
    ap_hungarian(ap);
    ap_assignment(ap, p);
    ap_free(ap);
}
