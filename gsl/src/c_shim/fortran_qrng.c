#include <gsl/gsl_qrng.h>
#include <stddef.h>

void *fgsl_qrng_alloc(int type_id, unsigned int dim) {
    const gsl_qrng_type *t = NULL;
    switch (type_id) {
        case 0: t = gsl_qrng_niederreiter_2; break;
        case 1: t = gsl_qrng_sobol; break;
        default: return NULL;
    }
    return (void *)gsl_qrng_alloc(t, dim);
}
void fgsl_qrng_free(void *p) { if (p) gsl_qrng_free((gsl_qrng *)p); }
void *fgsl_qrng_clone(void *p) { return (void *)gsl_qrng_clone((const gsl_qrng *)p); }
void fgsl_qrng_init(void *p) { gsl_qrng_init((gsl_qrng *)p); }
const char *fgsl_qrng_name(void *p) { return gsl_qrng_name((const gsl_qrng *)p); }
size_t fgsl_qrng_size(void *p) { return gsl_qrng_size((const gsl_qrng *)p); }
int fgsl_qrng_get(void *p, double *x) { return gsl_qrng_get((gsl_qrng *)p, x); }
