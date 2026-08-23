#include <gsl/gsl_rng.h>
#include <stddef.h>

void *fgsl_rng_alloc(int type_id) {
    const gsl_rng_type *t = NULL;
    switch (type_id) {
        case 0: t = gsl_rng_mt19937; break;
        case 1: t = gsl_rng_ranlxs0; break;
        case 2: t = gsl_rng_ranlxs1; break;
        case 3: t = gsl_rng_ranlxs2; break;
        case 4: t = gsl_rng_ranlxd1; break;
        case 5: t = gsl_rng_ranlxd2; break;
        case 6: t = gsl_rng_ranlux; break;
        case 7: t = gsl_rng_ranlux389; break;
        case 8: t = gsl_rng_cmrg; break;
        case 9: t = gsl_rng_mrg; break;
        case 10: t = gsl_rng_taus; break;
        case 11: t = gsl_rng_taus2; break;
        case 12: t = gsl_rng_gfsr4; break;
        case 13: t = gsl_rng_minstd; break;
        default: return NULL;
    }
    return (void *)gsl_rng_alloc(t);
}

void fgsl_rng_free(void *p) { if (p) gsl_rng_free((gsl_rng *)p); }
void fgsl_rng_set(void *p, unsigned long seed) { gsl_rng_set((gsl_rng *)p, seed); }
void *fgsl_rng_clone(void *p) { return (void *)gsl_rng_clone((const gsl_rng *)p); }
const char *fgsl_rng_name(void *p) { return gsl_rng_name((const gsl_rng *)p); }
unsigned long fgsl_rng_min(void *p) { return gsl_rng_min((const gsl_rng *)p); }
unsigned long fgsl_rng_max(void *p) { return gsl_rng_max((const gsl_rng *)p); }
unsigned long fgsl_rng_get(void *p) { return gsl_rng_get((const gsl_rng *)p); }
double fgsl_rng_uniform(void *p) { return gsl_rng_uniform((const gsl_rng *)p); }
double fgsl_rng_uniform_pos(void *p) { return gsl_rng_uniform_pos((const gsl_rng *)p); }
unsigned long fgsl_rng_uniform_int(void *p, unsigned long n) { return gsl_rng_uniform_int((const gsl_rng *)p, n); }
