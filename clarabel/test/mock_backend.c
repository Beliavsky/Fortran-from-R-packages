/* SPDX-License-Identifier: Apache-2.0
 * Test-only ABI backend.  It solves unconstrained/equality-constrained dense QPs.
 * It is not installed as the production Clarabel backend.
 */
#include "clarabel_bridge.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    size_t n, m, p_nnz, a_nnz;
    size_t *pp, *pi, *ap, *ai;
    double *px, *ax, *q, *b;
    int all_zero_cones;
    clarabel_settings_c settings;
} mock_solver;

static void err(char *buf, size_t cap, const char *msg) {
    if (!buf || cap == 0) return;
    snprintf(buf, cap, "%s", msg ? msg : "error");
}

void clarabel_settings_default(clarabel_settings_c *s) {
    memset(s, 0, sizeof(*s));
    s->max_iter = 200; s->time_limit = HUGE_VAL; s->verbose = 1;
    s->max_step_fraction = 0.99;
    s->tol_gap_abs = s->tol_gap_rel = s->tol_feas = 1e-8;
    s->tol_infeas_abs = s->tol_infeas_rel = 1e-8; s->tol_ktratio = 1e-6;
    s->reduced_tol_gap_abs = s->reduced_tol_gap_rel = 5e-5;
    s->reduced_tol_feas = 1e-4;
    s->reduced_tol_infeas_abs = s->reduced_tol_infeas_rel = 5e-5;
    s->reduced_tol_ktratio = 1e-4;
    s->equilibrate_enable = 1; s->equilibrate_max_iter = 10;
    s->equilibrate_min_scaling = 1e-4; s->equilibrate_max_scaling = 1e4;
    s->linesearch_backtrack_step = 0.8; s->min_switch_step_length = 0.1;
    s->min_terminate_step_length = 1e-4; s->direct_kkt_solver = 1;
    s->direct_solve_method = 1; s->static_regularization_enable = 1;
    s->static_regularization_constant = 1e-8;
    s->static_regularization_proportional = 4.930380657631324e-32;
    s->dynamic_regularization_enable = 1; s->dynamic_regularization_eps = 1e-13;
    s->dynamic_regularization_delta = 2e-7; s->iterative_refinement_enable = 1;
    s->iterative_refinement_reltol = 1e-13; s->iterative_refinement_abstol = 1e-12;
    s->iterative_refinement_max_iter = 10; s->iterative_refinement_stop_ratio = 5.0;
    s->presolve_enable = 1; s->chordal_decomposition_compact = 0;
    s->chordal_decomposition_complete_dual = 0; s->chordal_decomposition_merge_method = 0;
}

static void *dup_bytes(const void *p, size_t n) {
    if (n == 0) return NULL;
    void *q = malloc(n); if (q) memcpy(q, p, n); return q;
}

int32_t clarabel_solver_create(const clarabel_csc_c *P, const double *q, size_t q_len,
                               const clarabel_csc_c *A, const double *b, size_t b_len,
                               const clarabel_cone_c *cones, size_t ncones,
                               const clarabel_settings_c *settings, void **out,
                               char *error, size_t error_capacity) {
    size_t i;
    mock_solver *s;
    if (!P || !A || !settings || !out || q_len != P->ncols || b_len != A->nrows) {
        err(error, error_capacity, "bad mock-backend dimensions"); return -1;
    }
    s = calloc(1, sizeof(*s)); if (!s) return -1;
    s->n = P->ncols; s->m = A->nrows; s->p_nnz = P->nnz; s->a_nnz = A->nnz;
    s->pp = dup_bytes(P->colptr, (P->ncols + 1) * sizeof(size_t));
    s->pi = dup_bytes(P->rowind, P->nnz * sizeof(size_t));
    s->px = dup_bytes(P->values, P->nnz * sizeof(double));
    s->ap = dup_bytes(A->colptr, (A->ncols + 1) * sizeof(size_t));
    s->ai = dup_bytes(A->rowind, A->nnz * sizeof(size_t));
    s->ax = dup_bytes(A->values, A->nnz * sizeof(double));
    s->q = dup_bytes(q, q_len * sizeof(double));
    s->b = dup_bytes(b, b_len * sizeof(double));
    s->settings = *settings;
    s->all_zero_cones = 1;
    for (i = 0; i < ncones; ++i) if (cones[i].tag != 0) s->all_zero_cones = 0;
    *out = s; return 0;
}

static int dense_solve(double *a, double *rhs, size_t n) {
    size_t i,j,k,p; double best,tmp,f;
    for (k=0;k<n;k++) {
        p=k; best=fabs(a[k*n+k]);
        for (i=k+1;i<n;i++) if (fabs(a[i*n+k])>best) {best=fabs(a[i*n+k]);p=i;}
        if (best < 1e-12) return -1;
        if (p!=k) { for(j=k;j<n;j++){tmp=a[k*n+j];a[k*n+j]=a[p*n+j];a[p*n+j]=tmp;} tmp=rhs[k];rhs[k]=rhs[p];rhs[p]=tmp; }
        for(i=k+1;i<n;i++){f=a[i*n+k]/a[k*n+k];a[i*n+k]=0;for(j=k+1;j<n;j++)a[i*n+j]-=f*a[k*n+j];rhs[i]-=f*rhs[k];}
    }
    for (i=n;i-- > 0;) { for(j=i+1;j<n;j++)rhs[i]-=a[i*n+j]*rhs[j]; rhs[i]/=a[i*n+i]; }
    return 0;
}

static void dense_matrices(const mock_solver *s, double *P, double *A) {
    size_t j,k,i;
    memset(P,0,s->n*s->n*sizeof(double)); memset(A,0,s->m*s->n*sizeof(double));
    for(j=0;j<s->n;j++) for(k=s->pp[j];k<s->pp[j+1];k++){i=s->pi[k];P[i*s->n+j]=s->px[k];P[j*s->n+i]=s->px[k];}
    for(j=0;j<s->n;j++) for(k=s->ap[j];k<s->ap[j+1];k++){i=s->ai[k];A[i*s->n+j]=s->ax[k];}
}

int32_t clarabel_solver_solve(void *ptr, double *x, size_t x_len,
                              double *z, size_t z_len, double *slack, size_t s_len,
                              clarabel_result_c *r, char *error, size_t error_capacity) {
    mock_solver *s=ptr; size_t nsys,i,j; double *M,*rhs,*P,*A,obj=0;
    if(!s||!r||x_len!=s->n||z_len!=s->m||s_len!=s->m){err(error,error_capacity,"mock output size mismatch");return -1;}
    P=calloc(s->n*s->n,sizeof(double)); A=calloc(s->m*s->n,sizeof(double)); dense_matrices(s,P,A);
    nsys=s->n+(s->all_zero_cones?s->m:0); M=calloc(nsys*nsys,sizeof(double)); rhs=calloc(nsys,sizeof(double));
    for(i=0;i<s->n;i++){rhs[i]=-s->q[i];for(j=0;j<s->n;j++)M[i*nsys+j]=P[i*s->n+j];}
    if(s->all_zero_cones) for(i=0;i<s->m;i++){rhs[s->n+i]=s->b[i];for(j=0;j<s->n;j++){M[j*nsys+s->n+i]=A[i*s->n+j];M[(s->n+i)*nsys+j]=A[i*s->n+j];}}
    memset(r,0,sizeof(*r)); memset(z,0,s->m*sizeof(double)); memset(slack,0,s->m*sizeof(double));
    if(dense_solve(M,rhs,nsys)!=0){r->status=3;r->obj_val=NAN;r->obj_val_dual=NAN;free(P);free(A);free(M);free(rhs);return 0;}
    for (i = 0; i < s->n; ++i) x[i] = rhs[i];
    if (s->all_zero_cones) {
        for (i = 0; i < s->m; ++i) z[i] = rhs[s->n + i];
    }
    for(i=0;i<s->m;i++){slack[i]=s->b[i];for(j=0;j<s->n;j++)slack[i]-=A[i*s->n+j]*x[j];}
    for(i=0;i<s->n;i++){obj+=s->q[i]*x[i];for(j=0;j<s->n;j++)obj+=0.5*x[i]*P[i*s->n+j]*x[j];}
    r->status=1;r->iterations=1;r->obj_val=obj;r->obj_val_dual=obj;r->cost_primal=obj;r->cost_dual=obj;r->linear_solver_threads=1;
    free(P);free(A);free(M);free(rhs);return 0;
}

int32_t clarabel_solver_update(void *ptr, const double *pv,size_t pn,const double *av,size_t an,
                               const double *q,size_t qn,const double *b,size_t bn,
                               char *error,size_t cap){mock_solver*s=ptr;if(!s)return-1;
    if(!clarabel_solver_is_update_allowed(s)){err(error,cap,"updates disabled");return-1;}
    if(pn){if(pn!=s->p_nnz)return-1;memcpy(s->px,pv,pn*sizeof(double));}
    if(an){if(an!=s->a_nnz)return-1;memcpy(s->ax,av,an*sizeof(double));}
    if(qn){if(qn!=s->n)return-1;memcpy(s->q,q,qn*sizeof(double));}
    if(bn){if(bn!=s->m)return-1;memcpy(s->b,b,bn*sizeof(double));}return 0;}

int32_t clarabel_solver_is_update_allowed(const void *ptr){const mock_solver*s=ptr;if(!s)return 0;return !(s->settings.presolve_enable||s->settings.chordal_decomposition_enable||s->settings.input_sparse_dropzeros);}
void clarabel_solver_free(void *ptr){mock_solver*s=ptr;if(!s)return;free(s->pp);free(s->pi);free(s->px);free(s->ap);free(s->ai);free(s->ax);free(s->q);free(s->b);free(s);}
