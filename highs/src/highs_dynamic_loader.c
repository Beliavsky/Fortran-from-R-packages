/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "highs_fortran_abi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
typedef HMODULE hf_lib_handle;
typedef FARPROC hf_sym_handle;
#define HF_DLOPEN(path) LoadLibraryA(path)
#define HF_DLSYM(lib, name) GetProcAddress(lib, name)
#define HF_DLCLOSE(lib) FreeLibrary(lib)
#else
#include <dlfcn.h>
typedef void* hf_lib_handle;
typedef void* hf_sym_handle;
#define HF_DLOPEN(path) dlopen(path, RTLD_NOW | RTLD_LOCAL)
#define HF_DLSYM(lib, name) dlsym(lib, name)
#define HF_DLCLOSE(lib) dlclose(lib)
#endif

static hf_lib_handle g_lib = 0;
static char g_error[1024] = "backend not loaded";

typedef void* (*p_create)(void);
typedef void (*p_destroy)(void*);
typedef int (*p_pass_model)(void*, int, int, int, int, int, double,
  const double*, const double*, const double*, const double*, const double*,
  const int*, const int*, const double*, const int*);
typedef int (*p_pass_hessian)(void*, int, int, int, const int*, const int*, const double*);
typedef int (*p_run)(void*);
typedef int (*p_get_int)(void*);
typedef double (*p_get_double)(void*);
typedef int (*p_get_solution)(void*, double*, double*, double*, double*, int*, int*);
typedef int (*p_get_info)(void*, hf_info_c*);
typedef int (*p_copy_string)(void*, char*, int);
typedef int (*p_copy_global_string)(char*, int);
typedef int (*p_set_bool)(void*, const char*, int);
typedef int (*p_set_int)(void*, const char*, int);
typedef int (*p_set_double)(void*, const char*, double);
typedef int (*p_set_string)(void*, const char*, const char*);
typedef int (*p_simple)(void*);
typedef int (*p_change_indexed_double)(void*, int, const int*, const double*);
typedef int (*p_change_indexed_bounds)(void*, int, const int*, const double*, const double*);
typedef int (*p_change_coeff)(void*, int, int, double);
typedef int (*p_change_integrality)(void*, int, const int*, const int*);
typedef int (*p_change_sense)(void*, int);
typedef int (*p_change_offset)(void*, double);
typedef int (*p_file)(void*, const char*);
typedef int (*p_set_solution)(void*, const double*, const double*, const double*, const double*, int, int);
typedef int (*p_get_basis)(void*, int*, int*, int*);
typedef int (*p_set_basis)(void*, const int*, const int*);
typedef int (*p_get_ray)(void*, double*, int*);

static p_create f_create;
static p_destroy f_destroy;
static p_pass_model f_pass_model;
static p_pass_hessian f_pass_hessian;
static p_run f_run;
static p_get_int f_model_status, f_num_col, f_num_row, f_get_sense;
static p_get_double f_objective_value, f_infinity;
static p_get_solution f_get_solution;
static p_get_info f_get_info;
static p_copy_string f_status_message;
static p_copy_global_string f_version;
static p_set_bool f_set_bool;
static p_set_int f_set_int;
static p_set_double f_set_double;
static p_set_string f_set_string;
static p_simple f_reset_options, f_clear_model, f_clear_solver, f_presolve, f_clear_basis;
static p_change_indexed_double f_change_cost;
static p_change_indexed_bounds f_change_col_bounds, f_change_row_bounds;
static p_change_coeff f_change_coeff;
static p_change_integrality f_change_integrality;
static p_change_sense f_change_sense;
static p_change_offset f_change_offset;
static p_file f_write_model, f_read_model;
static p_set_solution f_set_solution;
static p_get_basis f_get_basis;
static p_set_basis f_set_basis;
static p_get_ray f_get_dual_ray, f_get_primal_ray;

static void set_error(const char* s) {
  if (!s) s = "unknown loader error";
  snprintf(g_error, sizeof(g_error), "%s", s);
}

#ifdef _WIN32
static void set_windows_error(const char* prefix) {
  DWORD code = GetLastError();
  char* msg = NULL;
  FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                 FORMAT_MESSAGE_IGNORE_INSERTS, NULL, code, 0, (LPSTR)&msg, 0, NULL);
  snprintf(g_error, sizeof(g_error), "%s: %s", prefix, msg ? msg : "unknown Windows error");
  if (msg) LocalFree(msg);
}
#else
static void set_dl_error(const char* prefix) {
  const char* e = dlerror();
  snprintf(g_error, sizeof(g_error), "%s: %s", prefix, e ? e : "unknown dlopen error");
}
#endif

#define RESOLVE(name, var, type) do { \
  hf_sym_handle sym = HF_DLSYM(g_lib, name); \
  (void)sizeof(type); \
  if (!sym) { set_error("backend is missing required symbol " name); goto fail; } \
  memcpy(&(var), &sym, sizeof(var)); \
} while (0)

static int resolve_all(void) {
  RESOLVE("hf_backend_create", f_create, p_create);
  RESOLVE("hf_backend_destroy", f_destroy, p_destroy);
  RESOLVE("hf_backend_pass_model", f_pass_model, p_pass_model);
  RESOLVE("hf_backend_pass_hessian", f_pass_hessian, p_pass_hessian);
  RESOLVE("hf_backend_run", f_run, p_run);
  RESOLVE("hf_backend_model_status", f_model_status, p_get_int);
  RESOLVE("hf_backend_objective_value", f_objective_value, p_get_double);
  RESOLVE("hf_backend_get_solution", f_get_solution, p_get_solution);
  RESOLVE("hf_backend_get_info", f_get_info, p_get_info);
  RESOLVE("hf_backend_status_message", f_status_message, p_copy_string);
  RESOLVE("hf_backend_version", f_version, p_copy_global_string);
  RESOLVE("hf_backend_infinity", f_infinity, p_get_double);
  RESOLVE("hf_backend_num_col", f_num_col, p_get_int);
  RESOLVE("hf_backend_num_row", f_num_row, p_get_int);
  RESOLVE("hf_backend_get_sense", f_get_sense, p_get_int);
  RESOLVE("hf_backend_set_bool_option", f_set_bool, p_set_bool);
  RESOLVE("hf_backend_set_int_option", f_set_int, p_set_int);
  RESOLVE("hf_backend_set_double_option", f_set_double, p_set_double);
  RESOLVE("hf_backend_set_string_option", f_set_string, p_set_string);
  RESOLVE("hf_backend_reset_options", f_reset_options, p_simple);
  RESOLVE("hf_backend_clear_model", f_clear_model, p_simple);
  RESOLVE("hf_backend_clear_solver", f_clear_solver, p_simple);
  RESOLVE("hf_backend_presolve", f_presolve, p_simple);
  RESOLVE("hf_backend_change_cost", f_change_cost, p_change_indexed_double);
  RESOLVE("hf_backend_change_col_bounds", f_change_col_bounds, p_change_indexed_bounds);
  RESOLVE("hf_backend_change_row_bounds", f_change_row_bounds, p_change_indexed_bounds);
  RESOLVE("hf_backend_change_coeff", f_change_coeff, p_change_coeff);
  RESOLVE("hf_backend_change_integrality", f_change_integrality, p_change_integrality);
  RESOLVE("hf_backend_change_sense", f_change_sense, p_change_sense);
  RESOLVE("hf_backend_change_offset", f_change_offset, p_change_offset);
  RESOLVE("hf_backend_write_model", f_write_model, p_file);
  RESOLVE("hf_backend_read_model", f_read_model, p_file);
  RESOLVE("hf_backend_set_solution", f_set_solution, p_set_solution);
  RESOLVE("hf_backend_get_basis", f_get_basis, p_get_basis);
  RESOLVE("hf_backend_set_basis", f_set_basis, p_set_basis);
  RESOLVE("hf_backend_clear_basis", f_clear_basis, p_simple);
  RESOLVE("hf_backend_get_dual_ray", f_get_dual_ray, p_get_ray);
  RESOLVE("hf_backend_get_primal_ray", f_get_primal_ray, p_get_ray);
  return 1;
fail:
  return 0;
}

static int try_path(const char* path) {
  if (!path || !*path) return 0;
  g_lib = HF_DLOPEN(path);
  if (!g_lib) {
#ifdef _WIN32
    set_windows_error(path);
#else
    set_dl_error(path);
#endif
    return 0;
  }
  if (!resolve_all()) {
    HF_DLCLOSE(g_lib);
    g_lib = 0;
    return 0;
  }
  set_error("");
  return 1;
}

int hf_api_load_backend(const char* explicit_path) {
  const char* env;
  const char* candidates[] = {
#ifdef _WIN32
    "highs_fortran_bridge.dll",
    "backend\\bin\\highs_fortran_bridge.dll",
    ".\\backend\\bin\\highs_fortran_bridge.dll",
#else
    "libhighs_fortran_bridge.so",
    "backend/lib/libhighs_fortran_bridge.so",
    "./backend/lib/libhighs_fortran_bridge.so",
    "backend/bin/libhighs_fortran_bridge.so",
#endif
    NULL
  };
  int i;
  if (g_lib) return 1;
  if (try_path(explicit_path)) return 1;
  env = getenv("HIGHS_FORTRAN_BRIDGE");
  if (try_path(env)) return 1;
  for (i = 0; candidates[i]; ++i) if (try_path(candidates[i])) return 1;
  return 0;
}

int hf_api_backend_available(void) { return g_lib ? 1 : hf_api_load_backend(NULL); }
void hf_api_unload_backend(void) { if (g_lib) { HF_DLCLOSE(g_lib); g_lib = 0; } }
int hf_api_last_error(char* buffer, int n) {
  if (buffer && n > 0) { snprintf(buffer, (size_t)n, "%s", g_error); buffer[n-1] = '\0'; }
  return (int)strlen(g_error);
}

#define NEED_BACKEND_RET(ret) do { if (!hf_api_backend_available()) return (ret); } while (0)
void* hf_api_create(void) { NEED_BACKEND_RET(NULL); return f_create(); }
void hf_api_destroy(void* p) { if (g_lib && p) f_destroy(p); }
int hf_api_pass_model(void* p, int nc, int nr, int nz, int fmt, int sense, double off,
 const double* c, const double* cl, const double* cu, const double* rl, const double* ru,
 const int* st, const int* ix, const double* av, const int* integ) {
 NEED_BACKEND_RET(-100); return f_pass_model(p,nc,nr,nz,fmt,sense,off,c,cl,cu,rl,ru,st,ix,av,integ); }
int hf_api_pass_hessian(void* p,int n,int nz,int fmt,const int* st,const int* ix,const double* v) { NEED_BACKEND_RET(-100); return f_pass_hessian(p,n,nz,fmt,st,ix,v); }
int hf_api_run(void* p) { NEED_BACKEND_RET(-100); return f_run(p); }
int hf_api_model_status(void* p) { NEED_BACKEND_RET(-100); return f_model_status(p); }
double hf_api_objective_value(void* p) { NEED_BACKEND_RET(0.0); return f_objective_value(p); }
int hf_api_get_solution(void* p,double* cv,double* cd,double* rv,double* rd,int* vv,int* dv) { NEED_BACKEND_RET(-100); return f_get_solution(p,cv,cd,rv,rd,vv,dv); }
int hf_api_get_info(void* p,hf_info_c* i) { NEED_BACKEND_RET(-100); return f_get_info(p,i); }
int hf_api_status_message(void* p,char* b,int n) { NEED_BACKEND_RET(-100); return f_status_message(p,b,n); }
int hf_api_version(char* b,int n) { NEED_BACKEND_RET(-100); return f_version(b,n); }
double hf_api_infinity(void* p) { NEED_BACKEND_RET(1e30); return f_infinity(p); }
int hf_api_num_col(void* p) { NEED_BACKEND_RET(-100); return f_num_col(p); }
int hf_api_num_row(void* p) { NEED_BACKEND_RET(-100); return f_num_row(p); }
int hf_api_get_sense(void* p) { NEED_BACKEND_RET(-100); return f_get_sense(p); }
int hf_api_set_bool_option(void* p,const char* n,int v) { NEED_BACKEND_RET(-100); return f_set_bool(p,n,v); }
int hf_api_set_int_option(void* p,const char* n,int v) { NEED_BACKEND_RET(-100); return f_set_int(p,n,v); }
int hf_api_set_double_option(void* p,const char* n,double v) { NEED_BACKEND_RET(-100); return f_set_double(p,n,v); }
int hf_api_set_string_option(void* p,const char* n,const char* v) { NEED_BACKEND_RET(-100); return f_set_string(p,n,v); }
int hf_api_reset_options(void* p) { NEED_BACKEND_RET(-100); return f_reset_options(p); }
int hf_api_clear_model(void* p) { NEED_BACKEND_RET(-100); return f_clear_model(p); }
int hf_api_clear_solver(void* p) { NEED_BACKEND_RET(-100); return f_clear_solver(p); }
int hf_api_presolve(void* p) { NEED_BACKEND_RET(-100); return f_presolve(p); }
int hf_api_change_cost(void* p,int n,const int* ix,const double* v) { NEED_BACKEND_RET(-100); return f_change_cost(p,n,ix,v); }
int hf_api_change_col_bounds(void* p,int n,const int* ix,const double* l,const double* u) { NEED_BACKEND_RET(-100); return f_change_col_bounds(p,n,ix,l,u); }
int hf_api_change_row_bounds(void* p,int n,const int* ix,const double* l,const double* u) { NEED_BACKEND_RET(-100); return f_change_row_bounds(p,n,ix,l,u); }
int hf_api_change_coeff(void* p,int r,int c,double v) { NEED_BACKEND_RET(-100); return f_change_coeff(p,r,c,v); }
int hf_api_change_integrality(void* p,int n,const int* ix,const int* v) { NEED_BACKEND_RET(-100); return f_change_integrality(p,n,ix,v); }
int hf_api_change_sense(void* p,int s) { NEED_BACKEND_RET(-100); return f_change_sense(p,s); }
int hf_api_change_offset(void* p,double v) { NEED_BACKEND_RET(-100); return f_change_offset(p,v); }
int hf_api_write_model(void* p,const char* f) { NEED_BACKEND_RET(-100); return f_write_model(p,f); }
int hf_api_read_model(void* p,const char* f) { NEED_BACKEND_RET(-100); return f_read_model(p,f); }
int hf_api_set_solution(void* p,const double* cv,const double* cd,const double* rv,const double* rd,int vv,int dv) { NEED_BACKEND_RET(-100); return f_set_solution(p,cv,cd,rv,rd,vv,dv); }
int hf_api_get_basis(void* p,int* cs,int* rs,int* valid) { NEED_BACKEND_RET(-100); return f_get_basis(p,cs,rs,valid); }
int hf_api_set_basis(void* p,const int* cs,const int* rs) { NEED_BACKEND_RET(-100); return f_set_basis(p,cs,rs); }
int hf_api_clear_basis(void* p) { NEED_BACKEND_RET(-100); return f_clear_basis(p); }
int hf_api_get_dual_ray(void* p,double* v,int* has) { NEED_BACKEND_RET(-100); return f_get_dual_ray(p,v,has); }
int hf_api_get_primal_ray(void* p,double* v,int* has) { NEED_BACKEND_RET(-100); return f_get_primal_ray(p,v,has); }
