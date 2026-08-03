/* SPDX-License-Identifier: Apache-2.0 */
#include "osqp_fortran_abi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
typedef HMODULE of_lib_handle;
typedef FARPROC of_sym_handle;
#define OF_DLOPEN(path) LoadLibraryA(path)
#define OF_DLSYM(lib, name) GetProcAddress(lib, name)
#define OF_DLCLOSE(lib) FreeLibrary(lib)
#else
#include <dlfcn.h>
typedef void* of_lib_handle;
typedef void* of_sym_handle;
#define OF_DLOPEN(path) dlopen(path, RTLD_NOW | RTLD_LOCAL)
#define OF_DLSYM(lib, name) dlsym(lib, name)
#define OF_DLCLOSE(lib) dlclose(lib)
#endif

static of_lib_handle g_lib = 0;
static char g_error[1024] = "backend not loaded";

typedef void* (*p_create)(int,int,int,const int*,const int*,const double*,int,const int*,const int*,const double*,const double*,const double*,const double*,const of_settings_c*,int*);
typedef void (*p_destroy)(void*);
typedef int (*p_solve)(void*);
typedef int (*p_get_dimensions)(void*,int*,int*);
typedef int (*p_get_solution)(void*,double*,double*,double*,double*,of_info_c*);
typedef int (*p_copy_handle_string)(void*,char*,int);
typedef int (*p_copy_global_string)(char*,int);
typedef int (*p_capabilities)(void);
typedef double (*p_infinity)(void);
typedef int (*p_error_message)(int,char*,int);
typedef int (*p_get_settings)(void*,of_settings_c*);
typedef int (*p_update_vec)(void*,const double*,int,const double*,int,const double*,int);
typedef int (*p_update_mat)(void*,const double*,const int*,int,int,const double*,const int*,int,int);
typedef int (*p_warm_start)(void*,const double*,int,const double*,int);
typedef int (*p_simple)(void*);
typedef int (*p_update_settings)(void*,const of_settings_c*);

static p_create f_create;
static p_destroy f_destroy;
static p_solve f_solve;
static p_get_dimensions f_get_dimensions;
static p_get_solution f_get_solution;
static p_copy_handle_string f_status_message;
static p_copy_global_string f_version;
static p_capabilities f_capabilities;
static p_infinity f_infinity;
static p_error_message f_error_message;
static p_get_settings f_get_settings;
static p_update_vec f_update_vec;
static p_update_mat f_update_mat;
static p_warm_start f_warm_start;
static p_simple f_cold_start;
static p_update_settings f_update_settings;

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

#define RESOLVE(name, var) do { \
  of_sym_handle sym = OF_DLSYM(g_lib, name); \
  if (!sym) { set_error("backend is missing required symbol " name); goto fail; } \
  memcpy(&(var), &sym, sizeof(var)); \
} while (0)

static int resolve_all(void) {
  RESOLVE("of_backend_create", f_create);
  RESOLVE("of_backend_destroy", f_destroy);
  RESOLVE("of_backend_solve", f_solve);
  RESOLVE("of_backend_get_dimensions", f_get_dimensions);
  RESOLVE("of_backend_get_solution", f_get_solution);
  RESOLVE("of_backend_status_message", f_status_message);
  RESOLVE("of_backend_version", f_version);
  RESOLVE("of_backend_capabilities", f_capabilities);
  RESOLVE("of_backend_infinity", f_infinity);
  RESOLVE("of_backend_error_message", f_error_message);
  RESOLVE("of_backend_get_settings", f_get_settings);
  RESOLVE("of_backend_update_data_vec", f_update_vec);
  RESOLVE("of_backend_update_data_mat", f_update_mat);
  RESOLVE("of_backend_warm_start", f_warm_start);
  RESOLVE("of_backend_cold_start", f_cold_start);
  RESOLVE("of_backend_update_settings", f_update_settings);
  return 1;
fail:
  return 0;
}

static int try_path(const char* path) {
  if (!path || !*path) return 0;
  g_lib = OF_DLOPEN(path);
  if (!g_lib) {
#ifdef _WIN32
    set_windows_error(path);
#else
    set_dl_error(path);
#endif
    return 0;
  }
  if (!resolve_all()) {
    OF_DLCLOSE(g_lib);
    g_lib = 0;
    return 0;
  }
  set_error("");
  return 1;
}

int of_api_load_backend(const char* explicit_path) {
  const char* env;
  const char* candidates[] = {
#ifdef _WIN32
    "osqp_fortran_bridge.dll",
    "backend\\bin\\osqp_fortran_bridge.dll",
    ".\\backend\\bin\\osqp_fortran_bridge.dll",
#else
    "libosqp_fortran_bridge.so",
    "backend/lib/libosqp_fortran_bridge.so",
    "./backend/lib/libosqp_fortran_bridge.so",
    "backend/bin/libosqp_fortran_bridge.so",
    "libosqp_fortran_bridge.dylib",
    "backend/lib/libosqp_fortran_bridge.dylib",
#endif
    NULL
  };
  int i;
  if (g_lib) return 1;
  if (try_path(explicit_path)) return 1;
  env = getenv("OSQP_FORTRAN_BRIDGE");
  if (try_path(env)) return 1;
  for (i = 0; candidates[i]; ++i) if (try_path(candidates[i])) return 1;
  return 0;
}

int of_api_backend_available(void) { return g_lib ? 1 : of_api_load_backend(NULL); }
void of_api_unload_backend(void) { if (g_lib) { OF_DLCLOSE(g_lib); g_lib = 0; } }
int of_api_last_error(char* buffer, int n) {
  if (buffer && n > 0) { snprintf(buffer, (size_t)n, "%s", g_error); buffer[n-1] = '\0'; }
  return (int)strlen(g_error);
}

#define NEED_BACKEND_RET(ret) do { if (!of_api_backend_available()) return (ret); } while (0)
void* of_api_create(int n,int m,int pn,const int* pp,const int* pi,const double* px,
                    int an,const int* ap,const int* ai,const double* ax,
                    const double* q,const double* l,const double* u,
                    const of_settings_c* s,int* ec) {
  NEED_BACKEND_RET(NULL);
  return f_create(n,m,pn,pp,pi,px,an,ap,ai,ax,q,l,u,s,ec);
}
void of_api_destroy(void* p) { if (g_lib && p) f_destroy(p); }
int of_api_solve(void* p) { NEED_BACKEND_RET(-100); return f_solve(p); }
int of_api_get_dimensions(void* p,int* n,int* m) { NEED_BACKEND_RET(-100); return f_get_dimensions(p,n,m); }
int of_api_get_solution(void* p,double* x,double* y,double* pc,double* dc,of_info_c* i) { NEED_BACKEND_RET(-100); return f_get_solution(p,x,y,pc,dc,i); }
int of_api_status_message(void* p,char* b,int n) { NEED_BACKEND_RET(-100); return f_status_message(p,b,n); }
int of_api_version(char* b,int n) { NEED_BACKEND_RET(-100); return f_version(b,n); }
int of_api_capabilities(void) { NEED_BACKEND_RET(0); return f_capabilities(); }
double of_api_infinity(void) { NEED_BACKEND_RET(1e30); return f_infinity(); }
int of_api_error_message(int c,char* b,int n) { NEED_BACKEND_RET(-100); return f_error_message(c,b,n); }
int of_api_get_settings(void* p,of_settings_c* s) { NEED_BACKEND_RET(-100); return f_get_settings(p,s); }
int of_api_update_data_vec(void* p,const double* q,int hq,const double* l,int hl,const double* u,int hu) { NEED_BACKEND_RET(-100); return f_update_vec(p,q,hq,l,hl,u,hu); }
int of_api_update_data_mat(void* p,const double* px,const int* pi,int pn,int hp,const double* ax,const int* ai,int an,int ha) { NEED_BACKEND_RET(-100); return f_update_mat(p,px,pi,pn,hp,ax,ai,an,ha); }
int of_api_warm_start(void* p,const double* x,int hx,const double* y,int hy) { NEED_BACKEND_RET(-100); return f_warm_start(p,x,hx,y,hy); }
int of_api_cold_start(void* p) { NEED_BACKEND_RET(-100); return f_cold_start(p); }
int of_api_update_settings(void* p,const of_settings_c* s) { NEED_BACKEND_RET(-100); return f_update_settings(p,s); }
