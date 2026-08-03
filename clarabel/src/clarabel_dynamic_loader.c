/* SPDX-License-Identifier: Apache-2.0
 * Runtime loader for the Clarabel Rust bridge.
 *
 * Keeping the bridge as a runtime-loaded shared library lets an ordinary
 * `fpm build` compile and link the Fortran package without requiring Cargo to
 * run inside FPM or requiring users to configure -L search paths manually.
 */
#include "../include/clarabel_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
typedef HMODULE bridge_handle_t;
#  define BRIDGE_HANDLE_NULL NULL
#else
#  include <dlfcn.h>
typedef void *bridge_handle_t;
#  define BRIDGE_HANDLE_NULL NULL
#endif

typedef void (*settings_default_fn)(clarabel_settings_c *);
typedef int32_t (*solver_create_fn)(const clarabel_csc_c *, const double *, size_t,
                                    const clarabel_csc_c *, const double *, size_t,
                                    const clarabel_cone_c *, size_t,
                                    const clarabel_settings_c *, void **,
                                    char *, size_t);
typedef int32_t (*solver_solve_fn)(void *, double *, size_t, double *, size_t,
                                   double *, size_t, clarabel_result_c *,
                                   char *, size_t);
typedef int32_t (*solver_update_fn)(void *, const double *, size_t,
                                    const double *, size_t, const double *,
                                    size_t, const double *, size_t,
                                    char *, size_t);
typedef int32_t (*solver_is_update_allowed_fn)(const void *);
typedef void (*solver_free_fn)(void *);

static bridge_handle_t bridge_handle = BRIDGE_HANDLE_NULL;
static settings_default_fn p_settings_default = NULL;
static solver_create_fn p_solver_create = NULL;
static solver_solve_fn p_solver_solve = NULL;
static solver_update_fn p_solver_update = NULL;
static solver_is_update_allowed_fn p_solver_is_update_allowed = NULL;
static solver_free_fn p_solver_free = NULL;
static int bridge_attempted = 0;
static char bridge_error[1024] = "Clarabel backend has not been loaded";

static void write_error(char *buffer, size_t capacity, const char *message) {
    size_t n;
    if (buffer == NULL || capacity == 0) return;
    if (message == NULL) message = "Clarabel backend error";
    n = strlen(message);
    if (n >= capacity) n = capacity - 1;
    memcpy(buffer, message, n);
    buffer[n] = '\0';
}

#ifdef _WIN32
#ifndef LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR
#  define LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR 0x00000100
#endif
#ifndef LOAD_LIBRARY_SEARCH_DEFAULT_DIRS
#  define LOAD_LIBRARY_SEARCH_DEFAULT_DIRS 0x00001000
#endif

static int library_file_exists(const char *path) {
    DWORD attrs;
    if (path == NULL || path[0] == '\0') return 0;
    attrs = GetFileAttributesA(path);
    return attrs != INVALID_FILE_ATTRIBUTES &&
           (attrs & FILE_ATTRIBUTE_DIRECTORY) == 0;
}

static bridge_handle_t open_library(const char *path) {
    bridge_handle_t handle;
    /* LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR makes dependencies beside the bridge
       discoverable. This matters for Rust's GNU Windows target, whose bridge
       can depend on MinGW runtime DLLs copied into rust_bridge\bin. */
    handle = LoadLibraryExA(path, NULL,
                            LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR |
                            LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
    if (handle == NULL && GetLastError() == ERROR_INVALID_PARAMETER) {
        /* Compatibility fallback for older Windows installations. */
        handle = LoadLibraryA(path);
    }
    return handle;
}

static void *load_symbol(bridge_handle_t handle, const char *name) {
    return (void *)(uintptr_t)GetProcAddress(handle, name);
}

static void close_library(bridge_handle_t handle) {
    if (handle != NULL) FreeLibrary(handle);
}

static void system_loader_error(char *buffer, size_t capacity) {
    DWORD code = GetLastError();
    char *text = NULL;
    DWORD n;
    if (buffer == NULL || capacity == 0) return;
    n = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER |
                       FORMAT_MESSAGE_FROM_SYSTEM |
                       FORMAT_MESSAGE_IGNORE_INSERTS,
                       NULL, code, 0, (LPSTR)&text, 0, NULL);
    if (n > 0 && text != NULL) {
        while (n > 0 && (text[n - 1] == '\r' || text[n - 1] == '\n')) {
            text[--n] = '\0';
        }
        snprintf(buffer, capacity, "%s", text);
        LocalFree(text);
    } else {
        snprintf(buffer, capacity, "Windows loader error %lu", (unsigned long)code);
    }
}
#else
static int library_file_exists(const char *path) {
    FILE *stream;
    if (path == NULL || path[0] == '\0') return 0;
    stream = fopen(path, "rb");
    if (stream == NULL) return 0;
    fclose(stream);
    return 1;
}

static bridge_handle_t open_library(const char *path) {
    return dlopen(path, RTLD_NOW | RTLD_LOCAL);
}

static void *load_symbol(bridge_handle_t handle, const char *name) {
    return dlsym(handle, name);
}

static void close_library(bridge_handle_t handle) {
    if (handle != NULL) dlclose(handle);
}

static void system_loader_error(char *buffer, size_t capacity) {
    const char *text = dlerror();
    snprintf(buffer, capacity, "%s", text != NULL ? text : "dynamic-loader error");
}
#endif

static int bind_symbols(bridge_handle_t handle) {
    void *symbol;

#define LOAD_BRIDGE_FUNCTION(target, name) \
    do { \
        symbol = load_symbol(handle, name); \
        _Static_assert(sizeof(target) == sizeof(symbol), \
                       "function and data pointers must have equal size"); \
        memcpy(&(target), &symbol, sizeof(target)); \
    } while (0)

    LOAD_BRIDGE_FUNCTION(p_settings_default, "clarabel_settings_default");
    LOAD_BRIDGE_FUNCTION(p_solver_create, "clarabel_solver_create");
    LOAD_BRIDGE_FUNCTION(p_solver_solve, "clarabel_solver_solve");
    LOAD_BRIDGE_FUNCTION(p_solver_update, "clarabel_solver_update");
    LOAD_BRIDGE_FUNCTION(p_solver_is_update_allowed,
                         "clarabel_solver_is_update_allowed");
    LOAD_BRIDGE_FUNCTION(p_solver_free, "clarabel_solver_free");

#undef LOAD_BRIDGE_FUNCTION

    return p_settings_default != NULL && p_solver_create != NULL &&
           p_solver_solve != NULL && p_solver_update != NULL &&
           p_solver_is_update_allowed != NULL && p_solver_free != NULL;
}

static int try_candidate(const char *path, char *last_error, size_t error_capacity) {
    bridge_handle_t handle;
    if (path == NULL || path[0] == '\0') return 0;
    handle = open_library(path);
    if (handle == BRIDGE_HANDLE_NULL) {
        system_loader_error(last_error, error_capacity);
        return 0;
    }
    if (!bind_symbols(handle)) {
        snprintf(last_error, error_capacity,
                 "library '%s' does not export the complete Clarabel C ABI", path);
        close_library(handle);
        p_settings_default = NULL;
        p_solver_create = NULL;
        p_solver_solve = NULL;
        p_solver_update = NULL;
        p_solver_is_update_allowed = NULL;
        p_solver_free = NULL;
        return 0;
    }
    bridge_handle = handle;
    return 1;
}

static int ensure_backend_loaded(void) {
    const char *explicit_path;
    char last_error[512] = "shared library not found";
    char existing_candidate[512] = "";
    int any_candidate_exists = 0;
    size_t i;
#ifdef _WIN32
    static const char *candidates[] = {
        "clarabel_fortran_bridge.dll",
        "rust_bridge\\bin\\clarabel_fortran_bridge.dll",
        "rust_bridge\\target\\release\\clarabel_fortran_bridge.dll"
    };
#elif defined(__APPLE__)
    static const char *candidates[] = {
        "libclarabel_fortran_bridge.dylib",
        "rust_bridge/bin/libclarabel_fortran_bridge.dylib",
        "rust_bridge/target/release/libclarabel_fortran_bridge.dylib"
    };
#else
    static const char *candidates[] = {
        "libclarabel_fortran_bridge.so",
        "rust_bridge/bin/libclarabel_fortran_bridge.so",
        "rust_bridge/target/release/libclarabel_fortran_bridge.so"
    };
#endif

    if (bridge_handle != BRIDGE_HANDLE_NULL) return 1;
    if (bridge_attempted) return 0;
    bridge_attempted = 1;

    explicit_path = getenv("CLARABEL_FORTRAN_BRIDGE");
    if (explicit_path != NULL && explicit_path[0] != '\0') {
        if (try_candidate(explicit_path, last_error, sizeof last_error)) return 1;
        snprintf(bridge_error, sizeof bridge_error,
                 "could not load CLARABEL_FORTRAN_BRIDGE='%s': %s",
                 explicit_path, last_error);
        return 0;
    }

    for (i = 0; i < sizeof candidates / sizeof candidates[0]; ++i) {
        if (library_file_exists(candidates[i])) {
            any_candidate_exists = 1;
            snprintf(existing_candidate, sizeof existing_candidate, "%s", candidates[i]);
        }
        if (try_candidate(candidates[i], last_error, sizeof last_error)) return 1;
    }

    if (any_candidate_exists) {
        snprintf(bridge_error, sizeof bridge_error,
                 "Clarabel backend library exists at '%.240s' but Windows/the dynamic loader "
                 "could not load it. A dependent runtime library may be missing. "
                 "Re-run scripts/build_backend.%s and use scripts/build_with_backend.%s run. "
                 "Last loader error: %.420s",
                 existing_candidate,
#ifdef _WIN32
                 "bat", "bat",
#else
                 "sh", "sh",
#endif
                 last_error);
    } else {
        snprintf(bridge_error, sizeof bridge_error,
                 "Clarabel backend library is not installed. A normal fpm build compiles "
                 "only the Fortran frontend. Run scripts/build_with_backend.%s run, or run "
                 "scripts/build_backend.%s once and then fpm run.",
#ifdef _WIN32
                 "bat", "bat"
#else
                 "sh", "sh"
#endif
                 );
    }
    return 0;
}

static void fallback_settings(clarabel_settings_c *s) {
    if (s == NULL) return;
    memset(s, 0, sizeof *s);
    s->max_iter = 200;
    s->time_limit = HUGE_VAL;
    s->verbose = 1;
    s->max_step_fraction = 0.99;
    s->tol_gap_abs = 1e-8;
    s->tol_gap_rel = 1e-8;
    s->tol_feas = 1e-8;
    s->tol_infeas_abs = 1e-8;
    s->tol_infeas_rel = 1e-8;
    s->tol_ktratio = 1e-6;
    s->reduced_tol_gap_abs = 5e-5;
    s->reduced_tol_gap_rel = 5e-5;
    s->reduced_tol_feas = 1e-4;
    s->reduced_tol_infeas_abs = 5e-5;
    s->reduced_tol_infeas_rel = 5e-5;
    s->reduced_tol_ktratio = 1e-4;
    s->equilibrate_enable = 1;
    s->equilibrate_max_iter = 10;
    s->equilibrate_min_scaling = 1e-4;
    s->equilibrate_max_scaling = 1e4;
    s->linesearch_backtrack_step = 0.8;
    s->min_switch_step_length = 0.1;
    s->min_terminate_step_length = 1e-4;
    s->direct_kkt_solver = 1;
    s->direct_solve_method = 1;
    s->static_regularization_enable = 1;
    s->static_regularization_constant = 1e-8;
    s->static_regularization_proportional = 4.930380657631324e-32;
    s->dynamic_regularization_enable = 1;
    s->dynamic_regularization_eps = 1e-13;
    s->dynamic_regularization_delta = 2e-7;
    s->iterative_refinement_enable = 1;
    s->iterative_refinement_reltol = 1e-13;
    s->iterative_refinement_abstol = 1e-12;
    s->iterative_refinement_max_iter = 10;
    s->iterative_refinement_stop_ratio = 5.0;
    s->presolve_enable = 1;
}

void clarabel_settings_default(clarabel_settings_c *settings) {
    if (ensure_backend_loaded()) {
        p_settings_default(settings);
    } else {
        fallback_settings(settings);
    }
}

int32_t clarabel_solver_create(const clarabel_csc_c *P, const double *q, size_t q_len,
                               const clarabel_csc_c *A, const double *b, size_t b_len,
                               const clarabel_cone_c *cones, size_t ncones,
                               const clarabel_settings_c *settings, void **out,
                               char *error, size_t error_capacity) {
    if (out != NULL) *out = NULL;
    if (!ensure_backend_loaded()) {
        write_error(error, error_capacity, bridge_error);
        return -100;
    }
    return p_solver_create(P, q, q_len, A, b, b_len, cones, ncones,
                           settings, out, error, error_capacity);
}

int32_t clarabel_solver_solve(void *solver, double *x, size_t x_len,
                              double *z, size_t z_len, double *s, size_t s_len,
                              clarabel_result_c *result, char *error,
                              size_t error_capacity) {
    if (!ensure_backend_loaded()) {
        write_error(error, error_capacity, bridge_error);
        return -100;
    }
    return p_solver_solve(solver, x, x_len, z, z_len, s, s_len,
                          result, error, error_capacity);
}

int32_t clarabel_solver_update(void *solver,
                               const double *p_values, size_t p_len,
                               const double *a_values, size_t a_len,
                               const double *q, size_t q_len,
                               const double *b, size_t b_len,
                               char *error, size_t error_capacity) {
    if (!ensure_backend_loaded()) {
        write_error(error, error_capacity, bridge_error);
        return -100;
    }
    return p_solver_update(solver, p_values, p_len, a_values, a_len,
                           q, q_len, b, b_len, error, error_capacity);
}

int32_t clarabel_solver_is_update_allowed(const void *solver) {
    if (!ensure_backend_loaded()) return 0;
    return p_solver_is_update_allowed(solver);
}

void clarabel_solver_free(void *solver) {
    if (solver == NULL) return;
    if (ensure_backend_loaded()) p_solver_free(solver);
}
