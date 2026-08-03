#ifndef PRINTING_H_
#define PRINTING_H_
#include "osqp_configure.h"
#ifdef __cplusplus
extern "C" {
#endif
#ifdef OSQP_USE_LONG
#define OSQP_INT_FMT "lld"
#else
#define OSQP_INT_FMT "d"
#endif
#if __STDC_VERSION__ >= 199901L
#define c_eprint(...) c_print("ERROR in %s: ", __func__); c_print(__VA_ARGS__); c_print("\n")
#else
#define c_eprint(...) c_print("ERROR in %s: ", __FUNCTION__); c_print(__VA_ARGS__); c_print("\n")
#endif
#ifdef OSQP_ENABLE_PRINTING
#include <stdio.h>
#define c_print printf
#else
#undef c_eprint
#define c_print(...)
#define c_eprint(...)
#endif
#ifdef __cplusplus
}
#endif
#endif
