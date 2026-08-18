#ifndef _RCUBA_HELPER
#define _RCUBA_HELPER

#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
void invoke_r_exit(void);
void R_print(char *x);
void safe_sprintf(char **poe, size_t *pavail, const char *fmt, ...);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
