#ifndef R_STREAMS_H
#define R_STREAMS_H

/*
 * R-compatibility declarations for vendored C/C++ sources.
 *
 * Printing: Rprintf, REprintf, Rvprintf, REvprintf
 *   from R_ext/Print.h — stable R API.
 *
 * Errors: Rf_error
 *   from R_ext/Error.h — stable R API. Used in place of abort().
 *
 * C++ streams: r_cout(), r_cerr()
 *   Defined in r_streams.cpp. Used in place of std::cout/std::cerr.
 */

#include <stdarg.h>

#ifdef __cplusplus
#include <ostream>
std::ostream& r_cout();
std::ostream& r_cerr();
extern "C" {
#endif

extern void Rprintf(const char *, ...);
extern void REprintf(const char *, ...);
extern void Rvprintf(const char *, va_list);
extern void REvprintf(const char *, va_list);
extern void Rf_error(const char *, ...) __attribute__((noreturn));

#ifdef __cplusplus
}
#endif

#endif
