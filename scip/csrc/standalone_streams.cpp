/* Apache-2.0 compatibility layer for the R-patched vendored sources. */
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <iostream>

extern "C" void Rprintf(const char *fmt, ...)
{
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
}
extern "C" void REprintf(const char *fmt, ...)
{
    va_list ap; va_start(ap, fmt); vfprintf(stderr, fmt, ap); va_end(ap);
}
extern "C" void Rvprintf(const char *fmt, va_list ap) { vprintf(fmt, ap); }
extern "C" void REvprintf(const char *fmt, va_list ap) { vfprintf(stderr, fmt, ap); }
extern "C" void Rf_error(const char *fmt, ...)
{
    va_list ap; va_start(ap, fmt); vfprintf(stderr, fmt, ap); va_end(ap);
    fputc('\n', stderr); abort();
}
std::ostream& r_cout() { return std::cout; }
std::ostream& r_cerr() { return std::cerr; }
