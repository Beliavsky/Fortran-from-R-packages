#ifndef HIGHS_FORTRAN_RCPP_COMPAT_H
#define HIGHS_FORTRAN_RCPP_COMPAT_H
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <stdexcept>
namespace Rcpp {
inline std::ostream& Rcout = std::cout;
[[noreturn]] inline void stop(const char* message) { throw std::runtime_error(message); }
}
inline int Rprintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  const int result = std::vprintf(format, args);
  va_end(args);
  return result;
}
#endif
