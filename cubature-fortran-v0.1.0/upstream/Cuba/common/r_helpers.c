#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <R.h>

void invoke_r_exit(void) {
    Rf_error("Error invoked from Cuba library");
}

/**
 * safe_sprintf:
 *   Append formatted data to a buffer, tracking the write pointer and
 *   available space. On error or truncation, calls Rf_error().
 *
 * @param poe    Pointer to current write pointer into your output buffer.
 * @param pavail Pointer to remaining space in your output buffer.
 * @param fmt    printf-style format string.
 * @param ...    Arguments matching fmt specifiers.
 */
void safe_sprintf(char **poe, size_t *pavail, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int n = vsnprintf(*poe, *pavail, fmt, ap);
  va_end(ap);
  
  /* vsnprintf returns <0 on encoding error, or the number it WOULD
     have written if avail had been large enough. */
  if (n < 0 || (size_t)n >= *pavail) {
    Rf_error("Error while doing safe_sprintf");
  }
  
  /* advance write pointer and shrink remaining space */
  *poe    += n;
  *pavail -= (size_t)n;
}

