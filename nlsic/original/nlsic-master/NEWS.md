## nlsic 1.2.0 (2026-01-27)

 * nlsic(): added reuse of jacobian in non linear iterations
   - con$reuse_jac=TRUE (do it by default)
   - con$max_reuse=5L (default maximal number of jacobian reuses)
 * nlsic(): con$btfrac default value is brought to 0.8 instead of 0.5.
 * lsi() and lsi_ln() return x with `aqr` attribute containing QR decomposition of
   matrix a. This is used to spare some calculations in non linear iterations

## nlsic 1.1.1 (2025-06-26)

 * ls_ln() and ls_0(): fixed case when rhs, b, is a multi-column matrix

## nlsic 1.1.0 (2025-05-14)

 * lsie_ln(): added mnorm and x0 parameters
 * lsi_ln(): rewritten from scratch to use lsie_ln() recursively
 * added pnull() and ls_0() functions
 * ls_ln() now gets mnorm and x0 parameters
 * complited unit tests
 * warning() is called in problematic situations. So user can intercept
   and debug them

## nlsic 1.0.4 (2023-06-26)

 * nlsic(): fixed error reporting when step has NA
 * nlsic(): fixed error reporting when inequalities are not satisfied

## nlsic 1.0.3 (2022-04-29)

 * nlsic(): fixed error when par is a matrix and u=NULL
 * nlsic(): added `paro` field in returned result with original structure of `par`

## nlsic 1.0.2 (2022-04-12)

 * now, lsi() returns NA vector if rhs is all NA.

## nlsic 1.0.1 (2020-01-10)

 * Added a `NEWS.md` file to track changes to the package.
 * minor fixes in documentation
