#ifndef __CONFIG_H__
#define __CONFIG_H__

#define SCIP_BUILD_TYPE "Release"
#define SCIP_VERSION_MAJOR 10
#define SCIP_VERSION_MINOR 0
#define SCIP_VERSION_PATCH 1
#define SCIP_VERSION_API 156
/* #undef BMS_NOBLOCKMEM */
/* #undef SCIP_NOBUFFERMEM */
/* #undef WITH_DEBUG_SOLUTION */
#ifdef _WIN32
#define SCIP_NO_SIGACTION
#define SCIP_NO_STRTOK_R
#else
/* #undef SCIP_NO_SIGACTION */
/* #undef SCIP_NO_STRTOK_R */
#endif
/* #undef TPI_NONE */
#define TPI_TNY
/* #undef TPI_OMP */
#define SCIP_THREADSAFE
#define WITH_SCIPDEF
/* #undef SCIP_WITH_LAPACK */
/* #undef SCIP_WITH_PAPILO */
/* #undef SCIP_WITH_ZLIB */
/* #undef SCIP_WITH_READLINE */
/* #undef SCIP_WITH_GMP */
/* #undef SCIP_WITH_MPFR */
/* #undef SCIP_WITH_QSOPTEX */
/* #undef SCIP_WITH_LPSCHECK */
/* #undef SCIP_WITH_ZIMPL */
/* #undef SCIP_WITH_AMPL */
#define SCIP_ROUNDING_FE
/* #undef SCIP_ROUNDING_FP */
/* #undef SCIP_ROUNDING_MS */
/* #undef SCIP_WITH_EXACTSOLVE */
/* #undef SCIP_WITH_EGLIB */
/* #undef SCIP_WITH_BOOST */
/* #undef SCIP_CHECK_STAGE */

#endif
