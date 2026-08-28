! SPDX-License-Identifier: GPL-2.0-or-later
module nleqslv_blas_lapack
   use la_blas_aux, only : idamax => la_idamax
   use la_blas_d, only : daxpy => la_daxpy
   use la_blas_d, only : dcopy => la_dcopy
   use la_blas_d, only : ddot => la_ddot
   use la_blas_d, only : dgemv => la_dgemv
   use la_blas_d, only : dnrm2 => la_dnrm2
   use la_blas_d, only : drot => la_drot
   use la_blas_d, only : dscal => la_dscal
   use la_blas_d, only : dtrmm => la_dtrmm
   use la_blas_d, only : dtrmv => la_dtrmv
   use la_blas_d, only : dtrsv => la_dtrsv
   use la_lapack_d, only : dgeqrf => la_dgeqrf
   use la_lapack_d, only : dlacpy => la_dlacpy
   use la_lapack_d, only : dlamch => la_dlamch
   use la_lapack_d, only : dlantr => la_dlantr
   use la_lapack_d, only : dlartg => la_dlartg
   use la_lapack_d, only : dorgqr => la_dorgqr
   use la_lapack_d, only : dormqr => la_dormqr
   use la_lapack_d, only : dtrcon => la_dtrcon
   implicit none
   private
   public :: daxpy, dcopy, ddot, dgemv, dgeqrf, dlacpy, dlamch, dlantr
   public :: dlartg, dnrm2, dorgqr, dormqr, drot, dscal, dtrcon, dtrmm
   public :: dtrmv, dtrsv, idamax
end module nleqslv_blas_lapack
