! SPDX-License-Identifier: GPL-2.0-or-later
! Public modern Fortran API for the computational core of R package svd 0.5.8.
module svd
   use svd_kinds, only : dp
   use rspectra, only : linear_operator
   use svd_extmat, only : ematmul, extmat_callbacks, extmat_ncol, extmat_nrow
   use svd_extmat, only : extmat_operator, extmat_vector, is_extmat, make_extmat, materialize_extmat
   use svd_extmat, only : svd_invalid_callback, svd_invalid_shape, svd_success
   use svd_solvers, only : propack_svd, propack_svd_dense, propack_svd_operator
   use svd_solvers, only : propack_svd_result, solver_invalid_input, solver_invalid_warm_start, solver_success
   use svd_solvers, only : svd_options, trlan_eigen, trlan_eigen_dense, trlan_eigen_operator, trlan_eigen_result
   use svd_solvers, only : trlan_svd, trlan_svd_dense, trlan_svd_operator, trlan_svd_result
   use svd_solvers, only : ztrlan_svd_complex, ztrlan_svd_result
   implicit none
   private

   public :: dp
   public :: linear_operator
   public :: extmat_callbacks, extmat_operator
   public :: svd_options
   public :: propack_svd_result, trlan_svd_result, ztrlan_svd_result, trlan_eigen_result
   public :: svd_success, svd_invalid_shape, svd_invalid_callback
   public :: solver_success, solver_invalid_input, solver_invalid_warm_start
   public :: make_extmat, is_extmat, extmat_nrow, extmat_ncol, ematmul, materialize_extmat, extmat_vector
   public :: propack_svd, propack_svd_dense, propack_svd_operator
   public :: trlan_svd, trlan_svd_dense, trlan_svd_operator
   public :: trlan_eigen, trlan_eigen_dense, trlan_eigen_operator
   public :: ztrlan_svd_complex

end module svd
