! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Public data structures for the modern Fortran translation of R package pan.
module pan_types
   use pan_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: PAN_OK = 0
   integer, parameter, public :: PAN_ERR_DIMENSION = 1
   integer, parameter, public :: PAN_ERR_ARGUMENT = 2
   integer, parameter, public :: PAN_ERR_LINALG = 3
   integer, parameter, public :: PAN_ERR_NUMERIC = 4

   type, public :: pan_prior
      real(dp) :: a = 0.0_dp
      real(dp), allocatable :: binv(:, :)
      real(dp) :: c = 0.0_dp
      real(dp), allocatable :: dinv(:, :)
   end type pan_prior

   type, public :: pan_bd_prior
      real(dp) :: a = 0.0_dp
      real(dp), allocatable :: binv(:, :)
      real(dp), allocatable :: c(:)
      real(dp), allocatable :: dinv(:, :, :)
   end type pan_bd_prior

   type, public :: pan_state
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: psi(:, :)
      real(dp), allocatable :: y(:, :)
   end type pan_state

   type, public :: pan_bd_state
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: psi(:, :, :)
      real(dp), allocatable :: y(:, :)
   end type pan_bd_state

   type, public :: pan_result
      real(dp), allocatable :: beta(:, :, :)
      real(dp), allocatable :: sigma(:, :, :)
      real(dp), allocatable :: psi(:, :, :)
      real(dp), allocatable :: y(:, :)
      type(pan_state) :: last
      integer :: status = PAN_OK
      character(len=:), allocatable :: message
   end type pan_result

   type, public :: pan_bd_result
      real(dp), allocatable :: beta(:, :, :)
      real(dp), allocatable :: sigma(:, :, :)
      real(dp), allocatable :: psi(:, :, :, :)
      real(dp), allocatable :: y(:, :)
      type(pan_bd_state) :: last
      integer :: status = PAN_OK
      character(len=:), allocatable :: message
   end type pan_bd_result

   type, public :: ecme_result
      real(dp), allocatable :: beta(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp), allocatable :: psi(:, :)
      real(dp), allocatable :: cov_beta(:, :)
      real(dp), allocatable :: bhat(:, :)
      real(dp), allocatable :: cov_b(:, :, :)
      real(dp), allocatable :: loglik(:)
      integer :: iter = 0
      logical :: converged = .false.
      integer :: status = PAN_OK
      character(len=:), allocatable :: message
   end type ecme_result

contains

end module pan_types
