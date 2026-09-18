! SPDX-License-Identifier: GPL-3.0-only
! Derived from dccmidas 0.1.3 by Vincenzo Candila; see PROVENANCE.md.
module dccmidas_types
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: DCCMIDAS_SUCCESS = 0
   integer, parameter, public :: DCCMIDAS_INVALID_INPUT = 1
   integer, parameter, public :: DCCMIDAS_LINALG_ERROR = 2
   integer, parameter, public :: DCCMIDAS_OPTIMIZATION_ERROR = 3
   integer, parameter, public :: DCCMIDAS_UNSUPPORTED = 4

   type, public :: dcc_matrices
      real(dp), allocatable :: h(:, :, :)
      real(dp), allocatable :: r(:, :, :)
      real(dp), allocatable :: r_bar(:, :, :)
   end type dcc_matrices

   type, public :: dcc_fit_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: conditional_sd(:, :)
      real(dp), allocatable :: standardized_residuals(:, :)
      type(dcc_matrices) :: matrices
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: status = DCCMIDAS_INVALID_INPUT
      logical :: converged = .false.
      character(len=96) :: message = 'not run'
   end type dcc_fit_result

   type, public :: bekk_fit_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: h(:, :, :)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: status = DCCMIDAS_INVALID_INPUT
      logical :: converged = .false.
      character(len=96) :: message = 'not run'
   end type bekk_fit_result

end module dccmidas_types
