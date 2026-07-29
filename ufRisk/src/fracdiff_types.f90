! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_types
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok
   implicit none
   private

   type, public :: fracdiff_model
      integer :: n = 0
      integer :: p = 0
      integer :: q = 0
      integer :: m_terms = 100
      integer :: status = fd_ok
      integer :: outer_iterations = 0
      integer :: function_evaluations = 0
      integer :: gradient_evaluations = 0
      real(dp) :: d = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: fnorm_min = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: d_tol = 0.0_dp
      real(dp) :: h = 0.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: correlation(:,:)
      real(dp), allocatable :: std_error(:)
      real(dp), allocatable :: hessian(:,:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fitted(:)
      character(len=:), allocatable :: message
   end type fracdiff_model

   type, public :: fracdiff_simulation
      integer :: n_start = 0
      integer :: status = fd_ok
      real(dp) :: d = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: series(:)
      character(len=:), allocatable :: message
   end type fracdiff_simulation

   type, public :: fractional_d_estimate
      integer :: status = fd_ok
      real(dp) :: d = 0.0_dp
      real(dp) :: sd_asymptotic = 0.0_dp
      real(dp) :: sd_regression = 0.0_dp
      character(len=:), allocatable :: message
   end type fractional_d_estimate

   type, public :: fracdiff_summary
      integer :: degrees_freedom = 0
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp), allocatable :: coefficients(:,:)
   end type fracdiff_summary

end module fracdiff_types
