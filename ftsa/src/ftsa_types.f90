! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_types
   use r_kinds, only : dp
   implicit none
   private

   type, public :: ftsm_fit
      real(dp), allocatable :: mean_curve(:)
      real(dp), allocatable :: basis(:, :)
      real(dp), allocatable :: scores(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: varprop(:)
      real(dp), allocatable :: weights(:)
      integer :: info = 0
   end type ftsm_fit

   type, public :: hurst_result
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: eigenfunctions(:, :)
      real(dp), allocatable :: prop_eigenvalue(:)
      real(dp) :: h_est = 0.0_dp
      real(dp) :: d_est = 0.0_dp
      real(dp) :: alpha_est = 0.0_dp
      integer :: ncomp = 0
      integer :: info = 0
   end type hurst_result

end module ftsa_types
