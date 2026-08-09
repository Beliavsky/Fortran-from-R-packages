! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_types
   use coneproj_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: coneproj_success = 0
   integer, parameter, public :: coneproj_max_iter = 1
   integer, parameter, public :: coneproj_invalid_input = 2
   integer, parameter, public :: coneproj_singular = 3


   type, public :: qr_result
      real(dp), allocatable :: q(:,:)
      integer :: rank = 0
   end type qr_result

   type, public :: cone_result
      real(dp), allocatable :: fit(:)
      real(dp), allocatable :: coefs(:)
      real(dp), allocatable :: xmat(:,:)
      integer, allocatable :: face(:)
      integer :: df = 0
      integer :: steps = 0
      integer :: status = coneproj_success
   end type cone_result

   type, public :: qprog_result
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: xmat(:,:)
      integer, allocatable :: face(:)
      real(dp) :: objective = 0.0_dp
      integer :: df = 0
      integer :: steps = 0
      integer :: status = coneproj_success
   end type qprog_result

   type, public :: constreg_result
      real(dp), allocatable :: coefs(:)
      real(dp), allocatable :: constrained_fit(:)
      real(dp), allocatable :: unconstrained_fit(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: tstat(:)
      real(dp), allocatable :: pvalues(:)
      real(dp), allocatable :: lower_coef(:)
      real(dp), allocatable :: upper_coef(:)
      real(dp), allocatable :: lower_fit(:)
      real(dp), allocatable :: upper_fit(:)
      integer, allocatable :: face(:)
      real(dp) :: pvalue_test = -1.0_dp
      real(dp) :: sse0 = 0.0_dp
      real(dp) :: sse1 = 0.0_dp
      real(dp) :: bstat = 0.0_dp
      integer :: df = 0
      integer :: status = coneproj_success
   end type constreg_result

   type, public :: shapereg_result
      real(dp), allocatable :: coefs(:)
      real(dp), allocatable :: constrained_fit(:)
      real(dp), allocatable :: linear_fit(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: pvalues(:)
      real(dp) :: pvalue_test = -1.0_dp
      real(dp) :: sse0 = 0.0_dp
      real(dp) :: sse1 = 0.0_dp
      integer :: shape = 0
      integer :: df = 0
      integer :: status = coneproj_success
   end type shapereg_result

end module coneproj_types
