! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_types
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok
   implicit none
   private
   public :: test_result, location_scatter_result, spatial_sign_result

   type :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      integer :: replications = 0
      integer :: status = icsnp_ok
      character(len=96) :: method = ""
   end type test_result

   type :: location_scatter_result
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: scatter(:,:)
      integer :: iterations = 0
      integer :: status = icsnp_ok
   end type location_scatter_result

   type :: spatial_sign_result
      real(dp), allocatable :: signs(:,:)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: shape(:,:)
      integer :: status = icsnp_ok
   end type spatial_sign_result
end module icsnp_types
