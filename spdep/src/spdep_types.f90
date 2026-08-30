! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_types
   use spdep_kinds, only : dp
   implicit none
   private

   type, public :: int_vector
      integer, allocatable :: values(:)
   end type int_vector

   type, public :: real_vector
      real(dp), allocatable :: values(:)
   end type real_vector

   type, public :: neighbor_list
      type(int_vector), allocatable :: neighbors(:)
      logical :: self_included = .false.
   contains
      procedure :: size => nb_size
   end type neighbor_list

   type, public :: spatial_weights
      type(neighbor_list) :: nb
      type(real_vector), allocatable :: weights(:)
      character(len=8) :: style = "B"
      logical :: zero_policy = .true.
   contains
      procedure :: size => listw_size
   end type spatial_weights

   type, public :: knn_result
      integer, allocatable :: index(:, :)
      real(dp), allocatable :: distance(:, :)
   end type knn_result

   type, public :: weights_constants
      integer :: n = 0
      real(dp) :: s0 = 0.0_dp
      real(dp) :: s1 = 0.0_dp
      real(dp) :: s2 = 0.0_dp
   end type weights_constants

   type, public :: spatial_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: expectation = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: z_score = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: kurtosis = 0.0_dp
   end type spatial_test_result

   type, public :: local_stat_result
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: expectation(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: z_score(:)
      real(dp), allocatable :: p_value(:)
   end type local_stat_result

   type, public :: eb_result
      real(dp), allocatable :: raw(:)
      real(dp), allocatable :: estimate(:)
      real(dp) :: global_rate = 0.0_dp
      real(dp) :: dispersion = 0.0_dp
   end type eb_result

   type, public :: mst_result
      integer, allocatable :: from(:)
      integer, allocatable :: to(:)
      real(dp), allocatable :: cost(:)
      real(dp) :: total_cost = 0.0_dp
   end type mst_result

   type, public :: spatial_delta_result
      real(dp) :: delta = 0.0_dp
      real(dp) :: rv = 0.0_dp
      real(dp) :: expectation = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: variance_product = 0.0_dp
      real(dp) :: z_score = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: skewness = 0.0_dp
      real(dp) :: excess_kurtosis = 0.0_dp
   end type spatial_delta_result

contains

   pure integer function nb_size(self) result(n)
      class(neighbor_list), intent(in) :: self !! Neighbor-list object whose number of regions is requested.

      if (allocated(self%neighbors)) then
         n = size(self%neighbors)
      else
         n = 0
      end if
   end function nb_size

   pure integer function listw_size(self) result(n)
      class(spatial_weights), intent(in) :: self !! Spatial-weights object whose number of regions is requested.

      n = self%nb%size()
   end function listw_size

end module spdep_types
