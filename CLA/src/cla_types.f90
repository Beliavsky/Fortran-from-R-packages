! SPDX-License-Identifier: GPL-3.0-or-later
module cla_types
   use kind_mod, only: dp
   implicit none
   private

   integer, parameter, public :: cla_success = 0
   integer, parameter, public :: cla_invalid_input = 1
   integer, parameter, public :: cla_infeasible_bounds = 2
   integer, parameter, public :: cla_singular_system = 3
   integer, parameter, public :: cla_no_improvement = 4
   integer, parameter, public :: cla_out_of_range = 5
   integer, parameter, public :: cla_garch_failure = 6

   type, public :: cla_result_t
      real(dp), allocatable :: weights(:,:)
      logical, allocatable :: free_mask(:,:)
      real(dp), allocatable :: lambdas(:)
      real(dp), allocatable :: gammas(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: mu(:)
      integer :: n_assets = 0
      integer :: n_turning = 0
      integer :: info = cla_success
      integer :: warnings = 0
   end type cla_result_t

   type, public :: cla_path_query_t
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: weights(:,:)
      integer :: info = cla_success
   end type cla_path_query_t

   type, public :: cla_garch_result_t
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: forecast_sigma(:)
      real(dp), allocatable :: conditional_sigma(:,:)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: shape(:)
      real(dp), allocatable :: alpha(:,:)
      real(dp), allocatable :: beta(:,:)
      logical, allocatable :: converged(:)
      integer :: info = cla_success
   end type cla_garch_result_t

end module cla_types
