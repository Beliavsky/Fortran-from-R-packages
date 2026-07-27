! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_types
  use pbo_kinds, only : dp
  implicit none
  private

  type, public :: pbo_result
    integer :: n_observations = 0
    integer :: n_strategies = 0
    integer :: n_subsets = 0
    integer :: n_cases = 0
    real(dp) :: threshold = 0.0_dp
    real(dp) :: inf_sub = 6.0_dp
    real(dp) :: phi = 0.0_dp
    real(dp) :: slope = 0.0_dp
    real(dp) :: intercept = 0.0_dp
    real(dp) :: adjusted_r2 = 0.0_dp
    real(dp) :: below_threshold = 0.0_dp
    real(dp) :: degradation_slope = 0.0_dp
    real(dp) :: degradation_intercept = 0.0_dp
    real(dp) :: degradation_r2 = 0.0_dp
    logical :: success = .false.
    character(len=:), allocatable :: message
    integer, allocatable :: combos(:,:)
    real(dp), allocatable :: performance_is(:,:)
    real(dp), allocatable :: performance_oos(:,:)
    integer, allocatable :: selected_is(:)
    integer, allocatable :: selected_oos(:)
    real(dp), allocatable :: oos_rank(:)
    real(dp), allocatable :: omega_bar(:)
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: selected_pairs(:,:)
  end type pbo_result

  type, public :: dominance_result
    logical :: success = .false.
    character(len=:), allocatable :: message
    real(dp), allocatable :: performance(:)
    real(dp), allocatable :: cdf_selected(:)
    real(dp), allocatable :: cdf_all(:)
    real(dp), allocatable :: sd2_difference(:)
    real(dp), allocatable :: integrated_difference(:)
  end type dominance_result

  type, public :: selection_result
    integer, allocatable :: strategy(:)
    integer, allocatable :: frequency(:)
  end type selection_result
end module pbo_types
