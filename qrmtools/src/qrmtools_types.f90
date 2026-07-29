! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_types
  use qrmtools_kinds, only : dp
  implicit none
  private
  public :: fit_result, hill_result, brownian_result, greeks_result
  public :: rearrangement_result, ra_bounds_result, allocation_result, test_result, garch_result

  type :: fit_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: standard_errors(:)
    real(dp) :: log_likelihood = 0.0_dp
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
  end type fit_result

  type :: hill_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    integer, allocatable :: k(:)
    real(dp), allocatable :: probability(:)
    real(dp), allocatable :: tail_index(:)
    real(dp), allocatable :: ci_low(:)
    real(dp), allocatable :: ci_high(:)
  end type hill_result

  type :: brownian_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp), allocatable :: paths(:,:,:)
  end type brownian_result

  type :: greeks_result
    real(dp) :: price = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp) :: vega = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: vanna = 0.0_dp
    real(dp) :: vomma = 0.0_dp
    logical :: ok = .false.
    character(len=256) :: message = ''
  end type greeks_result

  type :: rearrangement_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp) :: bound = 0.0_dp
    real(dp) :: tolerance = 0.0_dp
    logical :: converged = .false.
    integer :: iterations = 0
    real(dp), allocatable :: optimal_values(:)
    real(dp), allocatable :: rearranged(:,:)
    real(dp), allocatable :: optimal_row(:)
  end type rearrangement_result

  type :: ra_bounds_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp) :: bounds(2) = 0.0_dp
    real(dp) :: relative_gap = 0.0_dp
    real(dp) :: tolerances(3) = 0.0_dp
    logical :: converged(3) = .false.
    integer :: n_used = 0
    type(rearrangement_result) :: lower
    type(rearrangement_result) :: upper
  end type ra_bounds_result

  type :: allocation_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp), allocatable :: allocation(:)
    real(dp), allocatable :: standard_error(:)
    real(dp), allocatable :: conditional(:,:)
    integer :: n = 0
  end type allocation_result

  type :: test_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp), allocatable :: distances(:)
  end type test_result

  type :: garch_result
    logical :: ok = .false.
    character(len=256) :: message = ''
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: log_likelihood = 0.0_dp
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
  end type garch_result
end module qrmtools_types
