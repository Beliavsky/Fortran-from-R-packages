! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_types
  use pa_kinds, only: dp
  implicit none
  private

  type, public :: exposure_result
    integer, allocatable :: group(:)
    real(dp), allocatable :: portfolio(:)
    real(dp), allocatable :: benchmark(:)
    real(dp), allocatable :: difference(:)
  end type exposure_result

  type, public :: exposure_multi_result
    integer, allocatable :: group(:)
    integer, allocatable :: period(:)
    real(dp), allocatable :: portfolio(:, :)
    real(dp), allocatable :: benchmark(:, :)
    real(dp), allocatable :: difference(:, :)
  end type exposure_multi_result

  type, public :: brinson_period_result
    integer, allocatable :: category(:)
    real(dp), allocatable :: weight_portfolio(:)
    real(dp), allocatable :: weight_benchmark(:)
    real(dp), allocatable :: return_portfolio(:)
    real(dp), allocatable :: return_benchmark(:)
    real(dp), allocatable :: category_effect(:, :) ! allocation, selection, interaction
    real(dp) :: q(4) = 0.0_dp                 ! q4, q3, q2, q1
    real(dp) :: aggregate(4) = 0.0_dp         ! allocation, selection, interaction, active
    integer :: status = 0
  end type brinson_period_result

  type, public :: brinson_multi_result
    integer, allocatable :: period(:)
    integer, allocatable :: category(:)
    real(dp), allocatable :: weight_portfolio(:, :)
    real(dp), allocatable :: weight_benchmark(:, :)
    real(dp), allocatable :: return_portfolio(:, :)
    real(dp), allocatable :: return_benchmark(:, :)
    real(dp), allocatable :: category_effect(:, :, :)
    real(dp), allocatable :: q(:, :)            ! q4, q3, q2, q1 by period
    real(dp), allocatable :: raw(:, :)          ! allocation, selection, interaction, active
    integer :: status = 0
  end type brinson_multi_result

  type, public :: attribution_summary
    real(dp), allocatable :: raw(:, :)
    real(dp), allocatable :: aggregate(:)
    real(dp), allocatable :: linking_coefficient(:)
    integer :: status = 0
  end type attribution_summary

  type, public :: regression_period_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: active_exposure(:)
    real(dp), allocatable :: contribution(:)
    real(dp) :: portfolio_return = 0.0_dp
    real(dp) :: benchmark_return = 0.0_dp
    real(dp) :: active_return = 0.0_dp
    integer :: rank = 0
    integer :: status = 0
  end type regression_period_result

  type, public :: regression_multi_result
    integer, allocatable :: period(:)
    real(dp), allocatable :: coefficients(:, :)
    real(dp), allocatable :: active_exposure(:, :)
    real(dp), allocatable :: contribution(:, :)
    real(dp), allocatable :: portfolio_return(:)
    real(dp), allocatable :: benchmark_return(:)
    real(dp), allocatable :: active_return(:)
    integer, allocatable :: rank(:)
    integer :: status = 0
  end type regression_multi_result

  type, public :: regression_summary
    real(dp), allocatable :: raw(:, :)
    real(dp), allocatable :: aggregate(:)
    real(dp), allocatable :: linking_coefficient(:)
    integer :: status = 0
  end type regression_summary

end module pa_types
