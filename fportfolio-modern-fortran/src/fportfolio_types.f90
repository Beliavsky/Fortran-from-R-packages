! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_types
  use fportfolio_kinds, only: dp
  implicit none
  private

  type, public :: linear_constraints
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp) :: budget = 1.0_dp
    logical :: budget_equality = .true.
    real(dp) :: budget_min = 1.0_dp
    real(dp) :: budget_max = 1.0_dp
    logical :: has_target_return = .false.
    real(dp) :: target_return = 0.0_dp
    real(dp), allocatable :: a_eq(:,:)
    real(dp), allocatable :: b_eq(:)
    real(dp), allocatable :: a_ineq(:,:)
    real(dp), allocatable :: b_ineq(:)
  end type linear_constraints

  type, public :: optimizer_result
    real(dp), allocatable :: weights(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: risk = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = 1
    logical :: converged = .false.
    character(len=160) :: message = "not run"
  end type optimizer_result

  type, public :: frontier_result
    real(dp), allocatable :: target_return(:)
    real(dp), allocatable :: risk(:)
    real(dp), allocatable :: weights(:,:)
    logical, allocatable :: feasible(:)
  end type frontier_result

  type, public :: risk_report
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: volatility = 0.0_dp
    real(dp) :: var = 0.0_dp
    real(dp) :: es = 0.0_dp
    real(dp) :: max_drawdown = 0.0_dp
    real(dp) :: dar = 0.0_dp
    real(dp) :: cdar = 0.0_dp
    real(dp), allocatable :: marginal(:)
    real(dp), allocatable :: component(:)
    real(dp), allocatable :: budget(:)
  end type risk_report

  type, public :: backtest_result
    real(dp), allocatable :: portfolio_returns(:)
    real(dp), allocatable :: wealth(:)
    real(dp), allocatable :: turnover(:)
    real(dp), allocatable :: transaction_costs(:)
    real(dp), allocatable :: weights(:,:)
    real(dp) :: total_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: max_drawdown = 0.0_dp
  end type backtest_result
end module fportfolio_types
