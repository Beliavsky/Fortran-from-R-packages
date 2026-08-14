! SPDX-License-Identifier: GPL-3.0-only
module anmc_types
  use anmc_kinds, only : dp
  implicit none
  private

  type, public :: anmc_problem
    real(dp), allocatable :: mu_eq(:)
    real(dp), allocatable :: sigma_eq(:,:)
    real(dp) :: threshold = 0.0_dp
    real(dp), allocatable :: mu_emq(:)
    real(dp), allocatable :: ww_cond_q(:,:)
    ! Upper-triangular Cholesky factor U, covariance = transpose(U) * U.
    real(dp), allocatable :: sigma_cond_q_chol(:,:)
  end type anmc_problem

  type, public :: simulation_control
    integer :: max_outer = 200000
    integer :: max_inner = 10000
    integer :: max_rejection_batch = 1000000
    integer :: max_rejection_draws = 100000000
    real(dp) :: time_guard = 0.96_dp
    logical :: enforce_budget = .true.
  end type simulation_control

  type, public :: mc_params
    integer :: n = 0
    integer :: m = 1
    real(dp) :: cx = 0.0_dp
    real(dp) :: cx0 = 0.0_dp
    real(dp) :: alpha_cost = 0.0_dp
    real(dp) :: beta = 0.0_dp
    real(dp) :: eval_g = 0.0_dp
  end type mc_params

  type, public :: mc_result
    real(dp) :: estim = 0.0_dp
    real(dp) :: var_est = 0.0_dp
    type(mc_params) :: params
    integer :: n0 = 0
    integer :: m0 = 0
    real(dp) :: estim0 = 0.0_dp
    real(dp) :: var_est0 = 0.0_dp
    real(dp) :: ratio0 = 0.0_dp
    real(dp) :: time_part1 = 0.0_dp
    real(dp) :: time_total = 0.0_dp
    real(dp), allocatable :: exp_y_cond_x(:)
    real(dp), allocatable :: var_y_cond_x(:)
    logical :: ok = .true.
    character(len=256) :: message = ''
  end type mc_result

  type, public :: active_dims_result
    integer, allocatable :: ind_q(:)
    real(dp), allocatable :: eq(:,:)
    real(dp), allocatable :: mu_eq(:)
    real(dp), allocatable :: k_eq(:,:)
    real(dp) :: pq = 0.0_dp
    real(dp) :: error = 0.0_dp
    logical :: ok = .true.
    character(len=256) :: message = ''
  end type active_dims_result

  type, public :: probability_estimate
    real(dp) :: probability = 0.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: pq = 0.0_dp
    real(dp) :: rq = 0.0_dp
    real(dp) :: pq_error = 0.0_dp
    integer :: q = 0
    integer, allocatable :: ind_q(:)
    real(dp), allocatable :: eq(:,:)
    type(mc_result) :: remainder
    logical :: has_remainder = .false.
    logical :: ok = .true.
    character(len=256) :: message = ''
  end type probability_estimate

  type, public :: conservative_result
    logical, allocatable :: set(:)
    real(dp) :: level = 0.0_dp
    real(dp) :: probability = 0.0_dp
    ! Matches upstream `vars`: MC variance in the large-set branch, but
    ! mvtnorm absolute integration error in the direct-Genz branch.
    real(dp) :: uncertainty = 0.0_dp
    logical :: ok = .true.
    character(len=256) :: message = ''
  end type conservative_result

end module anmc_types
