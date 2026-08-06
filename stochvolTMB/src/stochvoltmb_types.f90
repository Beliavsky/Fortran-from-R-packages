! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_types
  use stochvoltmb_kinds, only : dp
  use stochvoltmb_status, only : sv_ok
  implicit none
  private

  integer, parameter, public :: sv_gaussian = 0
  integer, parameter, public :: sv_student_t = 1
  integer, parameter, public :: sv_skew_gaussian = 2
  integer, parameter, public :: sv_leverage = 3
  integer, parameter, public :: sv_skew_gaussian_leverage = 4

  type, public :: sv_parameters
    real(dp) :: sigma_y = 1.0_dp
    real(dp) :: sigma_h = 0.2_dp
    real(dp) :: phi = 0.95_dp
    real(dp) :: df = 8.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: rho = 0.0_dp
  end type sv_parameters

  type, public :: sv_control
    integer :: max_outer_iter = 250
    integer :: max_inner_iter = 80
    real(dp) :: outer_tolerance = 1.0e-5_dp
    real(dp) :: inner_tolerance = 1.0e-7_dp
    logical :: compute_covariance = .true.
    logical :: verbose = .false.
  end type sv_control

  type, public :: sv_fit_result
    type(sv_parameters) :: params
    integer :: model = sv_gaussian
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: theta_cov(:,:)
    real(dp), allocatable :: theta_se(:)
    real(dp), allocatable :: param_se(:)
    real(dp), allocatable :: h(:)
    real(dp), allocatable :: h_se(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: nobs = 0
    integer :: iterations = 0
    integer :: function_evaluations = 0
    logical :: converged = .false.
    integer :: status = sv_ok
    character(len=160) :: message = ''
  end type sv_fit_result

  type, public :: sv_simulation
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: h(:)
    integer :: model = sv_gaussian
    type(sv_parameters) :: params
    integer :: status = sv_ok
    character(len=160) :: message = ''
  end type sv_simulation

  type, public :: sv_prediction
    real(dp), allocatable :: y(:,:)
    real(dp), allocatable :: h(:,:)
    real(dp), allocatable :: h_exp(:,:)
    integer :: steps = 0
    integer :: nsim = 0
    integer :: status = sv_ok
    character(len=160) :: message = ''
  end type sv_prediction

  type, public :: sv_prediction_summary
    real(dp), allocatable :: probabilities(:)
    real(dp), allocatable :: y_quantiles(:,:)
    real(dp), allocatable :: h_quantiles(:,:)
    real(dp), allocatable :: h_exp_quantiles(:,:)
    real(dp), allocatable :: y_mean(:)
    real(dp), allocatable :: h_mean(:)
    real(dp), allocatable :: h_exp_mean(:)
    integer :: status = sv_ok
    character(len=160) :: message = ''
  end type sv_prediction_summary

end module stochvoltmb_types
