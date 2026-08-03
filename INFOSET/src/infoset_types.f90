! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_types
  use infoset_kinds, only : dp
  use infoset_status, only : infoset_success
  implicit none
  private
  public :: mixture_control, tail_mixture_result, information_set_result
  public :: window_collection, left_risk_result, portfolio_result, portfolio_summary

  type :: mixture_control
    real(dp) :: tolerance = 1.0e-8_dp
    integer :: max_iterations = 1000
    real(dp) :: minimum_scale = 1.0e-6_dp
  end type mixture_control

  type :: tail_mixture_result
    real(dp) :: change_point = 0.0_dp
    integer :: flag = 1
    real(dp) :: left_mean = 0.0_dp
    real(dp) :: left_sd = 0.0_dp
    real(dp) :: right_mean = 0.0_dp
    real(dp) :: right_sd = 0.0_dp
    real(dp) :: left_probability = 1.0_dp
    real(dp) :: first_type_error = 1.0_dp
    real(dp) :: second_type_error = 1.0_dp
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: iterations = 0
    logical :: converged = .false.
    integer :: status = infoset_success
  end type tail_mixture_result

  type :: information_set_result
    integer :: n_change_points = 0
    real(dp) :: change_points(2) = 0.0_dp
    real(dp) :: prior_probability(2) = 0.0_dp
    real(dp) :: first_type_error(2) = 0.0_dp
    real(dp) :: second_type_error(2) = 0.0_dp
    real(dp) :: left_mean(2) = 0.0_dp
    real(dp) :: left_sd(2) = 0.0_dp
    integer :: status = infoset_success
  end type information_set_result

  type :: window_collection
    real(dp), allocatable :: values(:,:,:)
    integer :: window_size = 0
    integer :: step = 0
    integer :: status = infoset_success
  end type window_collection

  type :: left_risk_result
    real(dp), allocatable :: values(:,:)
    real(dp), allocatable :: first_change_point(:)
    integer :: status = infoset_success
  end type left_risk_result

  type :: portfolio_result
    real(dp), allocatable :: weights(:,:)
    real(dp), allocatable :: oos_returns(:,:)
    real(dp), allocatable :: baseline_value(:)
    character(len=8) :: strategy = 'M'
    integer :: status = infoset_success
  end type portfolio_result

  type :: portfolio_summary
    integer :: count = 0
    real(dp) :: minimum = 0.0_dp
    real(dp) :: first_quartile = 0.0_dp
    real(dp) :: median = 0.0_dp
    real(dp) :: mean = 0.0_dp
    real(dp) :: third_quartile = 0.0_dp
    real(dp) :: maximum = 0.0_dp
    integer :: status = infoset_success
  end type portfolio_summary
end module infoset_types
