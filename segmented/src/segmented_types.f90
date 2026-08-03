! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_types
  use nlme_kinds, only : dp
  use nlme_types, only : lme_result, nlme_control, correlation_spec, variance_spec, pd_spec
  use segmented_status, only : SEG_SUCCESS
  implicit none
  private

  integer, parameter, public :: SEGMENTED_CONTINUOUS = 1
  integer, parameter, public :: SEGMENTED_STEP = 2
  integer, parameter, public :: FAMILY_GAUSSIAN = 1
  integer, parameter, public :: FAMILY_BINOMIAL = 2
  integer, parameter, public :: FAMILY_POISSON = 3

  type, public :: segmented_control
    integer :: max_iter = 40
    integer :: glm_max_iter = 100
    integer :: grid_points = 256
    integer :: max_line_search = 14
    real(dp) :: tolerance = 1.0e-7_dp
    real(dp) :: breakpoint_tolerance = 1.0e-6_dp
    real(dp) :: lower_quantile = 0.05_dp
    real(dp) :: upper_quantile = 0.95_dp
    real(dp) :: min_segment_fraction = 0.03_dp
    logical :: verbose = .false.
  end type segmented_control

  type, public :: segmented_result
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: breakpoints(:)
    real(dp), allocatable :: breakpoint_se(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: history_breakpoints(:,:)
    real(dp), allocatable :: history_objective(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: sigma = 0.0_dp
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = SEG_SUCCESS
    integer :: family = FAMILY_GAUSSIAN
    integer :: kind = SEGMENTED_CONTINUOUS
    integer :: n_base = 0
    integer :: n_break = 0
    logical :: converged = .false.
  end type segmented_result

  type, public :: segmented_lme_result
    type(lme_result) :: fit
    real(dp), allocatable :: breakpoints(:)
    real(dp), allocatable :: breakpoint_se(:)
    real(dp), allocatable :: history_breakpoints(:,:)
    real(dp), allocatable :: history_objective(:)
    integer :: iterations = 0
    integer :: status = SEG_SUCCESS
    logical :: converged = .false.
  end type segmented_lme_result

  type, public :: segmented_lme_options
    type(nlme_control) :: control
    type(correlation_spec) :: correlation
    type(variance_spec) :: variance
    type(pd_spec) :: random
  end type segmented_lme_options

  type, public :: test_result
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: breakpoint = 0.0_dp
    integer :: grid_evaluated = 0
    integer :: status = SEG_SUCCESS
  end type test_result
end module segmented_types
