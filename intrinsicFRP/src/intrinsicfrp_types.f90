! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_types
  use intrinsicfrp_kinds, only: dp, status_ok
  implicit none
  private

  type, public :: vector_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: standard_errors(:)
    integer, allocatable :: selected_indices(:)
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type vector_result

  type, public :: screening_result
    real(dp), allocatable :: sdf_coefficients(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: t_statistics(:)
    integer, allocatable :: selected_indices(:)
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type screening_result

  type, public :: hj_result
    real(dp) :: squared_distance = 0.0_dp
    real(dp) :: lower_bound = 0.0_dp
    real(dp) :: upper_bound = 0.0_dp
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type hj_result

  type, public :: rank_test_result
    integer :: rank = 0
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp), allocatable :: statistics(:)
    real(dp), allocatable :: p_values(:)
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type rank_test_result

  type, public :: pca_result
    real(dp), allocatable :: risk_premia(:)
    integer :: n_pca = 0
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type pca_result

  type, public :: oracle_control
    character(len=1) :: weighting_type = 'c'
    character(len=1) :: tuning_type = 'g'
    logical :: one_stddev_rule = .false.
    logical :: gcv_scaling_n_assets = .false.
    logical :: gcv_identification_check = .false.
    real(dp) :: target_level_kp2006 = 0.05_dp
    integer :: n_folds = 5
    integer :: n_train_observations = 120
    integer :: n_test_observations = 12
    integer :: roll_shift = 12
    logical :: relaxed = .false.
    logical :: include_standard_errors = .false.
    logical :: hac_prewhite = .false.
  end type oracle_control

  type, public :: oracle_result
    real(dp), allocatable :: risk_premia(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: model_score(:)
    real(dp) :: penalty_parameter = 0.0_dp
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type oracle_result

  type, public :: fgx_result
    real(dp), allocatable :: sdf_coefficients(:)
    real(dp), allocatable :: standard_errors(:)
    integer, allocatable :: controls_selected(:)
    integer :: status = status_ok
    character(len=160) :: message = ''
  end type fgx_result
end module intrinsicfrp_types
