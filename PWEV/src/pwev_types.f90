! SPDX-License-Identifier: GPL-3.0-only
module pwev_types
  use pwev_kinds, only : dp
  use pwev_status, only : PWEV_INVALID_INPUT, PWEV_GARCH_MEAN, PWEV_MEM_UPSTREAM_OOS
  implicit none
  private

  integer, parameter, public :: PWEV_N_BASE_MODELS = 4
  integer, parameter, public :: PWEV_N_MODELS = 5
  integer, parameter, public :: PWEV_N_METRICS = 9

  type, public :: pwev_control
    integer :: garch_output = PWEV_GARCH_MEAN
    integer :: mem_oos_mode = PWEV_MEM_UPSTREAM_OOS
    integer :: garch_max_iterations = 1800
    integer :: mem_max_iterations = 400
    integer :: mem_random_starts = 12
    integer :: pso_iterations = 1000
    integer :: pso_population = 0
    integer :: random_seed = 12345
    real(dp) :: pso_vmax = 2.0_dp
    real(dp) :: pso_individual = 1.49445_dp
    real(dp) :: pso_group = 1.49445_dp
    real(dp) :: pso_inertia = 0.729_dp
    logical :: round_accuracy = .true.
    logical :: fail_on_base_model_error = .false.
  end type pwev_control

  type, public :: pwev_result
    real(dp), allocatable :: train_fitted(:, :)
    real(dp), allocatable :: test_pred(:, :)
    real(dp), allocatable :: accuracy(:, :)
    real(dp), allocatable :: weights(:)
    integer :: train_size = 0
    integer :: test_size = 0
    integer :: status = PWEV_INVALID_INPUT
    integer :: model_status(PWEV_N_BASE_MODELS) = PWEV_INVALID_INPUT
    integer :: pso_iterations = 0
    real(dp) :: pso_objective = huge(1.0_dp)
    character(len=160) :: message = 'not run'
  end type pwev_result

  type, public :: pwev_pso_result
    real(dp), allocatable :: weights(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = PWEV_INVALID_INPUT
  end type pwev_pso_result

end module pwev_types
