module rmalschains_types
  use iso_fortran_env, only : int64
  use rmalschains_kinds, only : dp
  implicit none
  private
  public :: mals_control, mals_result, mals_objective, ls_state

  abstract interface
    function mals_objective(x) result(value)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function mals_objective
  end interface

  type :: mals_control
    integer :: popsize = 50
    character(len=24) :: ls = 'cmaes'
    integer :: istep = 500
    real(dp) :: effort = 0.5_dp
    real(dp) :: alpha = 0.5_dp
    real(dp) :: optimum = -huge(1.0_dp)
    real(dp) :: threshold = 1.0e-8_dp
    logical :: ls_only = .false.
    integer :: ls_param1 = 0
    integer :: ls_param2 = 0
    integer :: verbosity = 0
    integer(int64) :: seed = 1_int64
    logical :: legacy_ls_only_zero_start = .true.
  end type mals_control

  type :: mals_result
    real(dp), allocatable :: sol(:)
    real(dp) :: fitness = huge(1.0_dp)
    integer :: num_eval_ea = 0
    integer :: num_eval_ls = 0
    integer :: actual_nfe = 0
    integer :: num_improvement_ea = 0
    integer :: num_improvement_ls = 0
    integer :: num_total_ea = 0
    integer :: num_total_ls = 0
    real(dp) :: improvement_ea = 0.0_dp
    real(dp) :: improvement_ls = 0.0_dp
    real(dp) :: time_ms_ea = 0.0_dp
    real(dp) :: time_ms_ls = 0.0_dp
    real(dp) :: time_ms_ma = 0.0_dp
    integer :: generations_ea = 0
    integer :: local_search_calls = 0
    character(len=:), allocatable :: stop_reason
  end type mals_result

  type :: ls_state
    logical :: initialized = .false.
    integer :: kind = 0
    real(dp) :: delta = 0.0_dp
    real(dp), allocatable :: delta_vec(:), delta_init(:), bias(:)
    integer :: num_failed = 0, num_success = 0
    real(dp), allocatable :: simplex(:, :), simplex_fit(:)
    real(dp), allocatable :: xmean(:), pc(:), ps(:), c(:, :), b(:, :), d(:), bd(:, :), weights(:)
    real(dp) :: sigma = 1.0_dp
    integer :: lambda = 0, mu = 0, cma_evals = 0
  end type ls_state
end module rmalschains_types
