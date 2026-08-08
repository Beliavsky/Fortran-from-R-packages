! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_types
  use rgenoud_kinds, only : dp
  implicit none
  private

  type, public :: genoud_options
    integer :: pop_size = 1000
    integer :: max_generations = 100
    integer :: wait_generations = 10
    logical :: hard_generation_limit = .true.
    logical :: maximize = .false.
    logical :: memory_matrix = .true.
    integer :: max_memory_evals = 1000000
    real(dp) :: solution_tolerance = 1.0e-3_dp
    logical :: gradient_check = .true.
    logical :: use_bfgs = .true.
    integer :: boundary_enforcement = 0
    logical :: integer_parameters = .false.
    real(dp) :: operator_weights(9) = [50.0_dp, 50.0_dp, 50.0_dp, &
      50.0_dp, 50.0_dp, 50.0_dp, 50.0_dp, 50.0_dp, 0.0_dp]
    real(dp) :: p9_mix = -1.0_dp
    integer :: bfgs_burnin = 0
    integer :: bfgs_max_iter = 200
    real(dp) :: bfgs_gtol = 1.0e-8_dp
    integer :: seed = 53058
    integer :: print_level = 0
  end type genoud_options

  type, public :: genoud_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: fit(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    integer :: generations = 0
    integer :: peak_generation = 0
    integer :: pop_size = 0
    integer :: operators(9) = 0
    integer :: evaluations = 0
    integer :: unique_evaluations = 0
    logical :: converged = .false.
    integer :: status = 0
  end type genoud_result

  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn

    subroutine gradient_fn(x, g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_fn

    subroutine lexical_objective_fn(x, f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine lexical_objective_fn

    logical function lexical_better_fn(a, b)
      import dp
      real(dp), intent(in) :: a(:), b(:)
    end function lexical_better_fn
  end interface

  public :: objective_fn, gradient_fn, lexical_objective_fn, lexical_better_fn
end module rgenoud_types
