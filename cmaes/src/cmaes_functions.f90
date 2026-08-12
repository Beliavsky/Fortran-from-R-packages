! cmaes-fortran - GPL-2.0-only
module cmaes_functions
  use iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use cmaes_kinds, only : dp
  use cmaes_rng, only : rng_state, rng_seed, rng_uniform
  implicit none
  private
  public :: f_sphere, f_rand, seed_f_rand, f_rosenbrock, f_rastrigin
  public :: shifted_value, rotated_value, biased_value

  type(rng_state), save :: helper_rng
  logical, save :: helper_rng_ready = .false.

  abstract interface
    function scalar_objective(x) result(value)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function scalar_objective
  end interface
contains
  subroutine seed_f_rand(seed)
    integer(int64), intent(in) :: seed
    call rng_seed(helper_rng, seed)
    helper_rng_ready = .true.
  end subroutine seed_f_rand

  function f_rand(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (.not. helper_rng_ready) call seed_f_rand(1_int64)
    if (size(x) < 0) error stop "f_rand: unreachable"
    value = rng_uniform(helper_rng)
  end function f_rand

  function f_sphere(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = dot_product(x, x)
  end function f_sphere

  function f_rosenbrock(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp), allocatable :: z(:)
    integer :: j

    if (size(x) < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    allocate(z(size(x)))
    z = x + 1.0_dp
    value = 0.0_dp
    do j = 1, size(x) - 1
      value = value + 100.0_dp * (z(j)**2 - z(j + 1))**2 + (z(j) - 1.0_dp)**2
    end do
  end function f_rosenbrock

  function f_rastrigin(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)
    value = sum(x * x - 10.0_dp * cos(twopi * x) + 10.0_dp)
  end function f_rastrigin

  function shifted_value(f, x, offset) result(value)
    procedure(scalar_objective) :: f
    real(dp), intent(in) :: x(:), offset(:)
    real(dp) :: value
    if (size(x) /= size(offset)) error stop "shifted_value: size mismatch"
    value = f(x - offset)
  end function shifted_value

  function rotated_value(f, x, m) result(value)
    procedure(scalar_objective) :: f
    real(dp), intent(in) :: x(:), m(:, :)
    real(dp) :: value
    if (size(m, 1) /= size(x) .or. size(m, 2) /= size(x)) &
      error stop "rotated_value: size mismatch"
    ! R's rotate_function stores t(M), then evaluates that matrix times x.
    value = f(matmul(transpose(m), x))
  end function rotated_value

  function biased_value(f, x, bias) result(value)
    procedure(scalar_objective) :: f
    real(dp), intent(in) :: x(:), bias
    real(dp) :: value
    value = f(x) + bias
  end function biased_value
end module cmaes_functions
