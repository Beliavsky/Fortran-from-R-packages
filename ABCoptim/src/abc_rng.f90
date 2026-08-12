module abc_rng
  use iso_fortran_env, only : int64
  implicit none
  private

  integer(int64), parameter :: pm_mod = 2147483647_int64
  integer(int64), parameter :: pm_mul = 16807_int64

  type, public :: abc_rng_state
    integer(int64) :: state = 123456789_int64
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: randint => rng_randint
  end type abc_rng_state

contains

  subroutine rng_seed(self, seed_value)
    class(abc_rng_state), intent(inout) :: self
    integer(int64), intent(in) :: seed_value

    self%state = modulo(abs(seed_value), pm_mod - 1_int64) + 1_int64
  end subroutine rng_seed

  function rng_uniform(self) result(u)
    class(abc_rng_state), intent(inout) :: self
    real(kind(1.0d0)) :: u

    self%state = modulo(pm_mul * self%state, pm_mod)
    u = real(self%state, kind(1.0d0)) / real(pm_mod, kind(1.0d0))
  end function rng_uniform

  function rng_randint(self, n) result(k)
    class(abc_rng_state), intent(inout) :: self
    integer, intent(in) :: n
    integer :: k
    real(kind(1.0d0)) :: u

    if (n <= 0) error stop "abc_rng: randint requires n > 0"
    u = self%uniform()
    k = 1 + min(n - 1, int(u * real(n, kind(1.0d0))))
  end function rng_randint

end module abc_rng
