module dice_design_rng
  use iso_fortran_env, only : int64
  use dice_design_kinds, only : dp
  implicit none
  private

  integer(int64), parameter :: pm_m = 2147483647_int64
  integer(int64), parameter :: pm_a = 16807_int64
  integer(int64), parameter :: pm_q = 127773_int64
  integer(int64), parameter :: pm_r = 2836_int64

  type, public :: rng_state
    integer(int64) :: state = 1_int64
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: integer => rng_integer
    procedure :: shuffle => rng_shuffle
  end type rng_state

contains

  subroutine rng_seed(self, seed)
    class(rng_state), intent(inout) :: self
    integer(int64), intent(in) :: seed
    integer(int64) :: s

    s = modulo(seed, pm_m - 1_int64)
    if (s < 0_int64) s = s + pm_m - 1_int64
    self%state = s + 1_int64
  end subroutine rng_seed

  function rng_uniform(self) result(u)
    class(rng_state), intent(inout) :: self
    real(dp) :: u
    integer(int64) :: hi, lo, test

    hi = self%state / pm_q
    lo = modulo(self%state, pm_q)
    test = pm_a * lo - pm_r * hi
    if (test > 0_int64) then
      self%state = test
    else
      self%state = test + pm_m
    end if
    u = real(self%state, dp) / real(pm_m, dp)
  end function rng_uniform

  function rng_integer(self, lo, hi) result(k)
    class(rng_state), intent(inout) :: self
    integer, intent(in) :: lo, hi
    integer :: k
    real(dp) :: u

    if (hi < lo) then
      k = lo
      return
    end if
    u = self%uniform()
    k = lo + int(u * real(hi - lo + 1, dp))
    if (k > hi) k = hi
  end function rng_integer

  subroutine rng_shuffle(self, x)
    class(rng_state), intent(inout) :: self
    integer, intent(inout) :: x(:)
    integer :: i, j, tmp

    do i = size(x), 2, -1
      j = self%integer(1, i)
      tmp = x(i)
      x(i) = x(j)
      x(j) = tmp
    end do
  end subroutine rng_shuffle

end module dice_design_rng
