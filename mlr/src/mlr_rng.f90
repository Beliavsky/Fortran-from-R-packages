module mlr_rng
  use mlr_kinds, only : dp, i8
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_integer, rng_shuffle, rng_normal
  type :: rng_state
    integer(i8) :: state = 123456789_i8
  end type rng_state
contains
  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    rng%state = modulo(abs(seed), 2147483646_i8) + 1_i8
  end subroutine rng_seed

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(i8) :: hi, lo, test
    hi = rng%state / 127773_i8
    lo = modulo(rng%state, 127773_i8)
    test = 16807_i8 * lo - 2836_i8 * hi
    if (test <= 0_i8) test = test + 2147483647_i8
    rng%state = test
    u = real(test, dp) / 2147483647.0_dp
  end function rng_uniform

  integer function rng_integer(rng, lo, hi) result(v)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: lo, hi
    if (hi < lo) error stop "rng_integer: invalid interval"
    v = lo + min(hi-lo, int(rng_uniform(rng) * real(hi-lo+1,dp)))
  end function rng_integer

  subroutine rng_shuffle(rng, x)
    type(rng_state), intent(inout) :: rng
    integer, intent(inout) :: x(:)
    integer :: i, j, t
    do i = size(x), 2, -1
      j = rng_integer(rng, 1, i)
      t = x(i); x(i) = x(j); x(j) = t
    end do
  end subroutine rng_shuffle

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2
    u1 = max(rng_uniform(rng), tiny(1.0_dp))
    u2 = rng_uniform(rng)
    z = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function rng_normal
end module mlr_rng
