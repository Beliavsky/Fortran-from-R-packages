module rmalschains_rng
  use iso_fortran_env, only : int64
  use rmalschains_kinds, only : dp
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_int

  integer(int64), parameter :: modulus = 2147483647_int64
  integer(int64), parameter :: multiplier = 16807_int64

  type :: rng_state
    integer(int64) :: state = 1_int64
    logical :: have_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(int64), intent(in) :: seed
    rng%state = modulo(abs(seed), modulus - 1_int64) + 1_int64
    rng%have_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    rng%state = modulo(multiplier * rng%state, modulus)
    if (rng%state <= 0_int64) rng%state = 1_int64
    u = real(rng%state, dp) / real(modulus, dp)
  end function rng_uniform

  function rng_int(rng, lo, hi) result(k)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: lo, hi
    integer :: k
    if (hi < lo) error stop "rng_int: invalid interval"
    k = lo + int(rng_uniform(rng) * real(hi - lo + 1, dp))
    if (k > hi) k = hi
  end function rng_int

  function rng_normal(rng, sd) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: sd
    real(dp) :: z, u1, u2, scale
    real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp
    scale = 1.0_dp
    if (present(sd)) scale = sd
    if (rng%have_spare) then
      z = rng%spare
      rng%have_spare = .false.
    else
      u1 = max(rng_uniform(rng), tiny(1.0_dp))
      u2 = rng_uniform(rng)
      z = sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
      rng%spare = sqrt(-2.0_dp * log(u1)) * sin(twopi * u2)
      rng%have_spare = .true.
    end if
    z = scale * z
  end function rng_normal
end module rmalschains_rng
