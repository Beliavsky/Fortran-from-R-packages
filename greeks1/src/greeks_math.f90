! SPDX-License-Identifier: MIT
module greeks_math
  use greeks_kinds, only: dp
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)

  type, public :: rng_state
    integer(kind=8) :: state = 1_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

  public :: normal_pdf, normal_cdf, rng_seed, rng_uniform, rng_normal
  public :: rng_poisson, rng_student_t3, sample_mean, sample_stderr
  public :: regression_intercept

contains

  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(out) :: rng
    integer, intent(in) :: seed
    rng%state = int(max(1, abs(seed)), kind=8)
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  real(dp) function rng_uniform(rng) result(value)
    type(rng_state), intent(inout) :: rng
    integer(kind=8), parameter :: a = 48271_8
    integer(kind=8), parameter :: m = 2147483647_8
    integer(kind=8), parameter :: q = 44488_8
    integer(kind=8), parameter :: r = 3399_8
    integer(kind=8) :: hi, lo, test

    hi = rng%state / q
    lo = modulo(rng%state, q)
    test = a*lo - r*hi
    if (test > 0_8) then
      rng%state = test
    else
      rng%state = test + m
    end if
    value = real(rng%state, dp) / real(m, dp)
  end function rng_uniform

  real(dp) function rng_normal(rng) result(value)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2, radius

    if (rng%has_spare) then
      value = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = max(rng_uniform(rng), tiny(1.0_dp))
    u2 = rng_uniform(rng)
    radius = sqrt(-2.0_dp*log(u1))
    value = radius*cos(2.0_dp*pi*u2)
    rng%spare = radius*sin(2.0_dp*pi*u2)
    rng%has_spare = .true.
  end function rng_normal

  integer function rng_poisson(rng, lambda) result(value)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: lambda
    real(dp) :: product_value, limit_value, z
    integer :: k

    if (lambda <= 0.0_dp) then
      value = 0
    else if (lambda < 30.0_dp) then
      limit_value = exp(-lambda)
      product_value = 1.0_dp
      k = 0
      do
        k = k + 1
        product_value = product_value*rng_uniform(rng)
        if (product_value <= limit_value) exit
      end do
      value = k - 1
    else
      do
        z = lambda + sqrt(lambda)*rng_normal(rng)
        if (z >= 0.0_dp) exit
      end do
      value = nint(z)
    end if
  end function rng_poisson

  real(dp) function rng_student_t3(rng) result(value)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z, chi2
    integer :: i
    z = rng_normal(rng)
    chi2 = 0.0_dp
    do i = 1, 3
      chi2 = chi2 + rng_normal(rng)**2
    end do
    value = z/sqrt(max(chi2/3.0_dp, tiny(1.0_dp)))
  end function rng_student_t3

  pure real(dp) function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x)/real(size(x), dp)
    end if
  end function sample_mean

  pure real(dp) function sample_stderr(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_value
    integer :: n
    n = size(x)
    if (n <= 1) then
      value = 0.0_dp
      return
    end if
    mean_value = sample_mean(x)
    value = sqrt(sum((x - mean_value)**2)/real(n - 1, dp)/real(n, dp))
  end function sample_stderr

  pure real(dp) function regression_intercept(y, x) result(value)
    real(dp), intent(in) :: y(:), x(:)
    real(dp) :: mx, my, denom, beta
    mx = sample_mean(x)
    my = sample_mean(y)
    denom = sum((x - mx)**2)
    if (denom <= epsilon(1.0_dp)) then
      value = my
    else
      beta = sum((x - mx)*(y - my))/denom
      value = my - beta*mx
    end if
  end function regression_intercept

end module greeks_math
