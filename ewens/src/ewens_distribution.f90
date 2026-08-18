! SPDX-License-Identifier: MIT
module ewens_distribution
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_quiet_nan
  use ewens_kinds, only : dp
  use ewens_math, only : class_size_spectrum, log_unsigned_stirling1, log_rising_factorial
  implicit none
  private
  public :: dewens, dewens_log, dewens_counts, dewens_counts_log
  public :: dewens_k, dewens_k_log, ewens_k_exact

contains

  real(dp) function dewens(labels, theta) result(probability)
    integer, intent(in) :: labels(:)
    real(dp), intent(in), optional :: theta
    real(dp) :: th, lp

    th = 1.0_dp
    if (present(theta)) th = theta
    lp = dewens_log(labels, th)
    probability = exp(lp)
  end function dewens

  real(dp) function dewens_log(labels, theta) result(log_probability)
    integer, intent(in) :: labels(:)
    real(dp), intent(in), optional :: theta
    integer, allocatable :: mvec(:)
    real(dp) :: th

    th = 1.0_dp
    if (present(theta)) th = theta
    call class_size_spectrum(labels, mvec)
    log_probability = dewens_counts_log(mvec, th)
  end function dewens_log

  real(dp) function dewens_counts(mvec, theta) result(probability)
    integer, intent(in) :: mvec(:)
    real(dp), intent(in), optional :: theta
    real(dp) :: th

    th = 1.0_dp
    if (present(theta)) th = theta
    probability = exp(dewens_counts_log(mvec, th))
  end function dewens_counts

  real(dp) function dewens_counts_log(mvec, theta) result(log_probability)
    integer, intent(in) :: mvec(:)
    real(dp), intent(in), optional :: theta
    real(dp) :: th, neginf
    integer :: j, n, k

    th = 1.0_dp
    if (present(theta)) th = theta
    neginf = ieee_value(0.0_dp, ieee_negative_inf)

    if (th < 0.0_dp .or. any(mvec < 0)) then
      log_probability = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    n = 0
    k = 0
    do j = 1, size(mvec)
      n = n + j * mvec(j)
      k = k + mvec(j)
    end do

    if (n == 0) then
      log_probability = 0.0_dp
      return
    end if

    if (.not. (th > 0.0_dp)) then
      log_probability = neginf
      if (k == 1) then
        if (size(mvec) >= n) then
          if (mvec(n) == 1) log_probability = 0.0_dp
        end if
      end if
      return
    end if

    log_probability = log_gamma(real(n + 1, dp)) - log_rising_factorial(th, n)
    do j = 1, size(mvec)
      if (mvec(j) == 0) cycle
      log_probability = log_probability + real(mvec(j), dp) * log(th)
      log_probability = log_probability - real(mvec(j), dp) * log(real(j, dp))
      log_probability = log_probability - log_gamma(real(mvec(j) + 1, dp))
    end do
  end function dewens_counts_log

  real(dp) function dewens_k(k, n, theta) result(probability)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: theta
    probability = exp(dewens_k_log(k, n, theta))
  end function dewens_k

  real(dp) function dewens_k_log(k, n, theta) result(log_probability)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: theta
    real(dp) :: neginf

    neginf = ieee_value(0.0_dp, ieee_negative_inf)
    if (theta < 0.0_dp .or. n < 0 .or. k < 0 .or. k > n) then
      log_probability = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (n == 0) then
      if (k == 0) then
        log_probability = 0.0_dp
      else
        log_probability = neginf
      end if
      return
    end if
    if (.not. (theta > 0.0_dp)) then
      if (k == 1) then
        log_probability = 0.0_dp
      else
        log_probability = neginf
      end if
      return
    end if
    if (k == 0) then
      log_probability = neginf
      return
    end if

    log_probability = log_unsigned_stirling1(n, k) + real(k, dp) * log(theta) &
      - log_rising_factorial(theta, n)
  end function dewens_k_log

  real(dp) function ewens_k_exact(n, theta) result(mean_k)
    integer, intent(in) :: n
    real(dp), intent(in) :: theta
    integer :: j

    if (n <= 0) then
      mean_k = 0.0_dp
      return
    end if
    if (theta < 0.0_dp) then
      mean_k = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (.not. (theta > 0.0_dp)) then
      mean_k = 1.0_dp
      return
    end if

    mean_k = 0.0_dp
    do j = 0, n - 1
      mean_k = mean_k + theta / (theta + real(j, dp))
    end do
  end function ewens_k_exact

end module ewens_distribution
