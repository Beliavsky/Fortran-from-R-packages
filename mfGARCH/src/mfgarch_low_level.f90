! SPDX-License-Identifier: MIT
module mfgarch_low_level
  use mfgarch_kinds, only : dp
  use mfgarch_math, only : finite_value
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument, &
    mfgarch_dimension_error, mfgarch_numerical_error
  implicit none
  private

  public :: sum_tau, sum_tau_fcts
  public :: calculate_h_andersen, calculate_p, simulate_r

contains

  subroutine sum_tau(m, theta, phivar, covariate, k, exponential, status)
    real(dp), intent(in) :: m, theta
    real(dp), intent(in) :: phivar(:), covariate(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: exponential(:)
    integer, intent(out) :: status
    integer :: i, j, n

    status = mfgarch_success
    if (k < 0 .or. size(phivar) /= k .or. size(covariate) <= k) then
      status = mfgarch_dimension_error
      allocate(exponential(0))
      return
    end if
    n = size(covariate) - k
    allocate(exponential(n))
    do i = 1, n
      exponential(i) = m
      do j = 1, k
        exponential(i) = exponential(i) + theta * phivar(j) * covariate(k+i-j)
      end do
    end do
  end subroutine sum_tau

  function sum_tau_fcts(m, theta, phivar, covariate, k, status) result(exponential)
    real(dp), intent(in) :: m, theta
    real(dp), intent(in) :: phivar(:), covariate(:)
    integer, intent(in) :: k
    integer, intent(out) :: status
    real(dp) :: exponential
    integer :: j

    status = mfgarch_success
    exponential = m
    if (k < 0 .or. size(phivar) /= k .or. size(covariate) < k) then
      status = mfgarch_dimension_error
      exponential = 0.0_dp
      return
    end if
    do j = 1, k
      exponential = exponential + theta * phivar(j) * covariate(size(covariate)+1-j)
    end do
  end function sum_tau_fcts

  subroutine calculate_h_andersen(n_days, delta, theta, omega, lambda, z, h0, h, status)
    integer, intent(in) :: n_days, delta
    real(dp), intent(in) :: theta, omega, lambda, z(:), h0
    real(dp), allocatable, intent(out) :: h(:)
    integer, intent(out) :: status
    integer :: i, n

    status = mfgarch_success
    if (n_days <= 0 .or. delta <= 0 .or. theta <= 0.0_dp .or. omega <= 0.0_dp .or. &
        lambda <= 0.0_dp .or. h0 <= 0.0_dp) then
      status = mfgarch_invalid_argument
      allocate(h(0))
      return
    end if
    n = n_days * delta
    if (size(z) < n - 1) then
      status = mfgarch_dimension_error
      allocate(h(0))
      return
    end if
    allocate(h(n))
    h(1) = h0
    do i = 2, n
      h(i) = theta * omega / real(delta,dp) + h(i-1) * &
        (1.0_dp - theta / real(delta,dp) + &
        sqrt(2.0_dp*lambda*theta/real(delta,dp))*z(i-1))
      if (.not. finite_value(h(i))) then
        status = mfgarch_numerical_error
        return
      end if
    end do
  end subroutine calculate_h_andersen

  subroutine calculate_p(n_days, delta, zp, h, p0, p, status)
    integer, intent(in) :: n_days, delta
    real(dp), intent(in) :: zp(:), h(:), p0
    real(dp), allocatable, intent(out) :: p(:)
    integer, intent(out) :: status
    integer :: i, n

    status = mfgarch_success
    if (n_days <= 0 .or. delta <= 0) then
      status = mfgarch_invalid_argument
      allocate(p(0))
      return
    end if
    n = n_days * delta
    if (size(zp) < n .or. size(h) < n) then
      status = mfgarch_dimension_error
      allocate(p(0))
      return
    end if
    allocate(p(n))
    p(1) = p0
    do i = 2, n
      if (h(i) < 0.0_dp) then
        status = mfgarch_numerical_error
        return
      end if
      p(i) = p(i-1) + sqrt(h(i)/real(delta,dp)) * zp(i)
    end do
  end subroutine calculate_p

  subroutine simulate_r(n_days, n_intraday, alpha, beta, gamma, z, h0, ret_daily, &
      h_daily, ret_intraday, realized_variance, status)
    integer, intent(in) :: n_days, n_intraday
    real(dp), intent(in) :: alpha, beta, gamma, z(:), h0
    real(dp), allocatable, intent(out) :: ret_daily(:), h_daily(:), ret_intraday(:)
    real(dp), allocatable, intent(out) :: realized_variance(:)
    integer, intent(out) :: status
    real(dp) :: omega, shock
    integer :: day, first, last

    status = mfgarch_success
    omega = 1.0_dp - alpha - beta - 0.5_dp*gamma
    if (n_days <= 0 .or. n_intraday <= 0 .or. size(z) < n_days*n_intraday .or. &
        alpha < 0.0_dp .or. beta < 0.0_dp .or. gamma < 0.0_dp .or. omega <= 0.0_dp .or. h0 <= 0.0_dp) then
      status = mfgarch_invalid_argument
      allocate(ret_daily(0), h_daily(0), ret_intraday(0), realized_variance(0))
      return
    end if
    allocate(ret_daily(n_days), h_daily(n_days), ret_intraday(n_days*n_intraday), &
      realized_variance(n_days))
    h_daily(1) = h0
    do day = 1, n_days
      first = (day-1)*n_intraday + 1
      last = day*n_intraday
      ret_intraday(first:last) = z(first:last) * sqrt(h_daily(day)/real(n_intraday,dp))
      ret_daily(day) = sum(ret_intraday(first:last))
      realized_variance(day) = sum(ret_intraday(first:last)**2)
      if (day < n_days) then
        shock = alpha*ret_daily(day)**2
        if (ret_daily(day) < 0.0_dp) shock = shock + gamma*ret_daily(day)**2
        h_daily(day+1) = omega + shock + beta*h_daily(day)
        if (.not. finite_value(h_daily(day+1)) .or. h_daily(day+1) <= 0.0_dp) then
          status = mfgarch_numerical_error
          return
        end if
      end if
    end do
  end subroutine simulate_r

end module mfgarch_low_level
