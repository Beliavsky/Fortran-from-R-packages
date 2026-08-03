! SPDX-License-Identifier: Apache-2.0
module intraday_kalman
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use intraday_kinds, only : dp
  use intraday_types, only : volume_parameters, parameter_mask, kalman_output, &
                             intraday_ok, intraday_invalid_input, intraday_numerical_failure
  use intraday_utils, only : parameters_valid
  implicit none
  private

  public :: uniss_kalman, uniss_em_update

contains

  subroutine uniss_kalman(log_data, par, output, smooth)
    real(dp), intent(in) :: log_data(:, :)
    type(volume_parameters), intent(in) :: par
    type(kalman_output), intent(out) :: output
    logical, intent(in), optional :: smooth

    integer :: n_bin, n_day, n, t
    real(dp), allocatable :: y(:)
    real(dp) :: a(2, 2), q(2, 2), s, innovation
    real(dp) :: cv(2), inv_pred(2, 2), tmp(2, 2)
    logical :: do_smooth

    output%status = intraday_ok
    output%message = ''
    n_bin = size(log_data, 1)
    n_day = size(log_data, 2)
    n = n_bin * n_day
    do_smooth = .false.
    if (present(smooth)) do_smooth = smooth

    if (n_bin < 1 .or. n_day < 1 .or. .not. parameters_valid(par, n_bin)) then
      output%status = intraday_invalid_input
      output%message = 'invalid data dimensions or model parameters'
      return
    end if

    allocate(y(n))
    y = reshape(log_data, [n])
    if (.not. all(ieee_is_finite(y))) then
      output%status = intraday_invalid_input
      output%message = 'log data must be finite'
      return
    end if

    allocate(output%x_pred(2, n), output%v_pred(2, 2, n))
    allocate(output%x_filt(2, n), output%v_filt(2, 2, n), output%gain(2, n))
    output%x_pred(:, 1) = par%x0
    output%v_pred(:, :, 1) = par%v0
    call measurement_update(1)

    do t = 1, n - 1
      if (mod(t, n_bin) == 0) then
        a = 0.0_dp
        a(1, 1) = par%a_eta
        a(2, 2) = par%a_mu
        q = 0.0_dp
        q(1, 1) = par%var_eta
        q(2, 2) = par%var_mu
      else
        a = 0.0_dp
        a(1, 1) = 1.0_dp
        a(2, 2) = par%a_mu
        q = 0.0_dp
        q(2, 2) = par%var_mu
      end if
      output%x_pred(:, t + 1) = matmul(a, output%x_filt(:, t))
      output%v_pred(:, :, t + 1) = matmul(a, matmul(output%v_filt(:, :, t), transpose(a))) + q
      call symmetrize(output%v_pred(:, :, t + 1))
      call measurement_update(t + 1)
      if (output%status /= intraday_ok) return
    end do

    if (.not. do_smooth) return

    allocate(output%x_smooth(2, n), output%v_smooth(2, 2, n))
    allocate(output%smoother_gain(2, 2, max(0, n - 1)))
    output%x_smooth(:, n) = output%x_filt(:, n)
    output%v_smooth(:, :, n) = output%v_filt(:, :, n)

    do t = n - 1, 1, -1
      if (mod(t, n_bin) == 0) then
        a = 0.0_dp
        a(1, 1) = par%a_eta
        a(2, 2) = par%a_mu
      else
        a = 0.0_dp
        a(1, 1) = 1.0_dp
        a(2, 2) = par%a_mu
      end if
      call inverse2_regularized(output%v_pred(:, :, t + 1), inv_pred)
      output%smoother_gain(:, :, t) = matmul(output%v_filt(:, :, t), &
                                             matmul(transpose(a), inv_pred))
      output%x_smooth(:, t) = output%x_filt(:, t) + &
        matmul(output%smoother_gain(:, :, t), output%x_smooth(:, t + 1) - output%x_pred(:, t + 1))
      tmp = output%v_smooth(:, :, t + 1) - output%v_pred(:, :, t + 1)
      output%v_smooth(:, :, t) = output%v_filt(:, :, t) + &
        matmul(output%smoother_gain(:, :, t), matmul(tmp, transpose(output%smoother_gain(:, :, t))))
      call symmetrize(output%v_smooth(:, :, t))
    end do

  contains

    subroutine measurement_update(index)
      integer, intent(in) :: index
      integer :: bin

      bin = mod(index - 1, n_bin) + 1
      cv(1) = output%v_pred(1, 1, index) + output%v_pred(2, 1, index)
      cv(2) = output%v_pred(1, 2, index) + output%v_pred(2, 2, index)
      s = cv(1) + cv(2) + par%r
      if (.not. ieee_is_finite(s) .or. s <= epsilon(1.0_dp)) then
        output%status = intraday_numerical_failure
        output%message = 'nonpositive innovation variance in Kalman filter'
        return
      end if
      output%gain(1, index) = &
        (output%v_pred(1, 1, index) + output%v_pred(1, 2, index)) / s
      output%gain(2, index) = &
        (output%v_pred(2, 1, index) + output%v_pred(2, 2, index)) / s
      innovation = y(index) - par%phi(bin) - sum(output%x_pred(:, index))
      output%x_filt(:, index) = output%x_pred(:, index) + output%gain(:, index) * innovation
      output%v_filt(:, :, index) = output%v_pred(:, :, index) - &
        spread(output%gain(:, index), dim=2, ncopies=2) * &
        spread(cv, dim=1, ncopies=2)
      call symmetrize(output%v_filt(:, :, index))
    end subroutine measurement_update

  end subroutine uniss_kalman

  subroutine uniss_em_update(log_data, current, fixed, updated, status, message)
    real(dp), intent(in) :: log_data(:, :)
    type(volume_parameters), intent(in) :: current
    type(parameter_mask), intent(in) :: fixed
    type(volume_parameters), intent(out) :: updated
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    type(kalman_output) :: kf
    real(dp), allocatable :: y(:), pt(:, :, :), ptt1(:, :, :)
    real(dp) :: denom, numer, cpc, residual_mean
    integer :: n_bin, n_day, n, t, b, d, idx
    real(dp), parameter :: variance_floor = 1.0e-12_dp

    if (present(status)) status = intraday_ok
    if (present(message)) message = ''
    updated = current

    call uniss_kalman(log_data, current, kf, smooth=.true.)
    if (kf%status /= intraday_ok) then
      if (present(status)) status = kf%status
      if (present(message)) message = trim(kf%message)
      return
    end if

    n_bin = size(log_data, 1)
    n_day = size(log_data, 2)
    n = n_bin * n_day
    allocate(y(n), pt(2, 2, n), ptt1(2, 2, n))
    y = reshape(log_data, [n])
    ptt1 = 0.0_dp

    do t = 1, n
      pt(:, :, t) = kf%v_smooth(:, :, t) + outer2(kf%x_smooth(:, t), kf%x_smooth(:, t))
    end do
    do t = 2, n
      ptt1(:, :, t) = matmul(kf%v_smooth(:, :, t), &
                              transpose(kf%smoother_gain(:, :, t - 1))) + &
                        outer2(kf%x_smooth(:, t), kf%x_smooth(:, t - 1))
    end do

    if (.not. fixed%x0) updated%x0 = kf%x_smooth(:, 1)
    if (.not. fixed%v0) updated%v0 = kf%v_smooth(:, :, 1)

    if (.not. fixed%phi) then
      if (.not. allocated(updated%phi)) allocate(updated%phi(n_bin))
      do b = 1, n_bin
        residual_mean = 0.0_dp
        do d = 1, n_day
          idx = b + (d - 1) * n_bin
          residual_mean = residual_mean + y(idx) - sum(kf%x_smooth(:, idx))
        end do
        updated%phi(b) = residual_mean / real(n_day, dp)
      end do
      updated%phi = updated%phi - sum(updated%phi) / real(n_bin, dp)
    end if

    if (.not. fixed%r) then
      numer = 0.0_dp
      do t = 1, n
        b = mod(t - 1, n_bin) + 1
        cpc = sum(pt(:, :, t))
        numer = numer + y(t)**2 + cpc - 2.0_dp * y(t) * sum(kf%x_smooth(:, t)) + &
                updated%phi(b)**2 - 2.0_dp * y(t) * updated%phi(b) + &
                2.0_dp * updated%phi(b) * sum(kf%x_smooth(:, t))
      end do
      updated%r = max(variance_floor, numer / real(n, dp))
    end if

    if (.not. fixed%a_eta) then
      numer = 0.0_dp
      denom = 0.0_dp
      do t = n_bin + 1, n, n_bin
        numer = numer + ptt1(1, 1, t)
        denom = denom + pt(1, 1, t - 1)
      end do
      if (abs(denom) > tiny(1.0_dp)) updated%a_eta = numer / denom
    end if

    if (.not. fixed%a_mu) then
      numer = sum(ptt1(2, 2, 2:n))
      denom = sum(pt(2, 2, 1:n - 1))
      if (abs(denom) > tiny(1.0_dp)) updated%a_mu = numer / denom
    end if

    if (.not. fixed%var_eta) then
      numer = 0.0_dp
      d = 0
      do t = n_bin + 1, n, n_bin
        numer = numer + pt(1, 1, t) + updated%a_eta**2 * pt(1, 1, t - 1) - &
                2.0_dp * updated%a_eta * ptt1(1, 1, t)
        d = d + 1
      end do
      if (d > 0) updated%var_eta = max(variance_floor, numer / real(d, dp))
    end if

    if (.not. fixed%var_mu) then
      numer = 0.0_dp
      do t = 2, n
        numer = numer + pt(2, 2, t) + updated%a_mu**2 * pt(2, 2, t - 1) - &
                2.0_dp * updated%a_mu * ptt1(2, 2, t)
      end do
      updated%var_mu = max(variance_floor, numer / real(max(1, n - 1), dp))
    end if

    call symmetrize(updated%v0)
    updated%v0(1, 1) = max(updated%v0(1, 1), variance_floor)
    updated%v0(2, 2) = max(updated%v0(2, 2), variance_floor)

    if (.not. parameters_valid(updated, n_bin)) then
      if (present(status)) status = intraday_numerical_failure
      if (present(message)) message = 'EM update produced invalid parameters'
    end if
  end subroutine uniss_em_update

  pure function outer2(a, b) result(c)
    real(dp), intent(in) :: a(2), b(2)
    real(dp) :: c(2, 2)
    c = spread(a, dim=2, ncopies=2) * spread(b, dim=1, ncopies=2)
  end function outer2

  subroutine inverse2_regularized(a, ainv)
    real(dp), intent(in) :: a(2, 2)
    real(dp), intent(out) :: ainv(2, 2)
    real(dp) :: work(2, 2), det, ridge, scale
    integer :: attempt

    work = a
    scale = max(1.0_dp, maxval(abs(a)))
    ridge = 0.0_dp
    do attempt = 1, 8
      det = work(1, 1) * work(2, 2) - work(1, 2) * work(2, 1)
      if (abs(det) > 100.0_dp * epsilon(1.0_dp) * scale * scale) exit
      ridge = max(1.0e-14_dp * scale, merge(10.0_dp * ridge, 1.0e-14_dp * scale, ridge > 0.0_dp))
      work = a
      work(1, 1) = work(1, 1) + ridge
      work(2, 2) = work(2, 2) + ridge
    end do
    det = work(1, 1) * work(2, 2) - work(1, 2) * work(2, 1)
    if (abs(det) <= tiny(1.0_dp)) then
      ainv = 0.0_dp
      ainv(1, 1) = 1.0_dp / max(abs(work(1, 1)), sqrt(tiny(1.0_dp)))
      ainv(2, 2) = 1.0_dp / max(abs(work(2, 2)), sqrt(tiny(1.0_dp)))
    else
      ainv(1, 1) = work(2, 2) / det
      ainv(1, 2) = -work(1, 2) / det
      ainv(2, 1) = -work(2, 1) / det
      ainv(2, 2) = work(1, 1) / det
    end if
  end subroutine inverse2_regularized

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(2, 2)
    real(dp) :: off
    off = 0.5_dp * (a(1, 2) + a(2, 1))
    a(1, 2) = off
    a(2, 1) = off
  end subroutine symmetrize

end module intraday_kalman
