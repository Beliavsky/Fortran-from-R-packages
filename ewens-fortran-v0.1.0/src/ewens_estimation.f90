! SPDX-License-Identifier: MIT
module ewens_estimation
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
  use ewens_kinds, only : dp
  use ewens_math, only : number_of_classes
  implicit none
  private
  public :: ewens_mle, ewens_mle_nk, ewens_score

contains

  real(dp) function ewens_score(theta, n, k) result(score)
    real(dp), intent(in) :: theta
    integer, intent(in) :: n, k
    integer :: j

    if (theta <= 0.0_dp) then
      score = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    score = real(k - 1, dp) / theta
    do j = 1, n - 1
      score = score - 1.0_dp / (theta + real(j, dp))
    end do
  end function ewens_score

  real(dp) function ewens_mle(labels) result(theta_hat)
    integer, intent(in) :: labels(:)
    integer :: n, k

    n = size(labels)
    k = number_of_classes(labels)
    theta_hat = ewens_mle_nk(n, k)
  end function ewens_mle

  real(dp) function ewens_mle_nk(n, k) result(theta_hat)
    integer, intent(in) :: n, k
    real(dp) :: lo, hi, mid, flo, fhi, fmid
    integer :: iter

    if (n <= 1 .or. k < 1 .or. k > n) then
      theta_hat = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (k == 1) then
      theta_hat = 0.0_dp
      return
    end if
    if (k == n) then
      theta_hat = ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if

    lo = max(1.0e-12_dp, epsilon(1.0_dp))
    hi = max(real(n, dp), 1.0_dp)
    flo = ewens_score(lo, n, k)
    fhi = ewens_score(hi, n, k)

    do while (fhi > 0.0_dp .and. hi < 1.0e16_dp)
      hi = 2.0_dp * hi
      fhi = ewens_score(hi, n, k)
    end do
    if (flo <= 0.0_dp .or. fhi >= 0.0_dp) then
      theta_hat = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    do iter = 1, 200
      mid = 0.5_dp * (lo + hi)
      fmid = ewens_score(mid, n, k)
      if (abs(fmid) <= 1.0e-13_dp * max(1.0_dp, abs(real(k, dp) / mid))) exit
      if (fmid > 0.0_dp) then
        lo = mid
        flo = fmid
      else
        hi = mid
        fhi = fmid
      end if
      if (hi - lo <= 1.0e-13_dp * max(1.0_dp, mid)) exit
    end do
    theta_hat = 0.5_dp * (lo + hi)
  end function ewens_mle_nk

end module ewens_estimation
