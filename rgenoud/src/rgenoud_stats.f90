! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_stats
  use rgenoud_kinds, only : dp
  implicit none
  private
  public :: sample_moments
contains
  subroutine sample_moments(x, mean, variance, skewness, kurtosis, weights)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: mean(:), variance(:), skewness(:), kurtosis(:)
    real(dp), intent(in), optional :: weights(:)
    real(dp) :: w(size(x, 1)), sw, d, m2, m3, m4
    integer :: i, j

    if (size(mean) /= size(x, 2) .or. size(variance) /= size(x, 2) .or. &
        size(skewness) /= size(x, 2) .or. size(kurtosis) /= size(x, 2)) then
      error stop "sample_moments output size mismatch"
    end if
    if (present(weights)) then
      if (size(weights) /= size(x, 1)) error stop "sample_moments weight size mismatch"
      w = weights
    else
      w = 1.0_dp
    end if
    sw = sum(w)
    if (sw <= 0.0_dp) error stop "sample_moments requires positive total weight"

    do j = 1, size(x, 2)
      mean(j) = dot_product(w, x(:, j)) / sw
      m2 = 0.0_dp
      m3 = 0.0_dp
      m4 = 0.0_dp
      do i = 1, size(x, 1)
        d = x(i, j) - mean(j)
        m2 = m2 + w(i) * d**2
        m3 = m3 + w(i) * d**3
        m4 = m4 + w(i) * d**4
      end do
      variance(j) = m2 / sw
      if (variance(j) > 0.0_dp) then
        skewness(j) = (m3 / sw) / variance(j)**1.5_dp
        kurtosis(j) = (m4 / sw) / variance(j)**2
      else
        skewness(j) = 0.0_dp
        kurtosis(j) = 0.0_dp
      end if
    end do
  end subroutine sample_moments
end module rgenoud_stats
