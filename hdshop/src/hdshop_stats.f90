! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: symmetrize
  implicit none
  private
  public :: row_means, sample_covariance, normal_cdf, normal_quantile
  public :: chi_square1_survival

contains

  function row_means(x) result(means)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: means(:)
    integer :: i, j, count
    allocate(means(size(x,1)))
    means = 0.0_dp
    do i = 1, size(x,1)
      count = 0
      do j = 1, size(x,2)
        if (ieee_is_finite(x(i,j))) then
          means(i) = means(i) + x(i,j)
          count = count + 1
        end if
      end do
      if (count > 0) means(i) = means(i)/real(count,dp)
    end do
  end function row_means

  function sample_covariance(x, ok) result(covariance)
    real(dp), intent(in) :: x(:,:)
    logical, intent(out), optional :: ok
    real(dp), allocatable :: covariance(:,:), means(:)
    integer :: p, n, i, j, k, count
    p = size(x,1); n = size(x,2)
    allocate(covariance(p,p))
    covariance = 0.0_dp
    means = row_means(x)
    do j = 1, p
      do i = 1, j
        count = 0
        do k = 1, n
          if (ieee_is_finite(x(i,k)) .and. ieee_is_finite(x(j,k))) then
            covariance(i,j) = covariance(i,j) + &
              (x(i,k)-means(i))*(x(j,k)-means(j))
            count = count + 1
          end if
        end do
        if (count > 1) covariance(i,j) = covariance(i,j)/real(count-1,dp)
        covariance(j,i) = covariance(i,j)
      end do
    end do
    call symmetrize(covariance)
    if (present(ok)) ok = n > 1
  end function sample_covariance

  elemental real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  real(dp) function normal_quantile(probability) result(x)
    real(dp), intent(in) :: probability
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r, e, u
    if (probability <= 0.0_dp) then
      x = -huge(1.0_dp); return
    else if (probability >= 1.0_dp) then
      x = huge(1.0_dp); return
    else if (probability < plow) then
      q = sqrt(-2.0_dp*log(probability))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (probability <= phigh) then
      q = probability-0.5_dp; r=q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-probability))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    ! One Halley refinement.
    e = normal_cdf(x)-probability
    u = e*sqrt(2.0_dp*acos(-1.0_dp))*exp(0.5_dp*x*x)
    x = x-u/(1.0_dp+0.5_dp*x*u)
  end function normal_quantile

  elemental real(dp) function chi_square1_survival(x) result(value)
    real(dp), intent(in) :: x
    if (x <= 0.0_dp) then
      value = 1.0_dp
    else
      value = erfc(sqrt(0.5_dp*x))
    end if
  end function chi_square1_survival

end module hdshop_stats
