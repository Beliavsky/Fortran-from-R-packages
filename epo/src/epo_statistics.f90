! SPDX-License-Identifier: MIT
! Copyright (c) 2023 Bernardo Reckziegel
module epo_statistics
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use epo_kinds, only : dp
  implicit none
  private

  public :: sample_covariance
  public :: covariance_to_correlation
  public :: all_finite_vector
  public :: all_finite_matrix

contains

  subroutine sample_covariance(x, means, covariance, ok)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: means(:)
    real(dp), intent(out) :: covariance(:,:)
    logical, intent(out) :: ok

    real(dp), allocatable :: centered(:)
    integer :: i, j, k, n_obs, n_assets

    n_obs = size(x,1)
    n_assets = size(x,2)
    ok = .false.
    means = 0.0_dp
    covariance = 0.0_dp

    if (n_obs < 2) return
    if (size(means) /= n_assets) return
    if (size(covariance,1) /= n_assets .or. size(covariance,2) /= n_assets) return
    if (.not. all_finite_matrix(x)) return

    do j = 1, n_assets
      means(j) = sum(x(:,j)) / real(n_obs, dp)
    end do

    allocate(centered(n_assets))
    do i = 1, n_obs
      centered = x(i,:) - means
      do j = 1, n_assets
        do k = 1, j
          covariance(j,k) = covariance(j,k) + centered(j) * centered(k)
        end do
      end do
    end do

    covariance = covariance / real(n_obs - 1, dp)
    do j = 1, n_assets
      do k = 1, j - 1
        covariance(k,j) = covariance(j,k)
      end do
    end do

    ok = .true.
  end subroutine sample_covariance

  subroutine covariance_to_correlation(covariance, correlation, ok)
    real(dp), intent(in) :: covariance(:,:)
    real(dp), intent(out) :: correlation(:,:)
    logical, intent(out) :: ok

    real(dp), allocatable :: sd(:)
    integer :: i, j, n

    n = size(covariance,1)
    ok = .false.
    correlation = 0.0_dp

    if (size(covariance,2) /= n) return
    if (size(correlation,1) /= n .or. size(correlation,2) /= n) return
    if (.not. all_finite_matrix(covariance)) return

    allocate(sd(n))
    do i = 1, n
      if (covariance(i,i) <= 0.0_dp) return
      sd(i) = sqrt(covariance(i,i))
    end do

    do i = 1, n
      correlation(i,i) = 1.0_dp
      do j = 1, i - 1
        correlation(i,j) = covariance(i,j) / (sd(i) * sd(j))
        correlation(j,i) = correlation(i,j)
      end do
    end do

    ok = .true.
  end subroutine covariance_to_correlation

  pure logical function all_finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i

    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function all_finite_vector

  pure logical function all_finite_matrix(x) result(ok)
    real(dp), intent(in) :: x(:,:)
    integer :: i, j

    ok = .true.
    do j = 1, size(x,2)
      do i = 1, size(x,1)
        if (.not. ieee_is_finite(x(i,j))) then
          ok = .false.
          return
        end if
      end do
    end do
  end function all_finite_matrix

end module epo_statistics
