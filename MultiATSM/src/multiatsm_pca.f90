! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_pca
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : covariance_rows, symmetric_eigen
  implicit none
  private

  public :: pca_weights_one_country, spanned_factors, pca_variance_explained

contains

  subroutine pca_weights_one_country(yields, weights, eigenvalues, info)
    real(dp), intent(in) :: yields(:, :)
    real(dp), allocatable, intent(out) :: weights(:, :)
    real(dp), allocatable, intent(out) :: eigenvalues(:)
    integer, intent(out) :: info
    real(dp), allocatable :: values_asc(:), vectors_asc(:, :)
    real(dp) :: cov(size(yields, 1), size(yields, 1))
    integer :: j, i, mid

    j = size(yields, 1)
    cov = covariance_rows(yields)
    call symmetric_eigen(cov, values_asc, vectors_asc, info)
    if (info /= 0) then
      allocate(weights(0, 0), eigenvalues(0))
      return
    end if

    allocate(weights(j, j), eigenvalues(j))
    do i = 1, j
      eigenvalues(i) = values_asc(j - i + 1)
      weights(i, :) = vectors_asc(:, j - i + 1)
    end do

    if (j >= 1) then
      if (all(weights(1, :) < 0.0_dp)) weights(1, :) = -weights(1, :)
    end if
    if (j >= 2) then
      if (weights(2, 1) > weights(2, j)) weights(2, :) = -weights(2, :)
    end if
    if (j >= 3) then
      mid = nint(0.5_dp * real(j + 1, dp))
      if (weights(3, 1) > weights(3, mid) .and. weights(3, j) > weights(3, mid)) then
        weights(3, :) = -weights(3, :)
      end if
    end if
  end subroutine pca_weights_one_country

  subroutine spanned_factors(yields, n_countries, n_factors, factors, weights, info, scale_percent)
    real(dp), intent(in) :: yields(:, :)
    integer, intent(in) :: n_countries, n_factors
    real(dp), allocatable, intent(out) :: factors(:, :)
    real(dp), allocatable, intent(out) :: weights(:, :)
    integer, intent(out) :: info
    logical, intent(in), optional :: scale_percent
    real(dp), allocatable :: w_country(:, :), eig(:)
    real(dp) :: scale
    integer :: j, c, r0, r1, f0, f1

    info = 0
    if (n_countries <= 0 .or. n_factors <= 0) then
      info = -1
      allocate(factors(0, 0), weights(0, 0))
      return
    end if
    if (mod(size(yields, 1), n_countries) /= 0) then
      info = -2
      allocate(factors(0, 0), weights(0, 0))
      return
    end if
    j = size(yields, 1) / n_countries
    if (n_factors > j) then
      info = -3
      allocate(factors(0, 0), weights(0, 0))
      return
    end if

    scale = 100.0_dp
    if (present(scale_percent)) then
      if (.not. scale_percent) scale = 1.0_dp
    end if
    allocate(factors(n_countries * n_factors, size(yields, 2)))
    allocate(weights(n_countries * n_factors, n_countries * j))
    weights = 0.0_dp

    do c = 1, n_countries
      r0 = (c - 1) * j + 1
      r1 = c * j
      f0 = (c - 1) * n_factors + 1
      f1 = c * n_factors
      call pca_weights_one_country(yields(r0:r1, :), w_country, eig, info)
      if (info /= 0) return
      weights(f0:f1, r0:r1) = scale * w_country(1:n_factors, :)
      factors(f0:f1, :) = matmul(weights(f0:f1, r0:r1), yields(r0:r1, :))
    end do
  end subroutine spanned_factors

  subroutine pca_variance_explained(yields, n_countries, n_factors, explained, info)
    real(dp), intent(in) :: yields(:, :)
    integer, intent(in) :: n_countries, n_factors
    real(dp), allocatable, intent(out) :: explained(:)
    integer, intent(out) :: info
    real(dp), allocatable :: weights(:, :), eig(:)
    integer :: j, c, r0, r1, f0, f1

    info = 0
    if (n_countries <= 0 .or. mod(size(yields, 1), n_countries) /= 0) then
      info = -1
      allocate(explained(0))
      return
    end if
    j = size(yields, 1) / n_countries
    if (n_factors < 1 .or. n_factors > j) then
      info = -2
      allocate(explained(0))
      return
    end if
    allocate(explained(n_countries * n_factors))
    do c = 1, n_countries
      r0 = (c - 1) * j + 1
      r1 = c * j
      f0 = (c - 1) * n_factors + 1
      f1 = c * n_factors
      call pca_weights_one_country(yields(r0:r1, :), weights, eig, info)
      if (info /= 0) return
      if (sum(max(eig, 0.0_dp)) > 0.0_dp) then
        explained(f0:f1) = eig(1:n_factors) / sum(max(eig, 0.0_dp))
      else
        explained(f0:f1) = 0.0_dp
      end if
    end do
  end subroutine pca_variance_explained

end module multiatsm_pca
