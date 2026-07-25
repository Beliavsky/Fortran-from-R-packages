! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_linalg
  use chaos_kinds, only : dp
  implicit none
  private
  public :: ols_fit, pca_scores_standardized

  interface
    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      import dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      real(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine dgels

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*), work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface
contains
  subroutine ols_fit(x, y, coefficients, residuals, ssr, status)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), allocatable, intent(out) :: coefficients(:), residuals(:)
    real(dp), intent(out) :: ssr
    integer, intent(out) :: status
    real(dp), allocatable :: a(:, :), b(:, :), work(:)
    real(dp) :: work_query(1)
    integer :: nobs, npar, ldb, lwork, info

    nobs = size(x, 1)
    npar = size(x, 2)
    if (size(y) /= nobs .or. nobs < npar .or. npar < 1) then
      allocate(coefficients(0), residuals(0))
      ssr = huge(1.0_dp)
      status = 1
      return
    end if
    ldb = max(nobs, npar)
    allocate(a(nobs, npar), b(ldb, 1))
    a = x
    b = 0.0_dp
    b(1:nobs, 1) = y
    call dgels('N', nobs, npar, 1, a, nobs, b, ldb, work_query, -1, info)
    if (info /= 0) then
      allocate(coefficients(0), residuals(0))
      ssr = huge(1.0_dp)
      status = 2
      return
    end if
    lwork = max(1, ceiling(work_query(1)))
    allocate(work(lwork))
    a = x
    b = 0.0_dp
    b(1:nobs, 1) = y
    call dgels('N', nobs, npar, 1, a, nobs, b, ldb, work, lwork, info)
    if (info /= 0) then
      allocate(coefficients(0), residuals(0))
      ssr = huge(1.0_dp)
      status = 3
      return
    end if
    allocate(coefficients(npar), residuals(nobs))
    coefficients = b(1:npar, 1)
    residuals = y - matmul(x, coefficients)
    ssr = sum(residuals**2)
    status = 0
  end subroutine ols_fit

  subroutine pca_scores_standardized(data, scores, eigenvalues, status)
    real(dp), intent(in) :: data(:, :)
    real(dp), allocatable, intent(out) :: scores(:, :), eigenvalues(:)
    integer, intent(out) :: status
    real(dp), allocatable :: z(:, :), covariance(:, :), work(:), eval_ascending(:)
    real(dp) :: work_query(1), mean_value, sd_value
    integer :: nobs, nvar, i, j, lwork, info

    nobs = size(data, 1)
    nvar = size(data, 2)
    if (nobs < 2 .or. nvar < 1) then
      allocate(scores(0, 0), eigenvalues(0))
      status = 1
      return
    end if
    allocate(z(nobs, nvar))
    do j = 1, nvar
      mean_value = sum(data(:, j)) / real(nobs, dp)
      sd_value = sqrt(sum((data(:, j) - mean_value)**2) / real(nobs - 1, dp))
      if (sd_value <= sqrt(tiny(1.0_dp))) then
        allocate(scores(0, 0), eigenvalues(0))
        status = 2
        return
      end if
      z(:, j) = (data(:, j) - mean_value) / sd_value
    end do
    allocate(covariance(nvar, nvar), eval_ascending(nvar))
    covariance = matmul(transpose(z), z) / real(nobs - 1, dp)
    call dsyev('V', 'U', nvar, covariance, nvar, eval_ascending, work_query, -1, info)
    if (info /= 0) then
      allocate(scores(0, 0), eigenvalues(0))
      status = 3
      return
    end if
    lwork = max(1, ceiling(work_query(1)))
    allocate(work(lwork))
    covariance = matmul(transpose(z), z) / real(nobs - 1, dp)
    call dsyev('V', 'U', nvar, covariance, nvar, eval_ascending, work, lwork, info)
    if (info /= 0) then
      allocate(scores(0, 0), eigenvalues(0))
      status = 4
      return
    end if
    allocate(scores(nobs, nvar), eigenvalues(nvar))
    do i = 1, nvar
      eigenvalues(i) = eval_ascending(nvar - i + 1)
      scores(:, i) = matmul(z, covariance(:, nvar - i + 1))
      call normalize_sign(scores(:, i), covariance(:, nvar - i + 1))
    end do
    status = 0
  contains
    subroutine normalize_sign(score, loading)
      real(dp), intent(inout) :: score(:), loading(:)
      integer :: idx(1)
      idx = maxloc(abs(loading))
      if (loading(idx(1)) < 0.0_dp) score = -score
    end subroutine normalize_sign
  end subroutine pca_scores_standardized
end module fnonlinear_linalg
