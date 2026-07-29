! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_statistics
  use frapo_kinds, only : dp
  use frapo_types, only : frapo_ok, frapo_invalid_input
  implicit none
  private
  public :: sample_covariance, covariance_to_correlation, average_ranks, is_symmetric

contains

  pure logical function is_symmetric(a, tolerance) result(ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol
    integer :: i, j

    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance
    if (size(a, 1) /= size(a, 2)) then
      ok = .false.
      return
    end if
    ok = .true.
    do j = 1, size(a, 2)
      do i = j + 1, size(a, 1)
        if (abs(a(i, j) - a(j, i)) > tol * max(1.0_dp, abs(a(i, j)), abs(a(j, i)))) then
          ok = .false.
          return
        end if
      end do
    end do
  end function is_symmetric

  subroutine sample_covariance(x, covariance, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: centered(:, :), means(:)
    integer :: nobs, nvar, j

    nobs = size(x, 1)
    nvar = size(x, 2)
    allocate(covariance(nvar, nvar))
    if (nobs < 2 .or. nvar < 1) then
      covariance = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(centered(nobs, nvar), means(nvar))
    means = sum(x, dim=1) / real(nobs, dp)
    do j = 1, nvar
      centered(:, j) = x(:, j) - means(j)
    end do
    covariance = matmul(transpose(centered), centered) / real(nobs - 1, dp)
    if (present(status)) status = frapo_ok
  end subroutine sample_covariance

  subroutine covariance_to_correlation(covariance, correlation, status)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: correlation(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: sd(:)
    integer :: i, j, n

    n = size(covariance, 1)
    allocate(correlation(n, n), sd(n))
    if (size(covariance, 2) /= n .or. any([(covariance(i, i) <= 0.0_dp, i=1,n)])) then
      correlation = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    do i = 1, n
      sd(i) = sqrt(covariance(i, i))
    end do
    do j = 1, n
      do i = 1, n
        correlation(i, j) = covariance(i, j) / (sd(i) * sd(j))
      end do
    end do
    if (present(status)) status = frapo_ok
  end subroutine covariance_to_correlation

  subroutine average_ranks(x, ranks)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: ranks(:)
    integer :: i, j, less_count, equal_count, n

    n = size(x)
    allocate(ranks(n))
    do i = 1, n
      less_count = 0
      equal_count = 0
      do j = 1, n
        if (x(j) < x(i)) then
          less_count = less_count + 1
        else if (.not. (x(j) < x(i)) .and. .not. (x(j) > x(i))) then
          equal_count = equal_count + 1
        end if
      end do
      ranks(i) = real(less_count, dp) + 0.5_dp * real(equal_count + 1, dp)
    end do
  end subroutine average_ranks
end module frapo_statistics
