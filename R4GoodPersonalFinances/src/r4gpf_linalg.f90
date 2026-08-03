! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_linalg
  use r4gpf_kinds, only: dp
  use r4gpf_status, only: r4gpf_success, r4gpf_dimension_error, r4gpf_numerical_error
  implicit none
  private
  public :: cholesky_lower, covariance_from_sd_corr, solve_linear, project_simplex
  public :: quadratic_form, vector_norm2
contains

  pure real(dp) function vector_norm2(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm2

  pure real(dp) function quadratic_form(x, a) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: a(:, :)
    value = dot_product(x, matmul(a, x))
  end function quadratic_form

  subroutine covariance_from_sd_corr(sd, corr, covariance, status)
    real(dp), intent(in) :: sd(:)
    real(dp), intent(in) :: corr(:, :)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out), optional :: status
    integer :: i, j, n

    n = size(sd)
    if (size(corr, 1) /= n .or. size(corr, 2) /= n) then
      allocate(covariance(0, 0))
      if (present(status)) status = r4gpf_dimension_error
      return
    end if
    allocate(covariance(n, n))
    do j = 1, n
      do i = 1, n
        covariance(i, j) = sd(i) * corr(i, j) * sd(j)
      end do
    end do
    if (present(status)) status = r4gpf_success
  end subroutine covariance_from_sd_corr

  subroutine cholesky_lower(a, l, status, jitter)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: l(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: jitter
    real(dp) :: s, eps
    integer :: i, j, k, n

    n = size(a, 1)
    if (size(a, 2) /= n) then
      allocate(l(0, 0))
      status = r4gpf_dimension_error
      return
    end if
    eps = 0.0_dp
    if (present(jitter)) eps = max(0.0_dp, jitter)
    allocate(l(n, n))
    l = 0.0_dp
    do i = 1, n
      do j = 1, i
        s = a(i, j)
        if (i == j) s = s + eps
        do k = 1, j - 1
          s = s - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            status = r4gpf_numerical_error
            return
          end if
          l(i, j) = sqrt(s)
        else
          if (abs(l(j, j)) <= tiny(1.0_dp)) then
            status = r4gpf_numerical_error
            return
          end if
          l(i, j) = s / l(j, j)
        end if
      end do
    end do
    status = r4gpf_success
  end subroutine cholesky_lower

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:, :), rowtmp(:)
    real(dp) :: pivot, factor
    integer :: i, k, p, n

    n = size(b)
    if (size(a, 1) /= n .or. size(a, 2) /= n) then
      allocate(x(0))
      status = r4gpf_dimension_error
      return
    end if
    allocate(aug(n, n + 1), rowtmp(n + 1), x(n))
    aug(:, 1:n) = a
    aug(:, n + 1) = b
    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(aug(i, k)) > abs(aug(p, k))) p = i
      end do
      if (abs(aug(p, k)) <= 100.0_dp * epsilon(1.0_dp)) then
        status = r4gpf_numerical_error
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        rowtmp = aug(k, :)
        aug(k, :) = aug(p, :)
        aug(p, :) = rowtmp
      end if
      pivot = aug(k, k)
      aug(k, k:n + 1) = aug(k, k:n + 1) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i, k)
        if (abs(factor) <= tiny(1.0_dp)) cycle
        aug(i, k:n + 1) = aug(i, k:n + 1) - factor * aug(k, k:n + 1)
      end do
    end do
    x = aug(:, n + 1)
    status = r4gpf_success
  end subroutine solve_linear

  subroutine project_simplex(v, target_sum, projected)
    real(dp), intent(in) :: v(:)
    real(dp), intent(in) :: target_sum
    real(dp), intent(out) :: projected(:)
    real(dp), allocatable :: u(:)
    real(dp) :: cssv, theta
    integer :: i, j, n, rho

    n = size(v)
    if (size(projected) /= n) error stop "project_simplex: size mismatch"
    if (target_sum <= 0.0_dp) then
      projected = 0.0_dp
      return
    end if
    allocate(u(n))
    u = v
    do i = 2, n
      theta = u(i)
      j = i - 1
      do while (j >= 1)
        if (u(j) >= theta) exit
        u(j + 1) = u(j)
        j = j - 1
      end do
      u(j + 1) = theta
    end do
    cssv = 0.0_dp
    rho = 1
    do i = 1, n
      cssv = cssv + u(i)
      theta = (cssv - target_sum) / real(i, dp)
      if (u(i) - theta > 0.0_dp) rho = i
    end do
    theta = (sum(u(1:rho)) - target_sum) / real(rho, dp)
    projected = max(v - theta, 0.0_dp)
    if (sum(projected) > 0.0_dp) projected = projected * target_sum / sum(projected)
  end subroutine project_simplex

end module r4gpf_linalg
