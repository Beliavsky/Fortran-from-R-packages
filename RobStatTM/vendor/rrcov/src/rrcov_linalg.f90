! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_linalg
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success, rrcov_dimension_error, rrcov_singular, rrcov_no_convergence
  use rrcov_sort, only : sort_real_with_index
  implicit none
  private
  public :: identity_matrix, symmetrize, symmetric_eigen, symmetric_inverse
  public :: general_inverse, solve_linear, determinant, log_determinant
  public :: mahalanobis_squared, make_positive_definite, matrix_sqrt
  public :: orthonormalize, matrix_rank, outer_product
contains
  function identity_matrix(n) result(value)
    integer, intent(in) :: n
    real(dp) :: value(n, n)
    integer :: i
    value = 0.0_dp
    do i = 1, n
      value(i, i) = 1.0_dp
    end do
  end function identity_matrix

  pure function symmetrize(a) result(value)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: value(size(a, 1), size(a, 2))
    value = 0.5_dp * (a + transpose(a))
  end function symmetrize

  pure function outer_product(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: value(size(x), size(y))
    integer :: i
    do i = 1, size(x)
      value(i, :) = x(i) * y
    end do
  end function outer_product

  subroutine symmetric_eigen(a, values, vectors, status, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps
    real(dp), allocatable :: work(:, :), sorted(:)
    integer, allocatable :: order(:)
    real(dp) :: app, aqq, apq, tau, t, c, s, off, tol
    integer :: n, p, q, i, sweep, maxit

    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      allocate(values(0), vectors(0, 0))
      status = rrcov_dimension_error
      return
    end if
    allocate(work(n, n), values(n), vectors(n, n), sorted(n), order(n))
    work = symmetrize(a)
    vectors = identity_matrix(n)
    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = max(50, 20 * n * n)
    if (present(max_sweeps)) maxit = max(1, max_sweeps)

    do sweep = 1, maxit
      off = 0.0_dp
      do q = 2, n
        do p = 1, q - 1
          off = max(off, abs(work(p, q)))
          if (abs(work(p, q)) <= tol * max(1.0_dp, abs(work(p, p)) + abs(work(q, q)))) cycle
          app = work(p, p)
          aqq = work(q, q)
          apq = work(p, q)
          tau = (aqq - app) / (2.0_dp * apq)
          if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
          else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
          end if
          c = 1.0_dp / sqrt(1.0_dp + t * t)
          s = t * c
          work(p, p) = app - t * apq
          work(q, q) = aqq + t * apq
          work(p, q) = 0.0_dp
          work(q, p) = 0.0_dp
          do i = 1, n
            if (i /= p .and. i /= q) then
              app = work(i, p)
              aqq = work(i, q)
              work(i, p) = c * app - s * aqq
              work(p, i) = work(i, p)
              work(i, q) = s * app + c * aqq
              work(q, i) = work(i, q)
            end if
            app = vectors(i, p)
            aqq = vectors(i, q)
            vectors(i, p) = c * app - s * aqq
            vectors(i, q) = s * app + c * aqq
          end do
        end do
      end do
      if (off <= tol * max(1.0_dp, maxval(abs(work)))) exit
    end do

    values = [(work(i, i), i=1, n)]
    sorted = -values
    order = [(i, i=1, n)]
    call sort_real_with_index(sorted, order)
    values = -sorted
    vectors = vectors(:, order)
    if (sweep > maxit) then
      status = rrcov_no_convergence
    else
      status = rrcov_success
    end if
  end subroutine symmetric_eigen

  function symmetric_inverse(a, status, tolerance) result(inverse)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: inverse(:, :), values(:), vectors(:, :)
    real(dp) :: tol, vmax
    integer :: i, istat, n, rank

    n = size(a, 1)
    allocate(inverse(n, n))
    inverse = 0.0_dp
    if (n < 1 .or. size(a, 2) /= n) then
      if (present(status)) status = rrcov_dimension_error
      return
    end if
    call symmetric_eigen(a, values, vectors, istat)
    vmax = max(1.0_dp, maxval(abs(values)))
    tol = sqrt(epsilon(1.0_dp)) * vmax
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp)) * vmax
    rank = 0
    do i = 1, n
      if (abs(values(i)) > tol) then
        inverse = inverse + outer_product(vectors(:, i), vectors(:, i)) / values(i)
        rank = rank + 1
      end if
    end do
    inverse = symmetrize(inverse)
    if (present(status)) then
      if (rank < n) then
        status = rrcov_singular
      else
        status = istat
      end if
    end if
  end function symmetric_inverse

  function general_inverse(a, status, tolerance) result(inverse)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: inverse(:, :), aug(:, :), temp(:)
    real(dp) :: pivot, tol
    integer :: n, i, k, pivot_row

    n = size(a, 1)
    allocate(inverse(n, n))
    inverse = 0.0_dp
    if (n < 1 .or. size(a, 2) /= n) then
      if (present(status)) status = rrcov_dimension_error
      return
    end if
    allocate(aug(n, 2 * n), temp(2 * n))
    aug(:, 1:n) = a
    aug(:, n + 1:2 * n) = identity_matrix(n)
    tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(a)))
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(a)))
    do k = 1, n
      pivot_row = k - 1 + maxloc(abs(aug(k:n, k)), dim=1)
      pivot = aug(pivot_row, k)
      if (abs(pivot) <= tol) then
        if (present(status)) status = rrcov_singular
        return
      end if
      if (pivot_row /= k) then
        temp = aug(k, :)
        aug(k, :) = aug(pivot_row, :)
        aug(pivot_row, :) = temp
      end if
      aug(k, :) = aug(k, :) / aug(k, k)
      do i = 1, n
        if (i == k) cycle
        pivot = aug(i, k)
        aug(i, :) = aug(i, :) - pivot * aug(k, :)
      end do
    end do
    inverse = aug(:, n + 1:2 * n)
    if (present(status)) status = rrcov_success
  end function general_inverse

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: inverse(:, :)
    if (size(a, 1) /= size(a, 2) .or. size(a, 1) /= size(b)) then
      allocate(x(0))
      status = rrcov_dimension_error
      return
    end if
    inverse = general_inverse(a, status)
    allocate(x(size(b)))
    if (status == rrcov_success) then
      x = matmul(inverse, b)
    else
      x = 0.0_dp
    end if
  end subroutine solve_linear

  function determinant(a, status) result(value)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp) :: value
    real(dp), allocatable :: work(:, :), temp(:)
    real(dp) :: pivot
    integer :: n, i, k, pivot_row, sign
    n = size(a, 1)
    if (n < 1 .or. size(a, 2) /= n) then
      value = 0.0_dp
      if (present(status)) status = rrcov_dimension_error
      return
    end if
    allocate(work(n, n), temp(n))
    work = a
    value = 1.0_dp
    sign = 1
    do k = 1, n
      pivot_row = k - 1 + maxloc(abs(work(k:n, k)), dim=1)
      pivot = work(pivot_row, k)
      if (abs(pivot) <= tiny(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        value = 0.0_dp
        if (present(status)) status = rrcov_singular
        return
      end if
      if (pivot_row /= k) then
        temp = work(k, :)
        work(k, :) = work(pivot_row, :)
        work(pivot_row, :) = temp
        sign = -sign
      end if
      value = value * work(k, k)
      do i = k + 1, n
        work(i, k + 1:n) = work(i, k + 1:n) - work(i, k) / work(k, k) * work(k, k + 1:n)
      end do
    end do
    value = real(sign, dp) * value
    if (present(status)) status = rrcov_success
  end function determinant

  function log_determinant(a, status) result(value)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp) :: value
    real(dp), allocatable :: values(:), vectors(:, :)
    integer :: i, istat
    call symmetric_eigen(a, values, vectors, istat)
    value = 0.0_dp
    do i = 1, size(values)
      if (values(i) <= tiny(1.0_dp)) then
        value = -huge(1.0_dp)
        if (present(status)) status = rrcov_singular
        return
      end if
      value = value + log(values(i))
    end do
    if (present(status)) status = istat
  end function log_determinant

  subroutine mahalanobis_squared(x, center, covariance, distances, status)
    real(dp), intent(in) :: x(:, :), center(:), covariance(:, :)
    real(dp), allocatable, intent(out) :: distances(:)
    integer, intent(out) :: status
    real(dp), allocatable :: inverse(:, :), delta(:)
    integer :: i, istat, n, p
    n = size(x, 1)
    p = size(x, 2)
    allocate(distances(n), delta(p))
    if (size(center) /= p .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) then
      distances = huge(1.0_dp)
      status = rrcov_dimension_error
      return
    end if
    inverse = symmetric_inverse(covariance, istat)
    do i = 1, n
      delta = x(i, :) - center
      distances(i) = max(0.0_dp, dot_product(delta, matmul(inverse, delta)))
    end do
    status = istat
    if (istat == rrcov_singular) status = rrcov_success
  end subroutine mahalanobis_squared

  function make_positive_definite(a, floor_ratio, status) result(value)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: floor_ratio
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :), values(:), vectors(:, :)
    real(dp) :: floor_value, ratio, vmax
    integer :: i, istat, n
    n = size(a, 1)
    allocate(value(n, n))
    call symmetric_eigen(a, values, vectors, istat)
    ratio = 1.0e-8_dp
    if (present(floor_ratio)) ratio = max(floor_ratio, epsilon(1.0_dp))
    vmax = max(1.0_dp, maxval(abs(values)))
    floor_value = ratio * vmax
    value = 0.0_dp
    do i = 1, n
      value = value + max(values(i), floor_value) * outer_product(vectors(:, i), vectors(:, i))
    end do
    value = symmetrize(value)
    if (present(status)) status = istat
  end function make_positive_definite

  function matrix_sqrt(a, inverse, status) result(value)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in), optional :: inverse
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :), values(:), vectors(:, :)
    logical :: inv
    real(dp) :: v, floor_value
    integer :: i, istat, n
    n = size(a, 1)
    call symmetric_eigen(a, values, vectors, istat)
    allocate(value(n, n))
    value = 0.0_dp
    inv = .false.
    if (present(inverse)) inv = inverse
    floor_value = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(values)))
    do i = 1, n
      v = max(values(i), floor_value)
      if (inv) then
        v = 1.0_dp / sqrt(v)
      else
        v = sqrt(v)
      end if
      value = value + v * outer_product(vectors(:, i), vectors(:, i))
    end do
    value = symmetrize(value)
    if (present(status)) status = istat
  end function matrix_sqrt

  subroutine orthonormalize(vectors, rank, tolerance)
    real(dp), intent(inout) :: vectors(:, :)
    integer, intent(out), optional :: rank
    real(dp), intent(in), optional :: tolerance
    real(dp) :: norm, tol
    integer :: i, j, r
    tol = sqrt(epsilon(1.0_dp))
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    r = 0
    do j = 1, size(vectors, 2)
      do i = 1, j - 1
        vectors(:, j) = vectors(:, j) - dot_product(vectors(:, i), vectors(:, j)) * vectors(:, i)
      end do
      norm = sqrt(sum(vectors(:, j) ** 2))
      if (norm > tol) then
        vectors(:, j) = vectors(:, j) / norm
        r = r + 1
      else
        vectors(:, j) = 0.0_dp
      end if
    end do
    if (present(rank)) rank = r
  end subroutine orthonormalize

  function matrix_rank(a, tolerance) result(rank)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    integer :: rank
    real(dp), allocatable :: gram(:, :), values(:), vectors(:, :)
    real(dp) :: tol, vmax
    integer :: istat
    gram = matmul(transpose(a), a)
    call symmetric_eigen(gram, values, vectors, istat)
    vmax = max(1.0_dp, maxval(abs(values)))
    tol = sqrt(epsilon(1.0_dp)) * vmax
    if (present(tolerance)) tol = tolerance * vmax
    rank = count(values > tol)
  end function matrix_rank
end module rrcov_linalg
